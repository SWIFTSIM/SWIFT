#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location=${scripts_location:=../../GEAR_ICs_and_SCRIPTS}

# Make run.sh fail if a subcommand fails; swift's exit code is piped into tee.
set -eo pipefail

n_threads=${n_threads:=8}  #Number of threads to use
redshift=${redshift:=9} #Starting redshift
cmb_floor=${cmb_floor:=1} #1: CMB temperature floor of the cooling library on, 0: off
gas_density=${gas_density:=10} #Physical gas density in atom/cm^3
gas_mass=${gas_mass:=1e4} #Mass of the gas particles (Msun)
temperature=${temperature:=100} #Initial gas temperature (K)
expansion=${expansion:=1.25} #Ratio a_end / a_begin
n_steps=${n_steps:=400} #Requested number of time steps; SWIFT rounds it up to a power of two
n_snapshots=${n_snapshots:=8} #Number of snapshots after the first; keep it a power of two so all snapshots fall at the same point of a time step
swift_exe=${swift_exe:="$PWD/../../../swift"} #swift binary
run_dir=${run_dir:=run_z${redshift}_floor${cmb_floor}} #Isolated directory of the run
extra_args=${extra_args:=""} #Further swift arguments, e.g. "-P Snapshots:compression:4"

glass_file="glassCube_32.hdf5"
yields_table="POPIIsw.h5"

if [ ! -e "$glass_file" ]
then
    echo "Fetching initial glass file (${glass_file})..."
    ./getGlass.sh 32
fi

if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    "$scripts_location"/getGrackleCoolingTable.sh
fi

if [ ! -e "$yields_table" ]
then
    echo "Fetching the stellar yields table..."
    "$scripts_location"/getChemistryTable.sh
fi

here=$PWD
# The run directory is emptied below: accept only a plain name of a directory under this one.
case "$run_dir" in
    ""|.|..|*/*)
	echo "run_dir must be a plain directory name, not '$run_dir'." >&2
	exit 1
	;;
esac
if [ -d "$run_dir" ]
then
    echo "$run_dir exists. Its content will be removed."
    rm -r "$run_dir"
fi
mkdir -p "$run_dir/snap"
cd "$run_dir"
ln -s "$here/$glass_file" "$here/CloudyData_UVB=HM2012.h5" "$here/$yields_table" .
cp "$here/params.yml" .

echo "Generating initial conditions to run the example..."
python3 "$here/makeIC.py" --glass "$glass_file" --rho $gas_density --mass $gas_mass \
	--temperature $temperature --redshift $redshift -o ICs_cmb_floor.hdf5

# The time line: a_end, dt_max (a step in ln a) and delta_time (a ratio of scale factors).
eval "$(python3 -c "
import math
a = 1.0 / (1.0 + $redshift)
print(f'a_begin={a:.17g}')
print(f'a_end={a * $expansion:.17g}')
print(f'dt_max={math.log($expansion) / $n_steps:.17g}')
print(f'delta_time={$expansion ** (1.0 / $n_snapshots):.17g}')
")"

echo "Running simulation..."

"$swift_exe" --hydro --cooling --cosmology --sync --limiter --verbose=0 --threads=$n_threads \
	     -P Cosmology:a_begin:$a_begin -P Cosmology:a_end:$a_end \
	     -P Snapshots:scale_factor_first:$a_begin -P Statistics:scale_factor_first:$a_begin \
	     -P TimeIntegration:dt_max:$dt_max -P Snapshots:delta_time:$delta_time \
	     -P Statistics:delta_time:$delta_time \
	     -P GrackleCooling:cmb_temperature_floor:$cmb_floor \
	     -P GeajayCooling:metal_cooling_cmb_floor:$cmb_floor \
	     $extra_args params.yml 2>&1 | tee output.log
