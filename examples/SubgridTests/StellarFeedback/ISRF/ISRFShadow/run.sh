#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

# Make run.sh fail if a subcommand fails; swift's exit code is piped into tee.
set -eo pipefail

n_threads=${n_threads:=4}  #Number of threads to use
gas_density=${gas_density:=1e3} #Gas density in atom/cm^3
gas_mass=${gas_mass:=0.1} #Mass of the gas particles (Msun)
clump=${clump:=1} #1: absorbing clump, 0: the same box without it
star_mass=${star_mass:=29.7} #Star mass (Msun)
n_steps=${n_steps:=64} #Index of the last time step
run_analysis=${run_analysis:=1} #Plot the result when the run ends
radiation_table=${radiation_table:=PopII_parsec_spectral_with_pe_photon_energy.hdf5} #Table with MeanPhotonEnergyPE
swift_exe=${swift_exe:="$PWD/../../../../../swift"} #swift binary
run_dir=${run_dir:=run_clump${clump}} #Isolated directory of the run
extra_args=${extra_args:=""} #Further swift arguments, e.g. "-P Snapshots:delta_time:3e-5"

glass_file="glassCube_64.hdf5"
if [ ! -e "$glass_file" ]
then
    echo "Fetching initial glass file (${glass_file})..."
    ./getGlass.sh 64
fi

if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    "$scripts_location"/getGrackleCoolingTable.sh
fi

"$scripts_location"/getRadiationTable.sh "$radiation_table" || exit 1
"$scripts_location"/checkRadiationTable.sh "$radiation_table" --with-isrf || exit 1

# The packet transport needs the mean photon energy of the PE band.
python3 -c "
import sys, h5py
sys.exit(0 if 'MeanPhotonEnergyPE' in h5py.File(sys.argv[1], 'r')['Data/Radiation'] else 1)" "$radiation_table" || {
    echo "$radiation_table has no MeanPhotonEnergyPE dataset. See the ISRF documentation." >&2
    exit 1
}

here=$PWD
if [ -d "$run_dir" ]
then
    echo "$run_dir exists. Its content will be removed."
    rm -r "$run_dir"
fi
mkdir -p "$run_dir/snap"
cd "$run_dir"
ln -s "$here/$glass_file" "$here/CloudyData_UVB=HM2012.h5" "$here/$radiation_table" .
cp "$here/params.yml" .

echo "Generating initial conditions to run the example..."
python3 "$here/makeIC.py" --glass "$glass_file" --rho $gas_density --mass $gas_mass --star_mass $star_mass --table "$radiation_table" --clump $clump -o ICs_isrf_shadow.hdf5

echo "Running simulation..."

# --external-gravity with no Potential: block gives zero gravitational force.
"$swift_exe" --hydro --stars --external-gravity --feedback --cooling \
	     --sync --limiter --verbose=0 --threads=$n_threads \
	     -n $((n_steps - 1)) \
	     -P GEARFeedback:yields_table:$radiation_table \
	     -P GEARFeedback:yields_table_first_stars:$radiation_table \
	     $extra_args params.yml 2>&1 | tee output.log

if [ "$run_analysis" = "1" ]
then
    python3 "$here/plotShadow.py" 
fi
