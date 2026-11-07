#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location=${scripts_location:=../../../../GEAR_ICs_and_SCRIPTS}

# Make run.sh fail if a subcommand fails; swift's exit code is piped into tee.
set -eo pipefail

n_threads=${n_threads:=8}  #Number of threads to use
leg=${leg:=dust} #redshift: expanding box, no dust; dust: expanding box with dust; static: static box with dust
redshift=${redshift:=0} #Starting redshift of the expanding box
gas_density=${gas_density:=1} #Physical gas density in atom/cm^3
gas_mass=${gas_mass:=5.12e8} #Mass of the gas particles (Msun)
temperature=${temperature:=100} #Gas temperature (K), fixed: the cooling is off
star_mass=${star_mass:=3} #Star mass (Msun)
duration=${duration:=0.184} #Proper time of the run (internal units, 0.184 is 180 Myr)
n_steps=${n_steps:=2000} #Number of time steps
n_snapshots=${n_snapshots:=40} #Number of snapshots after the first
radiation_table=${radiation_table:=PopII_parsec_spectral_with_pe_photon_energy.hdf5} #Table with MeanPhotonEnergyPE
swift_exe=${swift_exe:="$PWD/../../../../../swift"} #swift binary
run_dir=${run_dir:=run_${leg}} #Isolated directory of the run
extra_args=${extra_args:=""} #Further swift arguments, e.g. "-P Snapshots:compression:4"

star_age=1e-6 #Age of the star at the start (internal units)
glass_file="glassCube_32.hdf5"

case "$leg" in
    redshift) cosmology=1; dust_args=(-P GrackleCooling:local_dust_to_gas_ratio:1e-30) ;;
    dust)     cosmology=1; dust_args=() ;;
    static)   cosmology=0; dust_args=() ;;
    *) echo "leg must be redshift, dust or static." >&2; exit 1 ;;
esac

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

# The table is not downloadable: build it with the pychem table generator.
if [ ! -f "$radiation_table" ]
then
    echo "Radiation table '$radiation_table' not found. Set radiation_table=/path/to/table.hdf5." >&2
    echo "It must carry MeanPhotonEnergyPE and Integrated_MeanPhotonEnergyPE in Data/Radiation (see README)." >&2
    exit 1
fi
radiation_path=$(realpath "$radiation_table")
radiation_table=$(basename "$radiation_path")

"$scripts_location"/checkRadiationTable.sh "$radiation_path" --with-isrf || exit 1

# The packet transport needs the mean photon energy of the PE band.
for dataset in MeanPhotonEnergyPE Integrated_MeanPhotonEnergyPE
do
    python3 -c "
import sys, h5py
sys.exit(0 if sys.argv[2] in h5py.File(sys.argv[1], 'r')['Data/Radiation'] else 1)" "$radiation_path" "$dataset" || {
        echo "$radiation_path has no $dataset dataset in Data/Radiation (see README)." >&2
        exit 1
    }
done

here=$PWD
if [ -d "$run_dir" ]
then
    echo "$run_dir exists. Its content will be removed."
    rm -r "$run_dir"
fi
mkdir -p "$run_dir/snap"
cd "$run_dir"
ln -s "$here/$glass_file" "$here/CloudyData_UVB=HM2012.h5" "$radiation_path" .
cp "$here/params.yml" .

echo "Generating initial conditions to run the example..."
python3 "$here/makeIC.py" --glass "$glass_file" --rho $gas_density --mass $gas_mass --temperature $temperature \
	--cosmology $cosmology --redshift $redshift --star_mass $star_mass --star_age $star_age \
	-o ICs_isrf_cosmology.hdf5

# The time line: a_begin, a_end, dt_max and delta_time (expanding box), or time_end, dt_max and delta_time.
eval "$(python3 "$here/cosmologyTimeline.py" --cosmology $cosmology --redshift $redshift \
	--duration $duration --snapshots $n_snapshots --steps $n_steps)"

if [ "$cosmology" = 1 ]
then
    time_args=(--cosmology -P Cosmology:a_begin:$a_begin -P Cosmology:a_end:$a_end
	       -P Snapshots:scale_factor_first:$a_begin -P Statistics:scale_factor_first:$a_begin)
else
    time_end=$(python3 -c "print(repr($star_age + $time_end))")
    time_args=(-P TimeIntegration:time_begin:$star_age -P TimeIntegration:time_end:$time_end
	       -P Snapshots:time_first:$star_age -P Statistics:time_first:$star_age)
fi

echo "Running simulation..."

# --external-gravity with no Potential: block gives zero gravitational force.
"$swift_exe" --hydro --stars --external-gravity --feedback --cooling \
	     --sync --limiter --verbose=0 --threads=$n_threads \
	     "${time_args[@]}" \
	     -P TimeIntegration:dt_max:$dt_max -P Snapshots:delta_time:$delta_time \
	     -P Statistics:delta_time:$delta_time \
	     -P GEARFeedback:yields_table:$radiation_table \
	     -P GEARFeedback:yields_table_first_stars:$radiation_table \
	     "${dust_args[@]}" $extra_args params.yml 2>&1 | tee output.log
