#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location=${scripts_location:=../../../../GEAR_ICs_and_SCRIPTS}

# Make run.sh fail if a subcommand fails; swift's exit code is piped into tee.
set -eo pipefail

n_threads=${n_threads:=4}  #Number of threads to use
arm=${arm:=one} #one: one time bin, capped light speed; two: two time bins, fixed light speed
gas_density=${gas_density:=1e3} #Gas density in atom/cm^3
gas_mass=${gas_mass:=0.1} #Mass of the gas particles (Msun)
speed_fraction=${speed_fraction:=3e-5} #Light speed as a fraction of c
energy_scale=${energy_scale:=1} #Factor on the energy of the pulse
snapshot_every=${snapshot_every:=8} #Snapshot interval in base time steps (9.375e-6 internal units)
radiation_table=${radiation_table:=PopII_parsec_spectral_with_pe_photon_energy.hdf5} #Table with MeanPhotonEnergyPE
swift_exe=${swift_exe:="$PWD/../../../../../swift"} #swift binary
run_dir=${run_dir:=run_${arm}} #Isolated directory of the run
extra_args=${extra_args:=""} #Further swift arguments, e.g. "-P Snapshots:compression:4"

# The two arms: the light speed mode, the margin of the time step bound and the number of
# swift steps. With the margin 0.4305 about half of the gas takes half of the step.
case "$arm" in
    one) light_speed=${light_speed:=cap};   margin=${margin:=0.5};    n_steps=${n_steps:=80} ;;
    two) light_speed=${light_speed:=fixed}; margin=${margin:=0.4305}; n_steps=${n_steps:=96} ;;
    *) echo "arm must be one or two." >&2; exit 1 ;;
esac

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
python3 "$here/makeIC.py" --glass "$glass_file" --rho $gas_density --mass $gas_mass --energy_scale $energy_scale -o ICs_isrf_propagation_speed.hdf5
mapfile -t beam_args < beam_args.txt
snapshot_dt=$(python3 -c "print('%.9e' % ($snapshot_every * 9.375e-6))")

echo "Running simulation..."

# --external-gravity with no Potential: block gives zero gravitational force.
"$swift_exe" --hydro --stars --external-gravity --feedback --cooling \
	     --sync --limiter --verbose=0 --threads=$n_threads \
	     -n $n_steps \
	     "${beam_args[@]}" \
	     -P GEARFeedback:ISRF_packet_light_speed:$light_speed \
	     -P GEARFeedback:ISRF_packet_light_speed_fraction_of_c:$speed_fraction \
	     -P GEARFeedback:ISRF_packet_light_speed_time_step_margin:$margin \
	     -P Snapshots:delta_time:$snapshot_dt \
	     -P GEARFeedback:yields_table:$radiation_table \
	     -P GEARFeedback:yields_table_first_stars:$radiation_table \
	     $extra_args params.yml 2>&1 | tee output.log
