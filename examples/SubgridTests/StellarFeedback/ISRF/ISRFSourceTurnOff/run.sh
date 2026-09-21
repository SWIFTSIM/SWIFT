#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

# make run.sh fail if a subcommand fails. pipefail matters here specifically:
# swift's own exit code is piped into `tee output.log`, and without it a
# crashed/errored swift run still lets the pipeline "succeed" (tee's own
# exit code), silently producing a run_name output directory that looks
# complete but only contains a startup-error log.
set -eo pipefail

n_threads=${n_threads:=8}  #Number of threads to use
gas_density=${gas_density:=20} #Gas density in atom/cm^3
gas_particle_mass=${gas_mass:=0.1} #Mass of the gas particles (Msun)
star_mass=${star_mass:=70} #Star mass (Msun)
star_type=${star_type:="single_star"}
level=${level:=5} #Resolution level: N = (2**level)**3 gas particles
time_end=${time_end:=1.239684e-2} #TimeIntegration:time_end override (internal units)
dt_max=${dt_max:=1e-5} #TimeIntegration:dt_max override (internal units)
initial_metallicity=${initial_metallicity:=1} #GEARChemistry:initial_metallicity override (Z/Zsun)
run_name=${run_name:=""}

# The shipped output_list_isrf_source_turn_off.txt is tied to star_mass and
# initial_metallicity (they set the death time) and to gas_density, gas_mass
# and dt_max (they set the post-death decay timescale, through c_hyp =
# C_hyp*h/dt: see README). level only sets the box size. Overriding any
# of the former without regenerating that file leaves the snapshot cadence
# mismatched to the new death/decay times; time_end must not be shortened
# below the list's last entry.

# Remove the ICs
if [ -e ICs_isrf_source_turn_off.hdf5 ]
then
    rm ICs_isrf_source_turn_off.hdf5
fi

glass_n=$((2**level))
glass_file="glassCube_${glass_n}.hdf5"
if [ ! -e "$glass_file" ]
then
    echo "Fetching initial glass file (${glass_file})..."
    ./getGlass.sh $glass_n
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    "$scripts_location"/getGrackleCoolingTable.sh
fi

if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    "$scripts_location"/getChemistryTable.sh
fi

# Stop here on a table the radiation reader cannot use, rather than
# aborting at start-up once the initial conditions are built.
"$scripts_location"/checkRadiationTable.sh POPIIsw.h5 --with-isrf || exit 1

echo "Generating initial conditions to run the example..."
ic_output=$(python3 makeIC.py --level $level --rho $gas_density \
	--mass $gas_particle_mass --star_mass $star_mass \
	--star_type $star_type \
    -o ICs_isrf_source_turn_off.hdf5)
echo "$ic_output"

# Create output directory
DIR=snap
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. Its content will be removed."
    rm -r $DIR
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

printf "Running simulation..."

# --external-gravity with no Potential: block gives zero gravitational
# force. --cooling is on so cooling_init() actually resolves chemistry_data
# (including local_dust_to_gas_ratio); the per-particle energy update
# itself is skipped via GrackleCooling:disable_cooling_for_debugging in
# params.yml, which isolates the ISRF field from gas thermal/dynamical
# response without leaving chemistry_data unresolved (see README).
../../../../../swift --hydro --stars --external-gravity --feedback --cooling \
		   --sync --limiter --verbose=0 --threads=$n_threads \
		   -P TimeIntegration:time_end:$time_end \
		   -P TimeIntegration:dt_max:$dt_max \
		   -P GEARChemistry:initial_metallicity:$initial_metallicity \
		   params.yml 2>&1 | tee output.log

moved=0
if [ -z "$run_name" ]; then
    echo "run_name is empty."
else
    if [ -d "$run_name" ]; then
	echo "$run_name directory exists. Nothing will be moved."
    else
	echo "$run_name directory does not exists. It will be created."
	mkdir -p $run_name
	mv snap $run_name
	mv output.log $run_name
	mv timesteps.txt $run_name
	mv statistics.txt $run_name
	mv unused_parameters.yml $run_name
	mv used_parameters.yml $run_name
	moved=1
    fi
fi

# Check that the ISRF field decays away after the star's death instead of
# persisting (see README). Run last, after any run_name move above, and
# point --snapshot/--timesteps-log/--used-parameters at the moved paths;
# running this check manually after a move needs the same repointing.
if [ "$moved" -eq 1 ]; then
    check_dir="$run_name/"
else
    check_dir=""
fi
rc=0
python3 isrf_source_turn_off_check.py \
    --snapshot "${check_dir}snap/snapshot_*.hdf5" \
    --timesteps-log "${check_dir}timesteps.txt" \
    --used-parameters "${check_dir}used_parameters.yml" \
    --output "${check_dir}isrf_source_turn_off_check.png" || rc=$?
echo "isrf_source_turn_off_check.py exit code: $rc"
exit $rc
