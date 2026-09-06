#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../GEAR_ICs_and_SCRIPTS"

# make run.sh fail if a subcommand fails. pipefail matters here specifically:
# swift's own exit code is piped into `tee output.log`, and without it a
# crashed/errored swift run still lets the pipeline "succeed" (tee's own
# exit code), silently producing a run_name output directory that looks
# complete but only contains a startup-error log.
set -eo pipefail

n_threads=${n_threads:=8}  #Number of threads to use
gas_density=${gas_density:=1e3} #Gas density in atom/cm^3
gas_particle_mass=${gas_mass:=0.1} #Mass of the gas particles (Msun)
star_mass=${star_mass:=29.7} #Star mass (Msun)
star_type=${star_type:="single_star"}
level=${level:=5} #Resolution level: N = (2**level)**3 gas particles
time_end=${time_end:=3e-5} #TimeIntegration:time_end override (internal units)
dt_max=${dt_max:=1e-5} #TimeIntegration:dt_max override (internal units)
delta_time=${delta_time:=1e-5} #Snapshots:delta_time override (internal units)
run_name=${run_name:=""}

# Remove the ICs
if [ -e ICs_isrf_injection_conservation.hdf5 ]
then
    rm ICs_isrf_injection_conservation.hdf5
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

echo "Generating initial conditions to run the example..."
ic_output=$(python3 makeIC.py --level $level --rho $gas_density \
	--mass $gas_particle_mass --star_mass $star_mass \
	--star_type $star_type \
    -o ICs_isrf_injection_conservation.hdf5)
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
# force. --cooling is on so cooling_init() resolves chemistry_data (needed
# for GEARChemistry:initial_metallicity:0's extinction=1 case to actually
# be exercised through the real code path, not skipped as unresolved); the
# per-particle energy update itself is skipped via
# GrackleCooling:disable_cooling_for_debugging in params.yml.
../../../../swift --hydro --stars --external-gravity --feedback --cooling \
		   --sync --limiter --verbose=0 --threads=$n_threads \
		   -P TimeIntegration:time_end:$time_end \
		   -P TimeIntegration:dt_max:$dt_max \
		   -P Snapshots:delta_time:$delta_time \
		   params.yml 2>&1 | tee output.log

# Check that the injected FUV/LW energy sums to Delta_t * L_band per star
# feedback pass (see README).
python3 isrf_injection_conservation_check.py

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
    fi
fi
