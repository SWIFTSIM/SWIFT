#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../GEAR_ICs_and_SCRIPTS"

set -eo pipefail

n_threads=${n_threads:=8}  #Number of threads to use
gas_density=${gas_density:=1e3} #Gas density in atom/cm^3
gas_particle_mass=${gas_mass:=0.1} #Mass of the gas particles (Msun)
star_mass=${star_mass:=29.7} #Star mass (Msun)
star_type=${star_type:="single_star"}
level=${level:=6} #Resolution level: N = (2**level)**3 gas particles
time_end=${time_end:=1.9e-4} #TimeIntegration:time_end override (internal units)
dt_max=${dt_max:=1e-5} #TimeIntegration:dt_max override (internal units)
delta_time=${delta_time:=9.5e-6} #Snapshots:delta_time override (internal units)
initial_metallicity=${initial_metallicity:=1} #GEARChemistry:initial_metallicity override (Z/Zsun)
bulk_temperature_K=${bulk_temperature_K:=1000} #Bulk gas InternalEnergy, via makeIC.py
hot_particle_temperature_K=${hot_particle_temperature_K:=3e4} #Pinned-neighbour variant (see README)
alpha_max=${alpha_max:=0.5} #GEARFeedback:LW_FUV_dissipation_alpha_max override
run_name=${run_name:=""}

# Remove the ICs
if [ -e ICs_isrf_dissipation.hdf5 ]
then
    rm ICs_isrf_dissipation.hdf5
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
	--star_type $star_type --bulk-temperature-K $bulk_temperature_K \
	--hot-particle-temperature-K $hot_particle_temperature_K \
    -o ICs_isrf_dissipation.hdf5)
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
# force. SPH:initial_temperature stays 0 (params.yml default) so
# makeIC.py's own per-particle InternalEnergy -- including the heated
# particle -- survives start-up: this example is always the
# pinned-neighbour variant (see README).
../../../../swift --hydro --stars --external-gravity --feedback --cooling \
		   --sync --limiter --verbose=0 --threads=$n_threads \
		   -P TimeIntegration:time_end:$time_end \
		   -P TimeIntegration:dt_max:$dt_max \
		   -P Snapshots:delta_time:$delta_time \
		   -P GEARChemistry:initial_metallicity:$initial_metallicity \
		   -P SPH:initial_temperature:0 \
		   -P GEARFeedback:LW_FUV_dissipation_alpha_max:$alpha_max \
		   params.yml 2>&1 | tee output.log

# Check sign closure and causal-reach transport fidelity (see README).
# hot_particle_id.txt only exists for the pinned-neighbour variant
# (hot_particle_temperature_K != 0).
if [ -e hot_particle_id.txt ]; then
    python3 isrf_dissipation_check.py --hot-particle-id $(cat hot_particle_id.txt)
else
    python3 isrf_dissipation_check.py
fi

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
	[ -e hot_particle_id.txt ] && mv hot_particle_id.txt $run_name
    fi
fi
