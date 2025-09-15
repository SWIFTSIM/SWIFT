#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=${n_threads:=8}  #Number of threads to use
level=${level:=6}  #Number of particles = 2^(3*level)
gas_density=${gas_density:=1} #Gas density in atom/cm^3
gas_particle_mass=${gas_particle_mass:=10} #Mass of the gas particles
star_mass=${star_mass:=29.7} #Mass of the gas particles
run_name=${run_name:=""}

# Remove the ICs
if [ -e ICs_homogeneous_box.hdf5 ]
then
    rm ICs_homogeneous_box.hdf5
fi

#Create the ICs if they do not exist
if [ ! -e ICs_homogeneous_box.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --level $level --rho $gas_density \
	    --mass $gas_particle_mass --star_mass $star_mass \
	    -o ICs_disc.hdf5
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ./getGrackleCoolingTable.sh
fi

if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    ./getChemistryTable.sh
fi

# Create output directory
DIR=snap #First test of units conversion
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. Its content will be removed."
    rm -r $DIR
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

printf "Running simulation..."
# NOTE: run without cooling
../../../swift --hydro --stars --external-gravity --feedback \
	       --sync --limiter --threads=$n_threads \
	       params.yml 2>&1 | tee output.log

python3 plot_momentum_conservation.py

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
	mv *.png $run_name
	mv *.mp4 $run_name
    fi
fi
