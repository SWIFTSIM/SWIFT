#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=${n_threads:=8}  #Number of threads to use
level=${level:=5}  #Number of particles = 2^(3*level)
gas_density=${gas_density:=0.1} #Gas density in atom/cm^3
gas_particle_mass=${gas_particle_mass:=10} #Mass of the gas particles
scale_height=${z_scale:=0.1} #Mass of the gas particles
star_mass=${star_mass:=29.7} #Mass of the gas particles
seed=${seed:=1} # Random seed for random numbers

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
	    --z_scale $scale_height --seed $seed \
	    -o ICs_homogeneous_box.hdf5
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
DIR="snap"
if [ -d "$DIR" ]; then
    echo "$DIR directory exists. Its content will be removed."
    rm -r "$DIR"
fi
mkdir "$DIR"
    
printf "Running simulation..."
../../../swift --hydro --stars --external-gravity --feedback \
	       --sync --limiter --threads=$n_threads \
	       params.yml 2>&1 | tee output.log
