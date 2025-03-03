#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=${n_threads:=8}  #Number of threads to use
level=${level:=5}  #Number of particles = 2^(3*level)
jeans_length=${jeans_length:=0.250}  #Jeans wavelenght in unit of the boxsize
gas_density=${gas_density:=0.1} #Gas density in atom/cm^3
gas_particle_mass=${gas_particle_mass:=50} #Mass of the gas particles
n_sinks=${n_sinks:=10}  #Number of sinks

# Remove the ICs
if [ -e ICs_homogeneous_box.hdf5 ]
then
    rm ICs_homogeneous_box.hdf5
fi

#Create the ICs if they do not exist
if [ ! -e ICs_homogeneous_box.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --level $level -o ICs_homogeneous_box.hdf5 --lJ $jeans_length --n_sink $n_sinks --sink_pos 0 0 0 --sinks_vel 10 10 0
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

../../../swift --hydro --sinks --stars --external-gravity --feedback --threads=$n_threads params.yml 2>&1 | tee output.log
