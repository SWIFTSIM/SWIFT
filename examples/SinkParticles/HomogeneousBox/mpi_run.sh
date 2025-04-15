#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_ranks=${n_ranks:=2}  #Number of MPI ranks
n_threads=${n_threads:=1}  #Number of threads to use
level=${level:=5}  #Number of particles = 2^(3*level)
jeans_length=${jeans_length:=0.250}  #Jeans wavelenght in unit of the boxsize

#Create the ICs if they do not exist
if [ ! -e ICs_homogeneous_box.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --level $level -o ICs_homogeneous_box.hdf5 --lJ $jeans_length
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

mpirun -n $n_ranks ../../../swift_mpi --pin --hydro --sinks --stars --self-gravity --feedback --cooling --sync --limiter --threads=1 --verbose 0 params.yml 2>&1 | tee output.log

#Do some data analysis to show what's in this box
python3 plot_gas_density.py -i 282 -s 'snap/snapshot'
python3 rhoTPlot.py -i 282 -s 'snap/snapshot'
python3 rhoTPlot.py -i 0 -f 282 -s 'snap/snapshot'
python3 plot_gas_density.py -i 0 -f 282 -s 'snap/snapshot'
