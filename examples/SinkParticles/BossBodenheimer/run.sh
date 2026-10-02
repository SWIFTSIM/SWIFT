#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

scripts_location="../../GEAR_ICs_and_SCRIPTS"

n_ranks=${n_ranks:=0}     # Number of ranks to use
n_threads=${n_threads:=8} # Number of threads to use
N=${N:=500000}            # Number of gas particles

# Create the ICs if they do not exist
if [ ! -e boss_bodenheimer.hdf5 ]; then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py -o boss_bodenheimer.hdf5 -N $N
fi

# Get the Grackle cooling table (high-density variant: gas here reaches
# densities well above the standard table's range)
if [ ! -e CloudyData_UVB=HM2012_high_density.h5 ]; then
    echo "Fetching the Cloudy tables required by Grackle..."
    $scripts_location/getGrackleCoolingTable.sh --high-density
fi

if [ ! -e POPIIsw.h5 ]; then
    echo "Fetching the chemistry tables..."
    $scripts_location/getChemistryTable.sh
fi

# Create output directory
DIR=snap
if [ -d "$DIR" ]; then
    echo "$DIR directory exists. Its content will be removed."
    rm -r $DIR
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

if [[ $n_ranks -gt 0 ]]; then
  swift="mpirun -n $n_ranks ../../../swift_mpi"
else
  swift="../../../swift"
fi

printf "Running simulation..."

$swift --hydro --sinks --stars --self-gravity --feedback --cooling --threads=$n_threads params.yml 2>&1 | tee output.log
