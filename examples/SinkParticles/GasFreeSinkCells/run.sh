#!/bin/bash

# make run.sh fail if a subcommand fails
set -eo pipefail

GEAR_SCRIPTS=../../GEAR_ICs_and_SCRIPTS # Shared GEAR table download scripts

# Force n_rank = 0 until the sink MPI branch is merged
n_ranks=0
# n_ranks=${n_ranks:=0}    # Number of ranks to use
n_threads=${n_threads:=8}  # Number of threads to use
level=${level:=5}  # Nominal gas resolution before carving out the corner: N = 2^(3*level)
jeans_length=${jeans_length:=0.250} # Jeans wavelength in unit of the boxsize
debug=${debug:=0}
run_name=${run_name:=""}

# The ICs depend on the level and Jeans length, so name them after both
IC_FILE=ICs_gas_free_sink_cells_level${level}_lJ${jeans_length}.hdf5

# Create the ICs if they do not exist
if [ ! -e $IC_FILE ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --level $level -o $IC_FILE --lJ $jeans_length
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    $GEAR_SCRIPTS/getGrackleCoolingTable.sh
fi

if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    $GEAR_SCRIPTS/getChemistryTable.sh
fi

# Create output directory
DIR=snap
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. Its content will be removed."
    rm -rf $DIR
else
    echo "$DIR directory does not exists. It will be created."
    mkdir -p $DIR
fi

if [[ $n_ranks -gt 0 ]]; then
  swift="mpirun -n $n_ranks ../../../swift_mpi"
else
  swift="../../../swift"
fi

printf "Running simulation..."
$swift --hydro --sinks --stars --self-gravity --feedback \
       --cooling --sync --limiter --threads=$n_threads \
       --param=InitialConditions:file_name:$IC_FILE \
       params.yml 2>&1 | tee output.log

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
