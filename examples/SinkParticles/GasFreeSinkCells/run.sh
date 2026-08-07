#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_ranks=${n_ranks:=0}      # Number of ranks to use
n_threads=${n_threads:=8}  # Number of threads to use
level=${level:=5}  # Nominal gas resolution before carving out the corner: N = 2^(3*level)
jeans_length=${jeans_length:=0.250} # Jeans wavelength in unit of the boxsize
debug=${debug:=0}
run_name=${run_name:=""}
with_star_formation=${with_star_formation=0}

# This example is sinks-only: GEAR star formation with sinks is untested and unsupported
if [[ "$with_star_formation" -ne 0 ]]; then
    echo "Error: This example is sinks-only. GEAR star formation combined with sinks has never been tested." >&2
    exit 1
fi

IC_FILE=ICs_gas_free_sink_cells.hdf5

# Remove the ICs
if [ -e $IC_FILE ]
then
    rm $IC_FILE
fi

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
    ./getGrackleCoolingTable.sh
fi

if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    ./getChemistryTable.sh
fi

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

if [[ $n_ranks -gt 0 ]]; then
  swift="mpirun -n $n_ranks ../../../swift_mpi"
else
  swift="../../../swift"
fi

# --sinks is mandatory: this example's ICs already contain sink particles.
runtime_params="--sinks"

printf "Running simulation..."
$swift --hydro $runtime_params --stars --self-gravity --feedback \
       --cooling --sync --limiter --threads=$n_threads \
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
