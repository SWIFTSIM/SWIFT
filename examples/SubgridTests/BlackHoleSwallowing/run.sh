#!/bin/bash

# This example requires SWIFT to be configured with:
#   --with-subgrid=EAGLE

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the black hole swallowing example..."
    ./getGlass.sh
fi
if [ ! -e bh_swallowing.hdf5 ]
then
    echo "Generating initial conditions for the black hole swallowing example..."
    python3 makeIC.py
fi

# Get the cooling tables (the cooling model is initialised by --temperature)
if [ ! -e coolingtables ]
then
    echo "Fetching EAGLE cooling tables..."
    ./getCoolingTable.sh
fi

# Run SWIFT
../../../swift --hydro --black-holes --self-gravity --temperature --limiter --sync --threads=4 swallowing.yml 2>&1 | tee output.log

# The same run over MPI, to compare with the serial one:
# mpirun -n 4 ../../../swift_mpi --hydro --black-holes --self-gravity --temperature --limiter --sync --threads=4 swallowing_mpi.yml 2>&1 | tee output_mpi.log
