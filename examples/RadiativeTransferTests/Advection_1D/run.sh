#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f advection_1D.hdf5 ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
# mpirun -n 4 ../../../swift_mpi \
../../../swift \
    --hydro \
    --threads=4 \
    --verbose=0  \
    --radiation \
    --stars \
    --feedback \
    --external-gravity \
    ./rt_advection1D.yml 2>&1 | tee output.log

python3 ./plotSolution.py
