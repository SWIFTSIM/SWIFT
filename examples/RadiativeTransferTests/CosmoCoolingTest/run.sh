#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f ./cooling_test.hdf5 ]; then
    echo "creating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swift \
    --hydro \
    --cosmology \
    --threads=4 \
    --verbose=0  \
    --radiation \
    --external-gravity \
    --stars \
    --feedback \
./rt_cooling_test.yml 2>&1 | tee output.log

# Wanna run with cooling, but no RT? This should do the trick
# ../../../swift \
#     --hydro \
#     --threads=4 \
#     --verbose=0  \
#     --cooling \
#     ./rt_cooling_test.yml 2>&1 | tee output.log

if [ ! -f "CosmoRTCoolingTestReference.txt" ]; then
    ./getReference.sh
fi
python3 plotSolution.py
