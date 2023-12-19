#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f advection_1D.hdf5 ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swift \
    --hydro \
    --cosmology \
    --threads=4 \
    --verbose=0  \
    --radiation \
    --stars \
    --feedback \
    --external-gravity \
    --fpe \
    $1 2>&1 | tee output.log
