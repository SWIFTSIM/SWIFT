#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

if [ ! -f 'randomized-sine.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../swift \
    --hydro \
    --threads=9 \
    --verbose=0  \
    --radiation \
    --self-gravity \
    --stars \
    --feedback \
    ./randomized-rt.yml 2>&1 | tee output.log

echo "running sanity checks"
python3 ../UniformBox_3D/rt_sanity_checks.py | tee sanity_check.log
