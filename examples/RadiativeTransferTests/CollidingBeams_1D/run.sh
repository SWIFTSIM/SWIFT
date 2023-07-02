#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f collision_1D.hdf5 ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swift \
    --hydro \
    --threads=1 \
    --verbose=0  \
    --radiation \
    --stars \
    --feedback \
    --external-gravity \
    ./rt_collision1D.yml 2>&1 | tee output.log

python3 ./plotEnergies.py
python3 ./plotSolution.py
