#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

if [ ! -f 'uniformBox-rt.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    uniform_rt_timestep_output_sync.yml 2>&1 | tee output.log

echo "running sanity checks"
python3 ./rt_sanity_checks.py
echo "running checks for uniform box test"
python3 ./rt_uniform_box_checks.py
