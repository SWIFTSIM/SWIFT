#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

if [ ! -e dSph_cusp_SIDM.hdf5 ]; then
    echo "Fetching initial conditions to run the example..."
    wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/IsolateddSph/dSph_cusp.hdf5
    python3 redo_ICs.py dSph_cusp.hdf5 dSph_cusp_SIDM.hdf5
fi

printf "Running simulation..."

../../../swift --sidm --self-gravity --threads=14 params.yml 2>&1 | tee output.log