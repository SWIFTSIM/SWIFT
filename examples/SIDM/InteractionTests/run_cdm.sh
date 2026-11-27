#!/bin/bash

# make run_cdm.sh fail if a subcommand fails
set -e

if [ ! -e dSph_cusp_fixed.hdf5 ]; then
    if [ ! -e dSph_cusp.hdf5 ]; then
        echo "Fetching initial conditions to run the example..."
        wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/IsolateddSph/dSph_cusp.hdf5
    fi
    python3 redo_ICs_CDM_PID.py # fix PID=0 issue
fi

printf "Running simulation...\n"
../../../swift --self-gravity --threads=14 params_cdm.yml 2>&1 | tee output_cdm.log
