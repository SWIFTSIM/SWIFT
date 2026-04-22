#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

if [ ! -e dSph_cusp.hdf5 ]; then
    echo "Fetching initial conditions to run the example..."
    wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/IsolateddSph/dSph_cusp.hdf5
fi

if [ ! -e dSph_cusp_gas.hdf5 ]; then
        python3 redo_ICs_gas.py dSph_cusp.hdf5 dSph_cusp_gas.hdf5
    
fi

printf "Running simulation..."

../../../../swift --hydro --external-gravity --self-gravity -n 1 --threads=14 params_gas.yml 2>&1 | tee output.log

printf "Comparing with SIDM densities..."
python3 compare_gas_SIDM_densities.py --sidm snap/snapshot_0000.hdf5 --gas snap_gas/snapshot_0000.hdf5
