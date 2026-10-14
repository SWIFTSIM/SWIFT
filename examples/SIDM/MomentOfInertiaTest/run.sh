#!/bin/bash

# Runs the EAGLE 6Mpc ICs once with the DM particles converted to SIDM particles
# and once with them converted to gas particles, and compares the
# InertiaTensors written in the snapshots.

IC_DIR=../../EAGLE_low_z/EAGLE_6

if [ ! -e ${IC_DIR}/EAGLE_ICs_6.hdf5 ]; then
    echo "Fetching initial conditions for the EAGLE 6Mpc example..."
    (cd ${IC_DIR} && ./getIC.sh)
fi

if [ ! -e EAGLE_ICs_6_gas.hdf5 ]; then
    echo "Converting DM particles to gas particles..."
    python3 make_ICs.py ${IC_DIR}/EAGLE_ICs_6.hdf5 EAGLE_ICs_6_gas.hdf5 gas
fi

if [ ! -e EAGLE_ICs_6_SIDM.hdf5 ]; then
    echo "Converting DM particles to SIDM particles..."
    python3 make_ICs.py ${IC_DIR}/EAGLE_ICs_6.hdf5 EAGLE_ICs_6_SIDM.hdf5 sidm
fi

echo "Running with SIDM..."
../../../swift --sidm --self-gravity --threads=4 -n 1 params_sidm.yml 2>&1 | tee output_sidm.log

echo "Running with gas..."
../../../swift --hydro --self-gravity --threads=4 -n 1 params_gas.yml 2>&1 | tee output_gas.log

echo "Comparing inertia tensors..."
python3 compare_inertia.py --sidm snap_sidm/snapshot_0000.hdf5 --gas snap_gas/snapshot_0000.hdf5
