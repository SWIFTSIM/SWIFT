#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_L0050N0752_ICs.hdf5 ]
then
    echo "Fetching initial conditions for the quick Lyman-alpha 50Mpc example..."
    ./getIC.sh
fi

if [ ! -e coolingtables ]
then
    echo "Fetching quick Lyman-alpha cooling tables..."
    ../getQLACoolingTable.sh
fi

# The following run-time options are broken down by line as:
# Basic run-time options
# Threading options - run with threads and pinning (latter not required but improves performance)
# The corresponding parameter file for this run

../../swift \
    --cosmology --quick-lyman-alpha \
    --threads=16 --pin \
    qla_50.yml
