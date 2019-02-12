#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the disc patch example..."
    ./getGlass.sh
fi
if [ ! -e Disc-Patch.hdf5 ]
then
    echo "Generating initial conditions for the disc patch example..."
    python makeIC.py
fi

# Run SWIFT
../../../swift --external-gravity --hydro --threads=4 disc-patch-icc.yml 2>&1 | tee output.log

python plotSolution.py
