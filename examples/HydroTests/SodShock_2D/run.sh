#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the Sod shock example..."
    ./getGlass.sh
fi
if [ ! -e sodShock.hdf5 ]
then
    echo "Generating initial conditions for the Sod shock example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=1 sodShock.yml 2>&1 | tee output.log

python plotSolution.py 1
