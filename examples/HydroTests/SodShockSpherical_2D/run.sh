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

# Get the high resolution 1D reference solution if not present.
if [ ! -e sodShockSpherical2D_exact.txt ]
then
    echo "Fetching reference solution for 2D spherical Sod shock example..."
    ./getReference.sh
fi
python plotSolution.py 1
