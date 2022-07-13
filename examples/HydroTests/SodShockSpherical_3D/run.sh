#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the Sod shock example..."
    ./getGlass.sh
fi
if [ ! -e sodShock.hdf5 ]
then
    echo "Generating initial conditions for the Sod shock example..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=4 sodShock.yml 2>&1 | tee output.log

# Get the high resolution 1D reference solution if not present.
if [ ! -e sodShockSpherical3D_exact.txt ]
then
    echo "Fetching reference solution for 3D spherical Sod shock example..."
    ./getReference.sh
fi
python3 plotSolution.py 1
