#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the 2D vacuum expansion example..."
    ./getGlass.sh
fi
if [ ! -e vacuum.hdf5 ]
then
    echo "Generating initial conditions for the 2D vacuum expansion example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=4 vacuum.yml 2>&1 | tee output.log

# Get the 1D high resolution reference result if not present.
if [ ! -e vacuumSpherical2D_exact.txt ]
then
    echo "Fetching reference solution for the 2D vacuum expansion test..."
    ./getReference.sh
fi

# Plot the result
python plotSolution.py 1
