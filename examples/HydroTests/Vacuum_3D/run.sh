#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the 3D vacuum expansion example..."
    ./getGlass.sh
fi
if [ ! -e vacuum.hdf5 ]
then
    echo "Generating initial conditions for the 3D vacuum expansion example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=16 vacuum.yml 2>&1 | tee output.log

# Plot the result
python plotSolution.py 1
