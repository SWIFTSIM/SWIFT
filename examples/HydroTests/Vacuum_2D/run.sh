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
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=4 vacuum.yml 2>&1 | tee output.log

# Plot the result
python3 plotSolution.py 1
