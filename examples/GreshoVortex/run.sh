#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e glass_128.hdf5 ]
then
    echo "Fetching initial glass file for the Gresho-Chan vortex example..."
    ./getGlass.sh
fi
if [ ! -e greshoVortex.hdf5 ]
then
    echo "Generating initial conditions for the Gresho-Chan vortex example..."
    python makeIC.py
fi

# Run SWIFT
../swift -s -t 1 gresho.yml

# Plot the solution
python plotSolution.py 11
