#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the SmoothedMetallicity example..."
    ./getGlass.sh
fi
if [ ! -e smoothed_metallicity.hdf5 ]
then
    echo "Generating initial conditions for the SmoothedMetallicity example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --steps=1 --hydro --threads=4 smoothed_metallicity.yml 2>&1 | tee output.log

# Plot the solution
python plotSolution.py 1
