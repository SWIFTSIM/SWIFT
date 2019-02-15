#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the Sedov blast example..."
    ./getGlass.sh
fi
if [ ! -e sedov.hdf5 ]
then
    echo "Generating initial conditions for the Sedov blast example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --limiter --threads=4 sedov.yml 2>&1 | tee output.log

# Plot the solution
python plotSolution.py 5
