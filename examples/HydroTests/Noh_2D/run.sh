#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the Noh problem..."
    ./getGlass.sh
fi
if [ ! -e noh.hdf5 ]
then
    echo "Generating initial conditions for the Noh problem..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=2 noh.yml 2>&1 | tee output.log

# Plot the solution
python plotSolution.py 12
