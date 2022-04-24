#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the Noh problem..."
    ./getGlass.sh
fi
if [ ! -e noh.hdf5 ]
then
    echo "Generating initial conditions for the Noh problem..."
    python3 makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=2 noh.yml 2>&1 | tee output.log

# Plot the solution
python3 plotSolution.py 12
