#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e square.hdf5 ]
then
    echo "Generating initial conditions for the square test ..."
    python3 makeICDifferentMasses.py
fi

# Run SWIFT
../../swift --hydro --threads=4 square.yml 2>&1 | tee output.log

# Plot the solution
python3 plotSolution.py 40
python3 makeMovie.py
