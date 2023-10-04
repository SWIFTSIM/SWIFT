#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e noh.hdf5 ]
then
    echo "Generating initial conditions for the Noh problem..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=1 noh.yml 2>&1 | tee output.log

# Plot the solution
python3 plotSolution.py 12
