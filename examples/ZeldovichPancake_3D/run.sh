#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e zeldovichPancake.hdf5 ]
then
    echo "Generating initial conditions for the 3D Zeldovich pancake example..."
    python makeIC.py
fi

# Run SWIFT
../swift -s -c -G -v 1 -t 1 zeldovichPancake.yml 2>&1 | tee output.log

# Plot the result
#python plotSolution.py 1
