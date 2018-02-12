#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e interactingBlastWaves.hdf5 ]
then
    echo "Generating initial conditions for the Sedov blast example..."
    python makeIC.py
fi

# Run SWIFT
../swift -s -t 1 interactingBlastWaves.yml 2>&1 | tee output.log

# Plot the solution
python plotSolution.py 4
