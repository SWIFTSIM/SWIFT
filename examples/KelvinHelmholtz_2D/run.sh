#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e kelvinHelmholtz.hdf5 ]
then
    echo "Generating initial conditions for the Kelvin-Helmholtz example..."
    python makeIC.py
fi

# Run SWIFT
../swift -s -t 1 kelvinHelmholtz.yml

# Plot the solution
python plotSolution.py 6
