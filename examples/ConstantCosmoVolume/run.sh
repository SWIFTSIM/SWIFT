#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e constantBox.hdf5 ]
then
    echo "Generating initial conditions for the uniform cosmo box example..."
    python makeIC.py
fi

# Run SWIFT
../swift -s -c -G -t 8 constant_volume.yml 2>&1 | tee output.log

# Plot the result
python plotSolution.py $i
