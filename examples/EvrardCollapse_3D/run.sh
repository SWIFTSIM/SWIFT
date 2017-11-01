#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e evrard.hdf5 ]
then
    echo "Generating initial conditions for the Evrard collapse example..."
    python makeIC.py
fi

# Run SWIFT
../swift -s -G -t 4 evrard.yml 2>&1 | tee output.log

# Plot the solution
python plot_density_profile.py 4
