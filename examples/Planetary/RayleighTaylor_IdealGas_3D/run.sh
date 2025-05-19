#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e rayleigh_taylor.hdf5 ]
then
    echo "Generating initial conditions for the 3D Rayleigh--Taylor test..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --external-gravity --threads=4 rayleigh_taylor.yml 2>&1 | tee output_rayleigh_taylor.log

# Plot the solutions
python3 ./plotSnapshots.py
