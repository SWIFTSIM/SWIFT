#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e stellar_evolution.hdf5 ]
then
    echo "Generating initial conditions for the 3D stellar evolution example..."
    python makeIC.py
fi

# Run SWIFT
../swift --limiter --feedback --stars --hydro --external-gravity --threads=8 stellar_evolution.yml 2>&1 | tee output.log

