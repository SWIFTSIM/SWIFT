#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e evrard_disc.hdf5 ]
then
    echo "Generating initial conditions for the Evrard disc example..."
    python makeIC.py
fi

# Run SWIFT
../swift -s -G -t 4 evrard.yml 2>&1 | tee output.log
