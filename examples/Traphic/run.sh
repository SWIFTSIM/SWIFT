#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the Traphic example..."
    ./getGlass.sh
fi
if [ ! -e sodShock.hdf5 ]
then
    echo "Generating initial conditions for the Traphic example..."
    python makeIC.py
fi

# Run SWIFT
../swift -R -t 4 traphic.yml 2>&1 | tee output.log
