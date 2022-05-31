#!/bin/bash

set -e

# Generate the initial conditions if they are not present.
if [ ! -e rayleigh_taylor.hdf5 ]
then
    echo "Generating initial conditions for the Rayleigh Taylor example..."
    python3 makeIC.py
fi

# Run SWIFT
../../swift --hydro --external-gravity  --threads=8 rayleigh_taylor.yml 2>&1 | tee output.log

python3 makeMovie.py -i 0 -f 1001
