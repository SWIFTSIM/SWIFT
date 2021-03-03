#!/bin/bash

# Fetch the initial conditions if they are not present.
if [ ! -e particles.hdf5 ]
then
    echo "Fetching initial conditions for the neutrino example..."
    ./getIC.sh
fi

# Run SWIFT
../../swift -c -G --threads=8 neutrino_cosmo.yml 2>&1 | tee output.log
