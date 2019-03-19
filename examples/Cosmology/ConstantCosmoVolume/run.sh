#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e gravity_glassCube_32.hdf5 ]
then
    echo "Fetching initial grvity glass file for the constant cosmological box example..."
    ./getGlass.sh
fi
if [ ! -e constantBox.hdf5 ]
then
    echo "Generating initial conditions for the uniform cosmo box example..."
    python3 makeIC.py
fi

# Run SWIFT
../../swift --hydro --cosmology --self-gravity --threads=8 constant_volume.yml 2>&1 | tee output.log

# Plot the result
python3 plotSolution.py $i
