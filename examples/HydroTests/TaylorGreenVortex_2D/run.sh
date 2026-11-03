#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the Taylor-Green vortex example..."
    ./getGlass.sh
fi
if [ ! -e taylorGreen.hdf5 ]
then
    echo "Generating initial conditions for the Taylor-Green vortex example..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=4 taylorGreen.yml 2>&1 | tee output.log

# Plot the decay of the kinetic energy
python3 plotSolution.py
