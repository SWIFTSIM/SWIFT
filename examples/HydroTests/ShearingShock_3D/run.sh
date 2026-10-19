#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the shearing shock example..."
    ./getGlass.sh
fi
if [ ! -e shearingShock.hdf5 ]
then
    echo "Generating initial conditions for the shearing shock example..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=4 shearingShock.yml 2>&1 | tee output.log

# Plot the solution
python3 plotSolution.py 2
