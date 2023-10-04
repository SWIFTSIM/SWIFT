#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the Gresho-Chan vortex example..."
    ./getGlass.sh
fi
if [ ! -e greshoVortex.hdf5 ]
then
    echo "Generating initial conditions for the Gresho-Chan vortex example..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=1 gresho.yml 2>&1 | tee output.log

# Plot the solution
python3 plotSolution.py 11
