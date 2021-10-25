#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the Sod shock example..."
    ../../HydroTests/SodShock_3D/getGlass.sh
fi
if [ ! -e sodShock.hdf5 ]
then
    echo "Generating initial conditions for the Sod shock example..."
    python3 makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=4 sodShock.yml 2>&1 | tee output.log

python3 ../../HydroTests/SodShock_3D/plotSolution.py 1
