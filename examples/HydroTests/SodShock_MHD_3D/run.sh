#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e sodShock.hdf5 ]
then
    echo "Generating initial conditions for the Sod shock example..."
    python makeGlass.py -n 64 -o glassCube_64.hdf5
    python makeGlass.py -n 32 -o glassCube_32.hdf5
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=4 sodShock.yml 2>&1 | tee output.log

python plotSolution.py 1
