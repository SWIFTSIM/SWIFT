#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the Toro (2009) test 2 example..."
    ./getGlass.sh
fi
if [ ! -e toroTest2.hdf5 ]
then
    echo "Generating initial conditions for the Toro (2009) test 2 example..."
    python3 makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=1 toroTest2.yml 2>&1 | tee output.log

python3 plotSolution.py 1
