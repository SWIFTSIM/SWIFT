#!/bin/bash

set -e

if [ ! -f glassCube_32.hdf5 ]; then
    ./getGlass.sh
fi

if [ ! -f feedback.hdf5 ]; then
    echo "Generating ICs"
    python3 ./makeIC.py
fi

# Run SWIFT
../../../swift --hydro --cooling --limiter --threads=4 feedback.yml 2>&1 | tee output.log

# Plot the solution
python3 plotSolution.py 5
python3 plotEnergy.py
