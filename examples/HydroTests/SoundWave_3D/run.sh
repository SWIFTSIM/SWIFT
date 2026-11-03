#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the sound wave example..."
    ./getGlass.sh
fi
if [ ! -e soundWave.hdf5 ]
then
    echo "Generating initial conditions for the sound wave example..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=4 soundWave.yml 2>&1 | tee output.log

# Phase speed and damping of the mode
python3 plotSolution.py
