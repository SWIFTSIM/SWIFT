#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the BrioWu example..."
    ./getGlass.sh
fi
if [ ! -e BrioWu.hdf5 ]
then
    echo "Generating initial conditions for the Brio Wu shock example..."
    python3 makeIC.py
fi

# Run SWIFT
rm -I BrioWu_0???.hdf5 
../../../swift --hydro --threads=8 BrioWu.yml 2>&1 | tee output.log

#python plotSolution.py 1
#julia plotSolution.jl BrioWu_0008.hdf5
