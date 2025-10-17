#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the shocktube example..."
    ./getGlass.sh
fi
if [ ! -e RyuJones1A.hdf5 ]
then
    echo "Generating initial conditions for the shock example..."
    python makeIC.py
fi

# Run SWIFT
rm -I RJ1A_????.hdf5 
../../../swift --hydro --threads=8 rj1A.yml 2>&1 | tee output.log

#python plotSolution.py 1
julia plotSolution.jl RJ1A_0010.hdf5
