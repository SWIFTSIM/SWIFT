#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e OrszagTangVortex.hdf5 ]
then
    echo "Generating initial conditions for the OrszagTang " \
         "example..."
    python makeIC.py
fi

# Run SWIFT
rm -I OrszagTangVortex_0???.hdf5
../../swift --hydro --threads=8 OrszagTangVortex.yml 2>&1 | tee output.log

julia plotSolution.jl OrszagTangVortex_0005.hdf5
