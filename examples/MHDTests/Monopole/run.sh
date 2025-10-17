#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching glass cubes to be used in IC generation for the monopole advection test ..."
    ./getGlass.sh
fi

if [ ! -f ./Monopole.hdf5 ]
then
    echo "Generating initial conditions for the monopole advection test ..."
    python3 makeIC.py
fi

# Run the example with SWIFT 
../../../swift --hydro --threads=4 Monopole.yml 2>&1 | tee output.log

# Plot the time evolution of the average divergence error
python3 divergence_evolution.py statistics.txt 
