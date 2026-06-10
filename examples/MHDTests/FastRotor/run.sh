#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./glassCube_32.hdf5 ]
then
    echo "Fetching a unit cell, copies of which are to be stacked to generate the ICs ..."
    ./getGlass.sh
fi

if [ ! -f ./FastRotor.hdf5 ]
then
    echo "Generating initial conditions for the Fast Rotor example..."
    python3 makeIC.py -r 4
fi

# Run the example with SWIFT 
../../../swift --hydro --threads=4 FastRotor.yml 2>&1 | tee output.log

# Plot the calculated solution at time t=0.15
python3 plotSolution.py FastRotor_0015.hdf5 FastRotor_0015.png	
