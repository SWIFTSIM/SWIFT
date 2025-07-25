#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./FCC_glassCube_32.hdf5 ]
then
    echo "Generating a BCC unit cell, copies of which are to be stacked to generate the Ryu & Jones 1A shock tube ..."
    python3 makeBCC.py -n 32
fi
if [ ! -f ./RyuJones_1A.hdf5 ]
then
    echo "Generating initial conditions for the Ryu & Jones 1A shock tube example..."
    python3 makeIC.py
fi

# Run the example with SWIFT 
../../../swift --hydro --threads=4 RyuJones_1A.yml 2>&1 | tee output.log

# Plot the calculated solution at time t=0.08
python3 plotSolution.py RyuJones_0008.hdf5 RyuJones_0008.png	
