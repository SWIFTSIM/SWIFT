#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./BCCglassCube_48.hdf5 ]
then
    echo "Generating a BCC unit cell, copies of which are to be stacked to generate the left state of the Brio & Wu shock tube ..."
    python3 makeBCC.py -n 48
fi
if [ ! -f ./BCCglassCube_24.hdf5 ]
then
    echo "Generating a BCC unit cell, copies of which are to be stacked to generate the right state of the Brio & Wu shock tube ..."
    python3 makeBCC.py -n 24
fi
if [ ! -f ./BrioWu.hdf5 ]
then
    echo "Generating initial conditions for the Brio & Wu shock tube example..."
    python3 makeIC.py
fi

# Run the example with SWIFT 
../../../swift --hydro --threads=4 BrioWu.yml 2>&1 | tee output.log

# Plot the calculated solution at time t=0.1
python3 plotSolution.py BrioWu_0011.hdf5 BrioWu_0011.png	
