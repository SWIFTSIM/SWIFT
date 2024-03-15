#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./glassCube_32.hdf5 ]
then
    echo "Fetching a unit cell, copies of which are to be stacked to generate the ICs ..."
    ./getGlass.sh
fi

if [ ! -f ./LoopAdvection.hdf5 ]
then
    echo "Generating initial conditions for the Current Loop Advection example..."
    python3 makeIC.py -r 2
fi

# Run the example with SWIFT 
../../../swift --hydro --threads=4 LoopAdvection.yml 2>&1 | tee output.log

# Plot the calculated solution at time t=20.0
python3 plotSolution.py LoopAdvection_0020.hdf5 LoopAdvection_0020.png
