#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./glassCube_32.hdf5 ]
then
    echo "Fetching a unit cell, copies of which are to be stacked to generate the ICs ..."
    ./getGlass.sh
fi

if [ ! -f ./MagneticBlastWave.hdf5 ]
then
    echo "Generating initial conditions for the Magnetic Blast Wave example..."
    python3 makeIC.py -r 4
fi

# Run the example with SWIFT 
../../../swift --hydro --threads=4 MagneticBlastWave.yml 2>&1 | tee output.log

# Plot the calculated solution at time t=0.2
python3 plotSolution.py MagneticBlastWave_0009.hdf5 MagneticBlastWave_0009.png	
