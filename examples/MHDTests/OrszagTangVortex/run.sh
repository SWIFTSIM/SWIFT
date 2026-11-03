#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f glassCube_16.hdf5 ]
then
    echo "Fetching glass cubes to be used in IC generation for the Orszag-Tang Vortex test ..."
    ./getGlass.sh
fi

if [ ! -f ./BCCglassCube_16.hdf5 ]
then
    echo "Generating a BCC unit cell, copies of which are can be stacked to generate the ICs - if you do not want to use glass files ..."
    python3 makeBCC.py -n 16
fi

if [ ! -f ./OrszagTangVortex.hdf5 ]
then
    echo "Generating initial conditions for the Orszag-Tang Vortex example..."
    python3 makeIC.py -r 4
fi

# Run the example with SWIFT 
../../../swift --hydro --threads=4 OrszagTangVortex.yml 2>&1 | tee output.log

# Plot the calculated solution at time t=1.0
python3 plotSolution.py OrszagTangVortex_0011.hdf5 OrszagTangVortex_0011.png	
