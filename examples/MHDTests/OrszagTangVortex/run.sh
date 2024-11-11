#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -e glassCube_16.hdf5 ]
then
    echo "Generating initial conditions for the magnetised cloud collapse example..."
    ./getGlass.sh
fi

if [ ! -f ./BCCglassCube_16.hdf5 ]
then
    echo "Generating a BCC unit cell, copies of which are to be stacked to generate the ICs ..."
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
