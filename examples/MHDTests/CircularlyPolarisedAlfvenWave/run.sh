#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./FCCglassCube_32.hdf5 ]
then
    echo "Generating a close packed lattice ..."
    python3 makeGlassFCC.py -n 32 -o FCCglassCube_32.hdf5
fi
if [ ! -f ./CircularlyPolarisedAlfvenWave.hdf5 ]
then
    echo "Generating initial conditions for the Circularly Polarised Alfven Wave example..."
    python3 makeIC.py
fi

# Run the example with SWIFT 
../../../swift --hydro --limiter --threads=16  CircularlyPolarisedAlfvenWave.yml 2>&1 | tee output.log
