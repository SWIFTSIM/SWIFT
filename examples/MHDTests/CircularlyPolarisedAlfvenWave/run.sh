#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./CPLglassCube_32x_16y_16z.hdf5 ]
then
    echo "Generating a close packed lattice ..."
    python3 makeCPL.py -nx 32 -ny 16 -nz 16
fi
if [ ! -f ./CircularlyPolarisedAlfvenWave.hdf5 ]
then
    echo "Generating initial conditions for the Circularly Polarised Alfven Wave example..."
    python3 makeIC_CPL.py -p CPLglassCube_32x_16y_16z.hdf5
fi

# Run the example with SWIFT 
../../../swift --hydro --limiter --threads=4  circularlyPolarisedAlfvenWave.yml 2>&1 | tee output.log
