#!/bin/bash

# Generate the initial conditions file if not present.
if [ ! -f ./BCCglassCube_32.hdf5 ]
then
    echo "Generating a close packed lattice ..."
    python3 makeGlassBCC.py -n 32 -o BCCglassCube_32.hdf5
fi
if [ ! -f ./CircularlyPolarisedAlfvenWave.hdf5 ]
then
    echo "Generating initial conditions for the Circularly Polarised Alfven Wave example ..."
    python3 makeIC.py
fi

# Run the example with SWIFT
../../../swift --hydro --limiter --threads=16  CircularlyPolarisedAlfvenWave.yml 2>&1 | tee output.log

# Plot solution after five wave periods
python3 plot.py CircularlyPolarisedAlfvenWave_0010.hdf5
