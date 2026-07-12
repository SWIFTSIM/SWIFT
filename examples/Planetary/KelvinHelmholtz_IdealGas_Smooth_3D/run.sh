#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e kelvin_helmholtz.hdf5 ]
then
    echo "Generating initial conditions for the 3D Kelvin--Helmholtz test..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=4 kelvin_helmholtz.yml 2>&1 | tee output_kelvin_helmholtz.log

# Plot the solutions
python3 ./plotSnapshots.py
python3 ./plotModeGrowth.py
