#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e rayleigh_taylor.hdf5 ]
then
    echo "Generating initial conditions for the 3D Rayleigh--Taylor test..."
    python3 makeIC.py
fi

# Download the equation of state tables if not already present
if [ ! -e ../EoSTables/CD21_HHe.txt ]
then
    cd ../EoSTables
    ./get_eos_tables.sh
    cd -
fi

# Run SWIFT
../../../swift --hydro --external-gravity --threads=4 rayleigh_taylor.yml 2>&1 | tee output_rayleigh_taylor.log

# Plot the solutions
python3 ./plotSnapshots.py
