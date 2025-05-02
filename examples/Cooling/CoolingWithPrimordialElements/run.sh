#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e homogeneousCube_32.hdf5 ]
then
    echo "Fetching initial conditions..."
    ./getICS.sh
fi


# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ../getGrackleCoolingTable.sh
fi


# Create output directory
rm -rf snap
mkdir snap

# Run SWIFT
../../../swift --hydro --cooling --cosmology --threads=14  CoolingWithPrimordialElements.yml

# Check energy conservation and cooling rate
python3 plotAbundances.py snap/snapshot*.hdf5

