#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_32.hdf5 ]
then
    echo "Fetching initial glass file for the cooling box example..."
    ./getGlass.sh
fi
if [ ! -e coolingBox.hdf5 ]
then
    echo "Generating initial conditions for the cooling box example..."
    python3 makeIC.py
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ../getGrackleCoolingTable.sh
fi

# Get the results
if [ ! -e CoolingHeatingBox_results.txt]
then
    echo "Fetching the the results..."
    ./getResults.sh
fi


# Create output directory
rm -rf snap
mkdir snap

# Run SWIFT
../../../swift --hydro --cooling --threads=14  coolingBox.yml

# Check energy conservation and cooling rate
python3 plotResults.py
