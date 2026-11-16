#!/bin/bash

set -eo pipefail

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_16.hdf5 ]
then
    echo "Fetching initial glass file for the pressure floor example..."
    ./getGlass.sh
fi
if [ ! -e pressureFloor.hdf5 ]
then
    echo "Generating initial conditions for the pressure floor example..."
    python3 makeIC.py
fi

rm -f pressureFloor_*
# Run SWIFT
../../../swift --self-gravity --hydro --cooling --threads=8 pressureFloor.yml

# Check if the simulation collapsed
python3 plotDensity.py 80
