#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e EAGLE_L0012N0188_ICs.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 12Mpc example..."
    ./getIC.sh
fi

# Grab the cooling, yield, and photometry tables if they are not present.
if [ ! -e yieldtables ]
then
    echo "Fetching EAGLE yield tables..."
    ../getEagleYieldTable.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]
then
    echo "Fetching EAGLE-XL cooling tables..."
    ../getPS2020CoolingTables.sh
fi

if [ ! -e photometry ]
then
    echo "Fetching EAGLE photometry tables..."
    ../getEaglePhotometryTable.sh
fi

# The following run-time options are broken down by line as:
# Basic run-time options
# Create and run with stars
# Radiative options - run with cooling and stellar feedback
# Run with the time-step limiter required to capture feedback
# Run with black holes - fof is needed for the seeding
# Threading options - run with threads and pinning (latter not required but improves performance)
# The corresponding parameter file for this run

../../../swift \
    --cosmology --eagle \
    --threads=16 --pin \
    eagle_12.yml

