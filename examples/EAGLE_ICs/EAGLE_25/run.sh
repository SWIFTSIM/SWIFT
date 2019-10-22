#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_L0025N0376_ICs.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 25Mpc example..."
    ./getIC.sh
fi

# Grab the cooling and yield tables if they are not present.
if [ ! -e yieldtables ]
then
    echo "Fetching EAGLE yield tables..."
    ../getEagleYieldtable.sh
fi

if [ ! -e coolingtables ]
then
    echo "Fetching EAGLE cooling tables..."
    ../getEagleCoolingTable.sh
fi

# The following run-time options are broken down by line as:
# Basic run-time options
# Create and run with stars
# Radiative options - run with cooling and stellar feedback
# Run with the time-step limiter required to capture feedback
# Run with black holes - fof is needed for the seeding
# Threading options - run with threads and pinning (latter not required but improves performance)
# The corresponding parameter file for this run

../../swift \
    --cosmology --hydro --self-gravity \
    --stars --star-formation \
    --cooling --feedback \
    --limiter \
    --black-holes --fof \
    --threads=16 --pin \
    eagle_25.yml