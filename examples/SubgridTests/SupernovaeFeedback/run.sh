#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the Supernovae feedback example..."
    ./getGlass.sh
fi
if [ ! -e SN_feedback.hdf5 ]
then
    echo "Generating initial conditions for the Supernovae feedback example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --external-gravity --feedback --hydro --stars --threads=4 SN_feedback.yml 2>&1 | tee output.log

# Plot the solution
# TODO
