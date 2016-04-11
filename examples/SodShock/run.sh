#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e sodShock.hdf5 ]
then
    echo "Generating initial conditions for the SodShock example..."
    python makeIC.py
fi

../swift -s -t 16 sodShock.yml
