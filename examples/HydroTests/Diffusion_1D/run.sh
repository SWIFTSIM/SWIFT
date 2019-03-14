#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e diffusion.hdf5 ]
then
    echo "Generating initial conditions for the Sedov blast example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --limiter --threads=1 diffusion.yml 2>&1 | tee output.log
