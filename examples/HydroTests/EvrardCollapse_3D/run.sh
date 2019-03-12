#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e evrard.hdf5 ]
then
    echo "Generating initial conditions for the Evrard collapse example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --self-gravity --threads=4 evrard.yml 2>&1 | tee output.log

# Get the high resolution 1D reference result if not present.
if [ ! -e evrardCollapse3D_exact.txt ]
then
    echo "Fetching the reference result for the Evrard collapse example..."
    ./getReference.sh
fi

# Plot the solution
python plotSolution.py 8
