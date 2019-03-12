#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e interactingBlastWaves.hdf5 ]
then
    echo "Generating initial conditions for the Sedov blast example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=1 interactingBlastWaves.yml 2>&1 | tee output.log

# Get the high resolution reference solution if not present.
if [ ! -e interactingBlastWaves1D_exact.txt ]
then
    echo "Fetching reference solution for the interacting blast waves test..."
    ./getReference.sh
fi

# Plot the solution
python plotSolution.py 4
