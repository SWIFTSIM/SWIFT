#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e vacuum.hdf5 ]
then
    echo "Generating initial conditions for the 1D vacuum expansion example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=1 vacuum.yml 2>&1 | tee output.log

# Plot the result
python plotSolution.py 1
