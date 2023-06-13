#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e impact.hdf5 ]
then
    echo "Generating initial conditions ..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=4 impact.yml 2>&1 | tee output.log

# Plot the solution
python3 plotImpact.py 35
