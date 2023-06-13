#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e cylinders.hdf5 ]
then
    echo "Generating initial conditions for the square test ..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=4 cylinders.yml 2>&1 | tee output.log

# Plot the solution
python3 plotCylinders.py 100
