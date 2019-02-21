#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e kelvinHelmholtzGrowthRate.hdf5 ]
then
    echo "Generating initial conditions for the Kelvin-Helmholtz growth rate " \
         "example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=1 kelvinHelmholtzGrowthRate.yml 2>&1 | tee output.log

# Plot the solution
python plotSolution.py 100
