#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e KelvinHelmholtzMHD.hdf5 ]
then
    echo "Generating initial conditions for the Kelvin-Helmholtz growth rate " \
         "example..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=8 KelvinHelmholtzMHD.yml 2>&1 | tee output.log

# Plot the solution
python3 plotSolution.py 100
