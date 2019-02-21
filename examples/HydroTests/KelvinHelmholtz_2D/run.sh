#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e kelvinHelmholtz.hdf5 ]
then
    echo "Generating initial conditions for the Kelvin-Helmholtz example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=4 kelvinHelmholtz.yml 2>&1 | tee output.log

# Plot the solution
python plotSolution.py 6
python makeMovie.py
