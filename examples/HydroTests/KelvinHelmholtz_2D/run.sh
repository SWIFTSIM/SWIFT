#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e kelvinHelmholtz.hdf5 ]
then
    echo "Generating initial conditions for the Kelvin-Helmholtz example..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --threads=4 kelvinHelmholtz.yml 2>&1 | tee output.log


# Plot the solution
python3 plotSolution.py 100
python3 makeMovieSwiftsimIO.py
