#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the Gresho-Chan vortex example..."
    ../../HydroTests/GreshoVortex_3D/getGlass.sh
fi
if [ ! -e greshoVortex.hdf5 ]
then
    echo "Generating initial conditions for the Gresho-Chan vortex example..."
    python3 makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=4 gresho.yml 2>&1 | tee output.log

# Plot the solution
python3 ../../HydroTests/GreshoVortex_3D/plotSolution.py 11
