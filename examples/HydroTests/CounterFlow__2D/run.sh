#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e counterFlow.hdf5 ]
then
    echo "Generating initial conditions for the counter flow example..."
    python3 makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=4 counterFlow.yml 2>&1 | tee output.log


# Plot the solution
python3 makeMovieSwiftsimIO.py
