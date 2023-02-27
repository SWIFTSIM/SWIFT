#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e turb-box.hdf5 ]
then
    echo "Generating initial conditions for the Sedov blast example..."
    ### DO ### python3 makeIC.py
fi

# Run SWIFT
./swift --hydro --threads=4 turb-box.yml 2>&1 | tee output.log

# Plot the solution
### DO ###python3 plotSolution.py 5
