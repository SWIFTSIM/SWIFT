#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e magnetised_cloud.hdf5 ]
then
    echo "Generating initial conditions for the magnetised cloud collapse example..."
    python3 makeIC.py -n 32
fi

# Run SWIFT
../../../swift --hydro --self-gravity --limiter --threads=70 magnetised_cloud.yml 2>&1 | tee output.log

# Plot the solution
python3 plotSolution.py magnetised_cloud_0070.hdf5 magnetised_cloud.png
