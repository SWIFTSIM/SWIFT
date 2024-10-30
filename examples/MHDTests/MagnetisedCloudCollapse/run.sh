#!/bin/bash

if [ ! -e glassCube_128.hdf5 ]
then
    echo "Generating initial conditions for the magnetised cloud collapse example..."
    ./getGlass.sh 
fi

# Generate the initial conditions if they are not present.
if [ ! -e magnetised_cloud.hdf5 ]
then
    echo "Generating initial conditions for the magnetised cloud collapse example..."
    python3 makeIC.py
fi

# Run SWIFT
../../../swift --hydro --self-gravity --limiter --threads=16 magnetised_cloud.yml 2>&1 | tee output.log
