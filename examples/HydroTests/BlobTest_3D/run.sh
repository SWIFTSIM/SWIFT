#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e blob.hdf5 ]
then
    echo "Generating initial conditions for the 'blob test'..."
    python3 makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=2 blob.yml

python3 makeMovie.py
