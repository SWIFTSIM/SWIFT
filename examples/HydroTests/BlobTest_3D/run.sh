#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e blob.hdf5 ]
then
    echo "Generating initial conditions for the Sod shock example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=4 blob.yml

python makeMovie.py
