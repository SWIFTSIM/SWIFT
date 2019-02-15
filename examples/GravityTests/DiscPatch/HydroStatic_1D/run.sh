#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Disc-Patch.hdf5 ]
then
    echo "Generating initial conditions for the disc patch example..."
    python makeIC.py
fi

# Run SWIFT
../../../swift --external-gravity --hydro --threads=4 disc-patch-icc.yml 2>&1 | tee output.log

python plotSolution.py
