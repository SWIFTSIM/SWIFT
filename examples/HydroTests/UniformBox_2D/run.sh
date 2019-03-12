#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e uniformPlane.hdf5 ]
then
    echo "Generating initial conditions for the uniform box example..."
    python makeIC.py 100
fi

../../swift --hydro --threads=16 uniformPlane.yml 2>&1 | tee output.log
