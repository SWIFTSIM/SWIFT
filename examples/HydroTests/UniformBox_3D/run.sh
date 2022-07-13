#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e uniformBox.hdf5 ]
then
    echo "Generating initial conditions for the uniform box example..."
    python3 makeIC.py 100
fi

../../../swift --hydro --threads=16 uniformBox.yml 2>&1 | tee output.log
