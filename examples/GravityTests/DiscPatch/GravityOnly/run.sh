#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Disc-Patch.hdf5 ]
then
    echo "Generating initial conditions for the disc-patch example..."
    python makeIC.py 1000
fi

../../../swift --external-gravity --threads=2 disc-patch.yml
