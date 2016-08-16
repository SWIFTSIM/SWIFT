#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Isothermal.hdf5 ]
then
    echo "Generating initial conditions for the disk-patch example..."
    python makeIC.py 1000
fi

../../swift -g -t 2 disk-patch.yml
