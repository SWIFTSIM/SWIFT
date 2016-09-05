#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Isothermal.hdf5 ]
then
    echo "Generating initial conditions for the isothermal potential box example..."
    python makeIC.py 1000 1
fi

../../swift -g -t 2 isothermal.yml 2>&1 | tee output.log
