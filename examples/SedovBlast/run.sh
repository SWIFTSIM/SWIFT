#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e sedov.hdf5 ]
then
    echo "Generating initial conditions for the SedovBlast example..."
    python makeIC_fcc.py
fi

../swift -s sedov.yml
