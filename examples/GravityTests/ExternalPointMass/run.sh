#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e PointMass.hdf5 ]
then
    echo "Generating initial conditions for the point mass potential box example..."
    python makeIC.py 10000
fi

rm -rf pointMass_*.hdf5
../../swift --external-gravity --threads=1 externalPointMass.yml 2>&1 | tee output.log

python energy_plot.py
