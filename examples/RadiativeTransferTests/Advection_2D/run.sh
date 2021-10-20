#!/bin/bash

# exit if anything fails
set -e
set -o pipefail

 # Generate the initial conditions if they are not present.
if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for the 2D RT advection example..."
    ./getGlass.sh
fi
if [ ! -e advection_2D.hdf5 ]
then
    echo "Generating initial conditions for the 2D RT advection example..."
    python makeIC.py
fi

# Run SWIFT with RT
../../swift \
    --hydro \
    --threads=4 \
    --verbose=0  \
    --radiation \
    --stars \
    --feedback \
    --external-gravity \
    -e \
    ./rt_advection2D.yml 2>&1 | tee output.log

python3 ./plotSolution.py
