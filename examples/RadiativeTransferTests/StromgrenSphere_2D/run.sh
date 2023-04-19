#!/bin/bash

#---------------------------------------
# Runs the Stromgren Sphere example
#---------------------------------------

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -e glassPlane_128.hdf5 ]
then
    echo "Fetching initial glass file for Strömgen Sphere 2D example ..."
    ./getGlass.sh
fi

if [ ! -f 'stromgrenSphere-2D.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../../swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation \
    stromgrenSphere-2D.yml 2>&1 | tee output.log

# Plot the Stromgren 2D checks.
python3 ./plotSolution.py

