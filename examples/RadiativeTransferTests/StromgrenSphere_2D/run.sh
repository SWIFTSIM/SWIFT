#!/bin/bash

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
../../swift \
    --hydro --threads=4 --stars --external-gravity \
    --feedback --radiation --fpe \
    stromgrenSphere-2D.yml 2>&1 | tee output.log

# Plot the photon propagation checks.
# Make sure you set the correct photon group to plot
# inside the script
python3 ./plotPhotonPropagationCheck.py

