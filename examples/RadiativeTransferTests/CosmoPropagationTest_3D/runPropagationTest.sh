#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for Strömgen Sphere 3D example ..."
    ./getGlass.sh
fi

if [ ! -f 'propagationTest-3D.hdf5' ]; then
    echo "Generating ICs"
    python3 makePropagationTestIC.py
fi

# Run SWIFT with RT
../../../swift \
    --hydro  \
    --cosmology \
    --threads=4 \
    --stars \
    --external-gravity \
    --feedback \
    --radiation \
    ./propagationTest-3D.yml 2>&1 | tee output.log

# Plot the photon propagation checks.
# Make sure you set the correct photon group to plot
# inside the script
python3 ./plotPhotonPropagationCheck.py
