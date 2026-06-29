#!/bin/bash

if [ ! -e 3e11-star-only-DM-halo-galaxy.hdf5 ]
then
    echo "Downloading initial conditions"
    ./getIC.sh
fi

../../../swift -n 1 --external-gravity --stars --threads=4 isolated_galaxy.yml  2>&1 | tee output.log

python tidaltensor.py