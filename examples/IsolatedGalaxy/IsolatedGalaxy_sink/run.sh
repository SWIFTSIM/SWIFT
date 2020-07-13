#!/bin/bash

if [ ! -e 3e11-star-only-DM-halo-galaxy.hdf5 ]
then
    echo "Fetching initial conditons for the isolated galaxy with an external potential ..."
    ./getIC.sh
fi

echo "Generate initial conditions"
python3 makeIC.py

../../swift --sinks --external-gravity --self-gravity --threads=16 isolated_galaxy.yml 2>&1 | tee output.log
