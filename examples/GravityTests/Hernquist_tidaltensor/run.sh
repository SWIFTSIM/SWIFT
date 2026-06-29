#!/bin/bash

# Generate IC if not present:
if [ ! -e  HernquistIC.hdf5 ]
then 
    echo "Generate initial conditions"
    python3 makeIC.py
fi


../../../swift --threads=8 --steps=1 --self-gravity --limiter --sync hernquist.yml 2>&1 | tee output.log

echo "Plot particle tidal tensor along with expected profile"

python3 tidaltensor.py