#!/bin/bash

if [ ! -e circularorbitshernquist.hdf5 ]
then 
    echo "Generate initial conditions for circular orbits"
    if command -v python3 &>/dev/null; then
        python3 makeIC.py
    else 
        python makeIC.py
    fi

fi

# self gravity G, external potential g, hydro s, threads t and high verbosity v
../../swift --external-gravity --threads=6 hernquistcirc.yml 2>&1 | tee output.log


echo "Save plots of the circular orbits"
if command -v python3 &>/dev/null; then
    python3 plotprog.py
else 
    python plotprog.py
fi
