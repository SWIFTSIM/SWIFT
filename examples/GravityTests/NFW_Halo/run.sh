#!/bin/bash

if [ ! -e test_nfw.hdf5 ]
then
    echo "Generate initial conditions for NFW example"	
    if command -v python3 &>/dev/null; then
        python3 makeIC.py
    else 
        python3 makeIC.py
    fi
fi

# self gravity G, external potential g, hydro s, threads t and high verbosity v
../../../swift --external-gravity --threads=6 test.yml 2>&1 | tee output.log

if command -v python3 &>/dev/null; then
    python3 makePlots.py
else 
    python3 makePlots.py
fi
