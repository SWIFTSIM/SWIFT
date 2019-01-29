#!/bin/bash

if [ ! -e fid.hdf5 ] 
then     
    echo "Fetching initial conditions for the isolated galaxy example..."
    ./getIC.sh
fi

../swift --threads=32 --external-gravity --self-gravity --stars --star-formation --cooling --temperature --hydro isolated_galaxy.yml 2>&1 | tee output.log

python3 plotSolution.py
