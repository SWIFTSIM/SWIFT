#!/bin/bash

if [ ! -e fid.hdf5 ] 
then     
    echo "Fetching initial conditions for the isolated galaxy example..."
    ./getIC.sh
fi

../../swift --threads=32 --external-gravity --self-gravity --stars --star-formation --cooling --temperature --hydro isolated_galaxy.yml 2>&1 | tee output.log

# Kennicutt-Schmidt law plot
python3 plotSolution.py

# Plot that the random star formation matches the expected SFH
python3 SFH.py
