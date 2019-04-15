#!/bin/bash

if [ ! -e lowres8.hdf5 ] 
then     
    echo "Fetching initial conditions for the isolated galaxy example..."
    ./getIC.sh
fi

if [ ! -e cooling_tables ] 
then     
    echo "Fetching EAGLE cooling tables for the isolated galaxy example..."
    ./getEagleCoolingTable.sh
fi

../../swift --threads=16 --feedback --external-gravity --self-gravity --stars --star-formation --cooling --hydro isolated_galaxy.yml 2>&1 | tee output.log

# Kennicutt-Schmidt law plot
python3 plotSolution.py

# Plot that the random star formation matches the expected SFH
python3 SFH.py
