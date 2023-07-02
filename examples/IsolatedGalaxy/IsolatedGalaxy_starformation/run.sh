#!/bin/bash

if [ ! -e fid.hdf5 ] 
then     
    echo "Fetching initial conditions for the isolated galaxy example..."
    ./getIC.sh
fi

if [ ! -e coolingtables ] 
then     
    echo "Fetching EAGLE cooling tables for the isolated galaxy example..."
    ./getEagleCoolingTable.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]
then
    echo "Fetching PS2020 cooling tables for the isolated galaxy example..."
    ./getPS2020CoolingTables.sh
fi

if [ ! -e yieldtables ] 
then     
    echo "Fetching EAGLE stellar yield tables for the isolated galaxy example..."
    ./getYieldTable.sh
fi

if [ ! -e photometry ]
then
    echo "Fetching EAGLE photometry tables..."
    ../getEaglePhotometryTable.sh
fi

../../../swift --threads=16 --external-gravity --self-gravity --stars --star-formation --cooling --hydro isolated_galaxy.yml 2>&1 | tee output.log

# Kennicutt-Schmidt law plot
python3 plotSolution.py

# Plot that the random star formation matches the expected SFH
python3 SFH.py
