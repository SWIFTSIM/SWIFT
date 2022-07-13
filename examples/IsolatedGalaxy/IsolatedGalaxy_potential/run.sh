#!/bin/bash

if [ ! -e 3e11-star-only-static-potential-galaxy.hdf5 ]
then
    echo "Fetching initial conditons for the isolated galaxy with an external potential ..."
    ./getIC.sh
fi 

../../../swift --external-gravity --self-gravity --stars --threads=16 isolated_galaxy.yml 2>&1 | tee output.log


echo "Make plots of conservation of total angular momentum" 
if command -v python3 &>/dev/null; then
    python3 angularmomentum.py 
else
    python3 angularmomentum.py 
fi

echo "Make plots of change of vertical and radial profile"
if command -v python3 &>/dev/null; then
    python3 profilefit.py 
else
    python3 profilefit.py 
fi
