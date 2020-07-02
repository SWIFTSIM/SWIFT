#!/bin/bash

if [ ! -e Isolated_NFW_MN.hdf5 ] 
then     
    echo "Fetching initial conditions for the isolated galaxy example..."
    ./getIC.sh
fi

../../swift -g -G -s --threads=16 isolated_galaxy.yml
