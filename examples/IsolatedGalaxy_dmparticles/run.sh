#!/bin/bash

if [ ! -e reddeathgalaxywithDM.hdf5 ]
then 
    echo "Fetch the isolated galaxy initial conditions with a DM halo"
    ./getIC.sh
fi

../swift -G -S -t 64 isolated_galaxy.yml 2>&1 | tee output.log

