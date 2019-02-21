#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_DMO_ICs_25.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE DMO 25Mpc example..."
    ./getIC.sh
fi

../../swift --cosmology --self-gravity --threads=16 eagle_25.yml 2>&1 | tee output.log

