#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_DMO_ICs_12.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE DMO 12Mpc example..."
    ./getIC.sh
fi

../../swift --cosmology --self-gravity --threads=16 eagle_12.yml 2>&1 | tee output.log

