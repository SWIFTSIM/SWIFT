#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_ICs_50.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 50Mpc example..."
    ./getIC.sh
fi

../../swift --cosmology --hydro --self-gravity --stars --threads=16 eagle_50.yml 2>&1 | tee output.log

