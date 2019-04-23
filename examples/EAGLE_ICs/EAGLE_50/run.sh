#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_L0050N0752_ICs.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 50Mpc example..."
    ./getIC.sh
fi

../../swift --cosmology --hydro --self-gravity --stars --threads=16 eagle_50.yml 2>&1 | tee output.log

