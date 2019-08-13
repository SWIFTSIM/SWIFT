#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_L0012N0188_ICs.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 12Mpc example..."
    ./getIC.sh
fi

../../swift --cosmology --hydro --self-gravity --stars --black-holes --cooling --star-formation --feedback --fof --threads=16 eagle_12.yml 2>&1 | tee output.log

