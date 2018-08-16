#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e EAGLE_ICs_50.hdf5 ]
then
    echo "Fetching initial conditions for the EAGLE 50Mpc example..."
    ./getIC.sh
fi

../swift -c -s -G -S -t 16 eagle_50.yml 2>&1 | tee output.log

