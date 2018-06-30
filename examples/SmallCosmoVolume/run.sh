#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e small_cosmo_volume.hdf5 ]
then
    echo "Fetching initial conditions for the small cosmological volume example..."
    ./getIC.sh
fi

../swift -a -c -G -t 16 small_cosmo_volume.yml 2>&1 | tee output.log

