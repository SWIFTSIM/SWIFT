#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e cosmoVolume.hdf5 ]
then
    echo "Fetching initial conditions for the cosmo volume example..."
    ./getIC.sh
fi

../swift -s cosmoVolume.yml
