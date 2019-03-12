#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e small_cosmo_volume.hdf5 ]
then
    echo "Fetching initial conditions for the small cosmological volume example..."
    ./getIC.sh
fi

# Run SWIFT
../../swift --cosmology --self-gravity --threads=8 small_cosmo_volume_dm.yml 2>&1 | tee output.log

