#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e small_cosmo_volume.hdf5 ]
then
    echo "Fetching initial conditions for the small cosmological volume example..."
    ./getIC.sh
fi

if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    ../../Cooling/getCoolingTable.sh 
fi

if [ ! -e coolingtables ]
then
    echo "Fetching cooling tables for the small cosmological volume example..."
    ./getCoolingTables.sh
fi

# Run SWIFT
../../swift --cosmology --hydro --self-gravity --cooling --threads=8 small_cosmo_volume.yml 2>&1 | tee output.log

# Plot the temperature evolution
python plotTempEvolution.py
