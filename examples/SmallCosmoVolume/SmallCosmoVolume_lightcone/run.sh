#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e small_cosmo_volume.hdf5 ]
then
    echo "Fetching initial conditions for the small cosmological volume example..."
    ./getIC.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]
then
    echo "Fetching cooling tables for the small cosmological volume example..."
    ./getPS2020CoolingTables.sh
fi

if [ ! -e photometry ]
then
    echo "Fetching photometry tables for the small cosmological volume example..."
    ./getEaglePhotometryTable.sh
fi

if [ ! -e yieldtables ]
then
    echo "Fetching yield tables for the small cosmological volume example..."
    ./getEagleYieldTable.sh
fi

if [ ! -e X_Ray_tables.13072021.hdf5 ]
then
    echo "Fetching X ray tables for the small cosmological volume example..."
    ./getXrayTables.sh
fi


# Run SWIFT
../../../swift --cosmology --eagle --lightcone --pin --threads=8 \
    small_cosmo_volume.yml 2>&1 | tee output.log

