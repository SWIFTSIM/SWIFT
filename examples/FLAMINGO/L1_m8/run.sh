#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e FLAMINGO_m8_ICs.hdf5 ]
then
    echo "Fetching initial conditions for the FLAMIGNO 100Mpc example..."
    ./getIC.sh
fi

# Grab the cooling, yield, X-ray, and photometry tables if they are not present.
if [ ! -e yieldtables ]
then
    echo "Fetching EAGLE yield tables..."
    ../getEagleYieldTable.sh
fi

if [ ! -e X_Ray_tables_flamingo.hdf5 ]
then
    echo "Fetching FLANINGO X-ray tables..."
    ../getXrayTables.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]
then
    echo "Fetching PS-2020 cooling tables..."
    ../getPS2020CoolingTables.sh
fi

if [ ! -e photometry ]
then
    echo "Fetching EAGLE photometry tables..."
    ../getEaglePhotometryTable.sh
fi


# Run the code
../../../swift --flamingo -t 16 flamingo.yml --cosmology  --power --lightcone -v 1
