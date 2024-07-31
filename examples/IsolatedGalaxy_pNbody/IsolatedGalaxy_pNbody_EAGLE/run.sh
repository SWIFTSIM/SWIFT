#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

if [ ! -e galaxy.hdf5 ]
then
    echo "Fetching initial conditions to run the example..."
    wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/pNbody_galaxy.hdf5
    mv pNbody_galaxy.hdf5 galaxy.hdf5
fi

if [ ! -e coolingtables ] 
then     
    echo "Fetching EAGLE cooling tables for the isolated galaxy example..."
    ./getEagleCoolingTable.sh
fi

if [ ! -e UV_dust1_CR1_G1_shield1.hdf5 ]
then
    echo "Fetching PS2020 cooling tables for the isolated galaxy example..."
    ./getPS2020CoolingTables.sh
fi

if [ ! -e yieldtables ] 
then     
    echo "Fetching EAGLE stellar yield tables for the isolated galaxy example..."
    ./getYieldTable.sh
fi

if [ ! -e photometry ]
then
    echo "Fetching EAGLE photometry tables..."
    ./getEaglePhotometryTable.sh
fi



printf "Running simulation..."

../../swift --threads=14 --feedback --self-gravity --stars --star-formation --cooling --hydro --limiter --sync params.yml 2>&1 | tee output.log
