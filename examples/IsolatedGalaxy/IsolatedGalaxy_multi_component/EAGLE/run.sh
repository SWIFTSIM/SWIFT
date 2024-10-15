#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

if [ ! -e galaxy_multi_component.hdf5 ]
then
    echo "Fetching initial conditions to run the example..."
    wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/IsolatedGalaxies/galaxy_multi_component.hdf5
fi

if [ ! -e coolingtables ] 
then     
    echo "Fetching EAGLE cooling tables for the isolated galaxy example..."
    ./getEagleCoolingTable.sh
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

../../../../swift --threads=14 --feedback --self-gravity --stars --star-formation --cooling --hydro --limiter --sync params.yml 2>&1 | tee output.log
