#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

if [ ! -e galaxy_multi_component.hdf5 ]
then
    echo "Fetching initial conditions to run the example..."
    wget http://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/IsolatedGalaxies/galaxy_multi_component.hdf5
fi


# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ./getGrackleCoolingTable.sh
fi


if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    ./getChemistryTable.sh
fi

printf "Running simulation..."

../../../../swift --hydro --stars --star-formation --self-gravity --feedback --cooling --threads=14 params.yml 2>&1 | tee output.log 
