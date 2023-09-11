#!/bin/bash

# Get the initial conditions
if [ ! -e agora_swift.hdf5 ]
then
    echo "Fetching the initial conditions"
    ./getIC.sh
fi


# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012_shielded.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ./getGrackleCoolingTable.sh
fi


if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    ./getChemistryTable.sh
fi




../../../swift --cooling --feedback --cosmology  --limiter --sync --self-gravity --hydro --stars --star-formation --threads=24 agora_cosmo.yml 2>&1 | tee output.log



