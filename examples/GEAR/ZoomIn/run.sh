#!/bin/bash

if [ ! -e h050.hdf5 ]
then 
    echo "Fetching initial conditions for the zoom in example..."
    ./getIC.sh
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







../../../swift --cooling --feedback --cosmology  --limiter --sync --self-gravity --hydro --stars --star-formation --threads=8 zoom_in.yml 2>&1 | tee output.log
