#!/bin/bash

# This example is based on the AGORA disk article (DOI: 10.3847/1538-4357/833/2/202)

# set the resolution (LOW or MED)
sim=LOW


# make run.sh fail if a subcommand fails
set -e

# Generate the initial conditions if they are not present.
if [ ! -e agora_disk.hdf5 ]
then
    echo "Fetching initial glass file for the Sedov blast example..."
    ./getIC.sh $sim
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ./getGrackleCoolingTable.sh
fi




# Run SWIFT
../../../swift --sync --limiter --cooling --hydro --self-gravity --star-formation --feedback --stars --threads=8 agora_disk.yml 2>&1 | tee output.log




