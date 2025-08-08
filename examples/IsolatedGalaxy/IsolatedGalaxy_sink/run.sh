#!/bin/bash

# filename=fid.hdf5
# filename=f10.hdf5
# filename=f90.hdf5
filename=lowres8.hdf5
# filename=lowres64.hdf5
# filename=lowres512.hdf5
# filename=highres6.hdf5

# make run.sh fail if a subcommand fails
set -e

if [ ! -e $filename ]
then
    echo "Fetching initial conditons for the isolated galaxy with an external potential ..."
    ./getIC.sh $filename
fi

if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    ./getChemistryTable.sh
fi


../../../swift --hydro --sinks --feedback --stars --external-gravity --self-gravity --threads=8 isolated_galaxy.yml 2>&1 | tee output.log
