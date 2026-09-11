#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

with_sinks=${with_sinks=0} # Run with sinks (1) or with SF (0)?
with_bh=${with_bh=0}       # Seed a central black hole (1) or not (0)?

scripts_location="../../../GEAR_ICs_and_SCRIPTS"

if [ ! -e galaxy_multi_component.hdf5 ]; then
    echo "Fetching initial conditions to run the example..."
    wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/IsolatedGalaxies/galaxy_multi_component.hdf5
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
    echo "Fetching the Cloudy tables required by Grackle..."
    $scripts_location/getGrackleCoolingTable.sh
fi

if [ ! -e POPIIsw.h5 ]; then
    echo "Fetching the chemistry tables..."
    $scripts_location/getChemistryTable.sh
fi

if [[ "$with_sinks" -eq 0 ]]; then
    runtime_params="--star-formation"
else
    runtime_params="--sinks"
fi

if [[ "$with_bh" -eq 1 ]]; then
    echo "Seeding a central black hole..."
    python3 add_black_hole.py galaxy_multi_component.hdf5
    runtime_params="$runtime_params --black-holes"
fi

printf "Running simulation..."

../../../../swift --hydro --stars $runtime_params --self-gravity --feedback \
		  --cooling --threads=8 params.yml 2>&1 | tee output.log
