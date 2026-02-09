#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

#if [ ! -e galaxy_multi_component.hdf5 ]; then
#    echo "Fetching initial conditions to run the example..."
#    wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/IsolatedGalaxies/galaxy_multi_component.hdf5
#fi


printf "Running simulation..."

../../../swift --hydro --self-gravity --threads=14 params.yml 2>&1 | tee output.log
