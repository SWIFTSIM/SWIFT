#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

#if [ ! -e galaxy_multi_component.hdf5 ]; then
#    echo "Fetching initial conditions to run the example..."
#    wget https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/IsolatedGalaxies/galaxy_multi_component.hdf5
#fi

# cleaning
rm -rf snap*

# generating initial conditions
printf "Generating ICs with 32768 parts"
./makeICs.py -n 32768

# running swift
printf "Running simulation for one step..."
../../../swift -n 1 --hydro --self-gravity --threads=14 params.yml 2>&1 | tee output.log

mv snap snap05

# generating initial conditions
printf "Generating ICs with 262144 parts"
./makeICs.py -n 262144

# running swift
printf "Running simulation for one step..."
../../../swift -n 1 --hydro --self-gravity --threads=14 params.yml 2>&1 | tee output.log

mv snap snap06


# generating initial conditions
printf "Generating ICs with 2097152 parts"
./makeICs.py -n 2097152

# running swift
printf "Running simulation for one step..."
../../../swift -n 1 --hydro --self-gravity --threads=14 params.yml 2>&1 | tee output.log

mv snap snap07

# plot
./plot.py --remove_ic_shift --nr 128 snap05/snapshot_0000.hdf5 snap06/snapshot_0000.hdf5 snap07/snapshot_0000.hdf5




