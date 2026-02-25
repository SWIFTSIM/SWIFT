#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

# cleaning
rm -rf snap*

# generating initial conditions
printf "Generating ICs"
./makeICs.py -n 1024

# running swift
printf "Running simulation for one step..."
../../../swift -n 1 --hydro --self-gravity --threads=1 params.yml 2>&1 | tee output.log


# plot result
./plot.py snap/snapshot_0000.hdf5







