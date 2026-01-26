#!/bin/bash

n_threads=${n_threads:=2}  # Number of threads to use
level=${level:=5}  # Number of particles = 2^(dimension*level)
uniform_mode=${uniform:=0}

if [ "$uniform_mode" -eq 1 ]; then
    echo "Mode: Pure Cartesian"
    carthesian_flag="--cartesian"
else
    echo "Mode: Dual Resolution"
    carthesian_flag=""
fi

# Remove the ICs
if [ -e advection.hdf5 ]
then
    rm advection.hdf5
fi

 # Generate the initial conditions if they are not present.
if [ ! -e advection.hdf5 ]
then
    echo "Generating initial conditions for the NonUniformrCarthesian_1D example..."
    python3 makeIC.py --dimension 1 --level $level $carthesian_flag
fi

# Create output directory
DIR="snap"
if [ -d "$DIR" ]; then
    echo "$DIR directory exists. Its content will be removed."
    rm -r "$DIR"
fi
mkdir "$DIR"

# Run SWIFT
../../../../swift --hydro --threads=$n_threads params.yml 2>&1 | tee output.log

# TODO: Add a python script
python3 plot_solution.py snap/snapshot_0000.hdf5 snap/snapshot_0050.hdf5
