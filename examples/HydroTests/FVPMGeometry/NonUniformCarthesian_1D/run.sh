#!/bin/bash

n_threads=${n_threads:=2}  # Number of threads to use
level=${level:=5}  # Number of particles = 2^(dimension*level)
uniform_mode=${uniform:=0}
dimension=${dim:=1} # Number of dimensions of the problem.
run_name=${run_name:=""}

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
    python3 makeIC.py --dimension $dimension --level $level $carthesian_flag
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

if [ -z "$run_name" ]; then
    echo "run_name is empty."
else
    if [ -d "$run_name" ]; then
	echo "$run_name directory exists. Nothing will be moved."
    else
	echo "$run_name directory does not exists. It will be created."
	mkdir $run_name
	mv timesteps.txt $run_name
	mv snap $run_name
	mv unused_parameters.yml $run_name
	mv used_parameters.yml $run_name
	mv *.png $run_name
	mv output.log $run_name
	mv statistics.txt $run_name
    fi
fi

