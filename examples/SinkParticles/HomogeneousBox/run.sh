#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_ranks=${n_ranks:=0}      # Number of ranks to use
n_threads=${n_threads:=8}  # Number of threads to use
level=${level:=5}  # Number of particles = 2^(3*level)
jeans_length=${jeans_length:=0.250} # Jeans wavelength in unit of the boxsize
debug=${debug:=0}
run_name=${run_name:=""}
with_star_formation=${with_star_formation=0}

scripts_location="../../GEAR_ICs_and_SCRIPTS"


echo "========================================"
echo "Preparing the simulation..."
echo "========================================"

# Remove the ICs
if [ -e ICs_homogeneous_box.hdf5 ]
then
    rm ICs_homogeneous_box.hdf5
fi

#Create the ICs if they do not exist
if [ ! -e ICs_homogeneous_box.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --level $level -o ICs_homogeneous_box.hdf5 --lJ $jeans_length
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012_high_density.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    $scripts_location/getGrackleCoolingTable.sh
fi

if [ ! -e POPII.h5 ]
then
    echo "Fetching the chemistry tables..."
    $scripts_location/getChemistryTable.sh --with-winds
fi

# Get the debugging ICs
if [[ -z "$debug" || "$debug" -eq 1 ]]; then
    echo "Preparing simulation for debug mode..."

    if [ ! -e snapshot_0003restart.hdf5 ]
    then
	echo "Fetching the debugging ICs..."
	./getDebuggingICs.sh
    fi

    parameter_file="params_debug.yml"
else
    parameter_file="params.yml"
fi

# Create output directory
DIR=snap #First test of units conversion
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. Its content will be removed."
    rm -r $DIR
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

if [[ $n_ranks -gt 0 ]]; then
  swift="mpirun -n $n_ranks ../../../swift_mpi"
else
  swift="../../../swift"
fi

if [[ "$with_star_formation" -eq 0 ]]; then
    runtime_params="--sinks"
else
    runtime_params="--star-formation"
fi

echo "========================================"
echo "Running simulation..."
echo "========================================"
$swift --hydro $runtime_params --stars --self-gravity --feedback \
       --cooling --sync --limiter --threads=$n_threads \
       $parameter_file 2>&1 | tee output.log


echo "========================================"
echo "Perfoming data analysis..."
echo "========================================"

#Do some data analysis to show what's in this box
python3 $scripts_location/plot_gas_density.py -i 282 -s 'snap/snapshot'
python3 $scripts_location/rhoTPlot.py -i 282 -s 'snap/snapshot'
python3 $scripts_location/rhoTPlot.py -i 0 -f 282 -s 'snap/snapshot'
python3 $scripts_location/plot_gas_density.py -i 0 -f 282 -s 'snap/snapshot'

echo "========================================"
echo "Saving the results..."
echo "========================================"
# Now save the results if you wanted to
if [ -z "$run_name" ]; then
    echo "run_name is empty."
else
    if [ -d "$run_name" ]; then
	echo "$run_name directory exists. Nothing will be moved."
    else
	echo "$run_name directory does not exists. It will be created."
	mkdir -p $run_name
	mv snap $run_name
	mv output.log $run_name
	mv timesteps.txt $run_name
	mv statistics.txt $run_name
	mv unused_parameters.yml $run_name
	mv used_parameters.yml $run_name
	mv *.png $run_name
	mv *.mp4 $run_name
    fi
fi
