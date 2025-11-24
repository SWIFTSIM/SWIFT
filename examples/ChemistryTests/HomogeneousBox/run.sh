#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=8 #Number of threads to use
level=5 #Number of particles = 2^(3*level)
jeans_length=${jeans_length:=0.250}  #Jeans wavelenght in unit of the boxsize
rho=${rho:=0.1} # Gas density
run_name=${run_name:=""}
with_cosmo=${with_cosmo:=0}


# Remove the ICs
if [ -e ICs_homogeneous_box.hdf5 ]
then
    rm ICs_homogeneous_box.hdf5
fi

#Create the ICs if they do not exist
if [ ! -e ICs_homogeneous_box.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --level $level -o ICs_homogeneous_box.hdf5 \
	    --lJ $jeans_length --rho $rho
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

printf "Running simulation..."

if [ "$with_cosmo" -eq 1 ]; then
../../../swift --hydro --star-formation --stars --self-gravity --feedback \
	       --cooling --sync --limiter --cosmology \
	       --threads=$n_threads params_cosmology.yml 2>&1 | tee output.log
else
../../../swift --hydro --star-formation --stars --self-gravity --feedback \
	       --cooling --sync --limiter \
	       --threads=$n_threads params.yml 2>&1 | tee output.log
fi

#Do some data analysis to show what's in this box
python3 plot_metal_mass_conservation_in_time.py snap/*.hdf5 --symlog
python3 metal_projection.py snap/snapshot_*0.hdf5 --vmin=-20 --vmax=-3 --log
python3 plot_projected_qty.py --qty "internal_energy" snap/snapshot_*0.hdf5 --log
python3 plot_projected_qty.py --qty "pressure" snap/snapshot_*0.hdf5 --log

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
