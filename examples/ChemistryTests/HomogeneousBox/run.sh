#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=${n_threads:=8}  #Number of threads to use
level=${level:=5}  #Number of particles = 2^(3*level)
jeans_length=${jeans_length:=0.250}  #Jeans wavelenght in unit of the boxsize
rho=${rho:=0.1} # Gas density
run_name=${run_name:=""}
with_cosmo=${with_cosmo:=0}
with_initial_metals=${with_initial_metals:=0} # Do we want some initial pollution?
with_restart=${with_restart:=0}

# Choose the right parameters
if [ "$with_initial_metals" -eq 1 ];
then
	echo "Running with an initial metal pollution (no SF and feedback)"

	# With initial pollution: no SF
	metal_mass_fraction=1e-4
    radius=0.05
	swift_run_parameters=""

	if [ "$with_cosmo" -eq 1 ];
	then
		parameter_file="params_cosmology_initial_metallicity.yml"
	else
		parameter_file="params_initial_metallicity.yml"
	fi
else
	echo "Running with SF and feedback (no initial pollution)"

	# Without initial pollution: do SF
	metal_mass_fraction=0.0
	radius=0.0
	swift_run_parameters="--star-formation --sync"

	if [ "$with_cosmo" -eq 1 ];
	then
		parameter_file="params_cosmology.yml"
	else
		parameter_file="params.yml"
	fi
fi

# Add cosmology runtime parameter to the executable
if [ "$with_cosmo" -eq 1 ];
then
	echo "Running with cosmology"
	swift_run_parameters="$swift_run_parameters --cosmology"
fi

# Add restart runtime parameter to the executable
if [ "$with_restart" -eq 1 ]; then
	echo "Running with cosmology"
	swift_run_parameters="$swift_run_parameters --restart"
fi

echo "Parameter file: $parameter_file"

# Remove the ICs
if [ "$with_restart" -eq 0 ]; then
	if [ -e ICs_homogeneous_box.hdf5 ]
	then
		echo "------------------------------"
		echo "Removing ICs..."
		rm ICs_homogeneous_box.hdf5
	fi

	#Create the ICs if they do not exist
	if [ ! -e ICs_homogeneous_box.hdf5 ]
	then
		echo "Generating initial conditions to run the example..."
		python3 makeIC.py --level $level -o ICs_homogeneous_box.hdf5 \
			--lJ $jeans_length --rho $rho --epsilon $radius --metal_mass_fraction $metal_mass_fraction
	fi
fi


# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    ./getGrackleCoolingTable.sh
fi

# Get the chemistry table
if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    ./getChemistryTable.sh
fi


# Create output directory
DIR=snap #First test of units conversion

if [ "$with_restart" -eq 0 ]; then
	if [ -d "$DIR" ];
	then
		echo "$DIR directory exists. Its content will be removed."
		rm -r $DIR
	else
		echo "$DIR directory does not exists. It will be created."
		mkdir $DIR
	fi
fi

printf "Running simulation..."


../../../swift --hydro --stars --self-gravity --feedback \
	       --cooling --limiter $swift_run_parameters --pin  \
	       --threads=$n_threads $parameter_file 2>&1 | tee output.log

#Do some data analysis to show what's in this box
python3 plot_metal_mass_conservation_in_time.py snap/*.hdf5 #--symlog
python3 metal_projection.py snap/snapshot_*0.hdf5 --vmin=-13 --vmax=-8 --log
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
