#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=${n_threads:=8}  #Number of threads to use
gas_density=${gas_density:=5} #Gas density in atom/cm^3
gas_particle_mass=${gas_mass:=0.1} #Mass of the gas particles
star_mass=${star_mass:=29.7} #Mass of the gas particles
star_type=${star_type:="single_star"}
with_cooling=${with_cooling:=1}
L=${boxsize:=0.05} #boxsize in kpc
run_name=${run_name:=""}
restart=${restart:=0}

# Remove the ICs
if [ -e ICs_homogeneous_box.hdf5 ]
then
    rm ICs_homogeneous_box.hdf5
fi

#Create the ICs if they do not exist
if [ ! -e ICs_homogeneous_box.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --boxsize $L --rho $gas_density \
		--mass $gas_particle_mass --star_mass $star_mass \
		--star_type $star_type \
	    -o ICs_homogeneous_box.hdf5
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
    if [ "$restart" -ne 1 ]; then
	echo "$DIR directory exists. Its content will be removed."
	rm -r $DIR
    fi
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

printf "Running simulation..."

if [ "$restart" -eq 1 ]; then
    runtime_param="--restart --verbose=0"
else
    runtime_param="--verbose=0"
fi

if [ "$with_cooling" -eq 1 ]; then
../../../swift --hydro --stars --external-gravity --feedback --cooling \
		   --sync --limiter $runtime_param --threads=$n_threads \
	       params.yml 2>&1 | tee output.log
else
../../../swift --hydro --stars --external-gravity --feedback \
		--sync --limiter $runtime_param --threads=$n_threads \
	       params.yml 2>&1 | tee output.log
fi

#Do some data analysis to show what's in this box

# Gas density projection
# python3 plot_gas_density.py -i 0 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 1 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 2 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 3 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 4 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 5 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 6 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 7 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 8 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 9 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 10 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 20 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 30 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 40 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 50 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 100 -s 'snap/snapshot'

# # Phase space diagram
# python3 rhoTPlot.py -i 0 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 1 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 2 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 3 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 4 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 5 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 6 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 7 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 8 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 9 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 10 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 20 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 30 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 40 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 50 -s 'snap/snapshot'
# python3 rhoTPlot.py -i 100 -s 'snap/snapshot'


# # Internal energy projection
# python3 plot_projected_qty.py --qty "internal_energy" snap/snapshot_000*.hdf5 --log --vmin 2 --vmax 6
# python3 plot_projected_qty.py --qty "internal_energy" snap/snapshot_*0.hdf5 --log --vmin 2 --vmax 6

# # Movies
# # python3 rhoTPlot.py -i 0 -f 100 -s 'snap/snapshot'
# python3 plot_gas_density.py -i 0 -f 100 -s 'snap/snapshot'


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
