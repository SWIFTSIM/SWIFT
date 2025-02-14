#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=8 #Number of threads to use
level=5 #Number of particles = 2^(3*level)
jeans_length=0.250 #Jeans wavelenght in unit of the boxsize

#Create the ICs if they do not exist
if [ ! -e ICs_homogeneous_box.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --level $level -o ICs_homogeneous_box.hdf5 --lJ $jeans_length
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

../../../swift --hydro --star-formation --stars --self-gravity --feedback \
	       --cooling --sync --limiter \
	       --threads=$n_threads params.yml 2>&1 | tee output.log

# ../../../swift --hydro --stars --self-gravity --cooling --feedback \
	       # --sync --limiter \
	       # --threads=$n_threads params.yml 2>&1 | tee output.log

#Do some data analysis to show what's in this box
python3 plot_metal_mass_conservation_in_time.py snap/*.hdf5 --symlog
python3 metal_projection.py snap/snapshot_*0.hdf5 --vmin=-20 --vmax=-3 --log
python3 plot_projected_qty.py --qty "internal_energy" snap/snapshot_*0.hdf5 --log
python3 plot_projected_qty.py --qty "pressure" snap/snapshot_*0.hdf5 --log
