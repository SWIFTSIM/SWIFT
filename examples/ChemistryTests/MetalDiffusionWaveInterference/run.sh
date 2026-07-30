#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

scripts_location="../../GEAR_ICs_and_SCRIPTS"

# Script parameters
n_threads=${n_threads:=8}  # Number of threads to use
level=${level:=5}  # Number of particles = 2^(3*level)
gas_density=${gas_density:=1}  # Gas density in atom/cm^3
box_mass=${box_mass:=10000000} # Mass of the gas particles
vx=${vx:=0.0}  # Default velocity x-component
vy=${vy:=0.0}  # Default velocity y-component
vz=${vz:=0.0}  # Default velocity z-component
random_positions=${random_positions:=0} # Use random positions instead of regular grid?
run_name=${run_name:=""}                # Name of the run


ICs_name="metal_diffusion_wave_interference.hdf5"

# Remove the ICs
if [ -e $ICs_name ]
then
    rm $ICs_name
fi

#Create the ICs if they do not exist
if [ ! -e $ICs_name ]
then
    echo "Generating initial conditions to run the example..."
    random_flag=""

    if [ "$random_positions" -eq 1 ]; then
        random_flag="--random_positions"
    fi

    python3 makeIC.py \
        --level "$level" \
        --rho "$gas_density" \
        --mass "$box_mass" \
        --velocity "$vx" "$vy" "$vz" \
        $random_flag \
        -o $ICs_name
fi

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    $scripts_location/getGrackleCoolingTable.sh --hm2012
fi


if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    $scripts_location/getChemistryTable.sh
fi


# Create output directory
DIR=snap
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. Its content will be removed."
    rm -r $DIR
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

printf "Running simulation..."

../../../swift --hydro --external-gravity --stars --feedback \
	     --threads=$n_threads params.yml 2>&1 | tee output.log

#Do some data analysis to show what's in this box
python3 ../plot_metal_mass_conservation_in_time.py snap/*.hdf5
python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30
python3 ../metal_projection.py snap/snapshot_*0.hdf5 --log --vmin -15 --vmax -9.5
python3 ../metal_projection.py snap/snapshot_*0.hdf5 --vmin "1e-15" --vmax "1e-9"
# The two seeds are antipodal on the periodic ring (README): their fronts
# cross simultaneously at the midpoint and at the periodic seam. This is
# the exact-solution comparison that actually shows it, see the README.
python3 ../analyze_line_vs_exact.py "$(pwd)"

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
    fi
fi
