#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

# Script parameters
n_threads=${n_threads:=8}  #Number of threads to use
level=${level:=5}  #Number of particles = 2^(3*level)
gas_density=${gas_density:=1} #Gas density in atom/cm^3
box_mass=${box_mass:=10000000} #Mass of the gas particles
vx=${vx:=0.0}  # Default velocity x-component
vy=${vy:=0.0}  # Default velocity y-component
vz=${vz:=0.0}  # Default velocity z-component
ZR=${ZR=0.0001}  # Metalicity on the right on the discontinuity
ZL=${ZL=0}       # Metalicity on the left on the discontinuity
rhoR=${rhoR=1}   # Density on the right on the discontinuity
rhoL=${rhoL=1}   # Density on the rightleft on the discontinuity
with_shear=${with_shear:=0} # Add a velocity shearing effect
with_hydro_MFM=${with_hydro_MFM:=0}
run_name=${run_name:=""}

# Remove the ICs
if [ -e ICs_homogeneous_box.hdf5 ]
then
    rm ICs_homogeneous_box.hdf5
fi

#Create the ICs if they do not exist
if [ ! -e ICs_homogeneous_box.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    if [ "$with_shear" -eq 0 ]; then
	python3 makeIC.py \
		--level "$level" \
		--rho "$gas_density" \
		--mass "$box_mass" \
		--velocity "$vx" "$vy" "$vz" \
		--Z_R "$ZR" \
		--Z_L "$ZL" \
		--rho_R "$rhoR" \
		--rho_L "$rhoL" \
		-o ICs_homogeneous_box.hdf5
    else
	python3 makeIC.py \
		--level "$level" \
		--rho "$gas_density" \
		--mass "$box_mass" \
		--velocity "$vx" "$vy" "$vz" \
		--Z_R "$ZR" \
		--Z_L "$ZL" \
		--rho_R "$rhoR" \
		--rho_L "$rhoL" \
		--add_shear \
		-o ICs_homogeneous_box.hdf5
    fi
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

# Run the appropriate command based on with_hydro_MFM value
if [ "$with_hydro_MFM" -eq 1 ]; then
    # ./configure --with-hydro=gizmo-mfm --with-chemistry=GEAR-MFM-DIFFUSION_10 --with-stars=GEAR --with-kernel=wendland-C2 --with-grackle=$GRACKLE_ROOT --with-tbbmalloc --enable-compiler-warnings --enable-debug --enable-debugging-checks --with-riemann-solver=hllc && make clean && make -j12
    echo "Running with MFM hydro solver..."
    ~/swiftsim/swift --hydro --external-gravity --stars --threads=$n_threads params.yml 2>&1 | tee output.log
else
    # ./configure --with-chemistry=GEAR-MFM-DIFFUSION_10 --with-cooling=grackle_0 --with-stars=GEAR --with-star-formation=GEAR --with-feedback=GEAR --with-sink=GEAR --with-kernel=wendland-C2 --with-adaptive-softening=yes --with-grackle=$GRACKLE_ROOT --with-tbbmalloc --enable-compiler-warnings --enable-debug --enable-debugging-checks
    echo "Running with SPH hydro solver"
    ~/swiftsim/swift --hydro --external-gravity --stars --feedback --threads=$n_threads params.yml 2>&1 | tee output.log
fi

#Do some data analysis to show what's in this box
# python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30 --r_min 1e-1
# --r_max=1.1
python3 plot_metal_mass_conservation_in_time.py snap/*.hdf5
python3 metal_projection.py snap/snapshot_*0.hdf5 --log
