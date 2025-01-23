#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=${n_threads:=8}  #Number of threads to use
level=${level:=5}  #Number of particles = 2^(3*level)
gas_density=${gas_density:=1} #Gas density in atom/cm^3
gas_particle_mass=${gas_particle_mass:=10000000} #Mass of the gas particles
with_hydro_MFM=${with_hydro_MFM:=0}

# Remove the ICs
# if [ -e ICs_homogeneous_box.hdf5 ]
# then
#     rm ICs_homogeneous_box.hdf5
# fi

#Create the ICs if they do not exist
if [ ! -e ICs_homogeneous_box.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --level $level --rho $gas_density --mass $gas_particle_mass -o ICs_homogeneous_box.hdf5
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
# python3 compare_num_to_sol.py snap/snapshot_*0.hdf5 --x_min 0 --x_max 1 --y_min 0 --y_max 1
python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30 --r_min 1e-1 --r_max 1.1
python3 metal_projection.py snap/snapshot_*0.hdf5 --log
