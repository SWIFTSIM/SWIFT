#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=${n_threads:=8}          # Number of threads to use
gas_density=${gas_density:=100}    # Diffuse background gas density in atom/cm^3
gas_particle_mass=${gas_mass:=0.001} # Mass of the gas particles (diffuse and clump)
star_mass=${star_mass:=29.7}       # Mass of the star
star_type=${star_type:="single_star"}
density_factor=${density_factor:=100.0} # Clump density, as a multiple of gas_density
clump_radius_pc=${clump_radius_pc:=""}  # Optional: clump radius in parsec, for a size sweep
                                         # (see README); empty = use the default half-angle size
with_cooling=${with_cooling:=1}
L=${boxsize:=0.004}                # boxsize in kpc
n_cells=${n_cells:=3}               # must match Scheduler:max_top_level_cells in params.yml
nside=${nside:=0}                   # GEARFeedback:HII_angular_nside override: 0 (spherical) or 1 (12 pixels)
run_name=${run_name:=""}
restart=${restart:=0}

if [ -e ICs_stromgren_clump.hdf5 ]
then
    rm ICs_stromgren_clump.hdf5
fi

if [ ! -e ICs_stromgren_clump.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    clump_radius_arg=""
    if [ -n "$clump_radius_pc" ]; then
        clump_radius_arg="--clump_radius_pc $clump_radius_pc"
    fi
    python3 makeIC_clump.py --boxsize $L --rho $gas_density \
                --mass $gas_particle_mass --star_mass $star_mass \
                --density_factor $density_factor $clump_radius_arg \
                --n_cells $n_cells --star_type $star_type \
            -o ICs_stromgren_clump.hdf5
fi

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

DIR=snap
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
               -P GEARFeedback:HII_angular_nside:$nside \
               params.yml 2>&1 | tee output.log
else
../../../swift --hydro --stars --external-gravity --feedback \
                --sync --limiter $runtime_param --threads=$n_threads \
               -P GEARFeedback:HII_angular_nside:$nside \
               params.yml 2>&1 | tee output.log
fi

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
    fi
fi
