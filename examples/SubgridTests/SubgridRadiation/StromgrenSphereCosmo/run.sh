#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location="../../../GEAR_ICs_and_SCRIPTS"
set -eo pipefail

n_threads=${n_threads:=8}      #Number of threads to use
tier=${tier:="identity"}       #"identity" (a~=1, catches frame bugs) or "highz" (z~=9)
gas_density=${gas_density:=5}  #Target PHYSICAL gas density at a_begin, in atom/cm^3
star_mass=${star_mass:=29.7}   #Mass of the star in Msun
run_name=${run_name:=""}
restart=${restart:=0}

# boxsize/gas_mass defaults are tier-dependent: identity mirrors
# ../StromgrenSphere exactly (0.05 kpc, 0.1 Msun/particle); highz needs a
# bigger physical box for its own R_st to fit (see params_highz.yml), with
# gas_mass scaled up to keep the particle count/runtime unchanged.
case "$tier" in
    identity)
        a_begin=0.99865
        params=params_identity.yml
        ics_file=ICs_stromgren_cosmo_identity.hdf5
        default_boxsize=0.05
        default_gas_mass=0.1
        ;;
    highz)
        a_begin=0.1
        params=params_highz.yml
        ics_file=ICs_stromgren_cosmo_highz.hdf5
        default_boxsize=0.15
        default_gas_mass=2.7
        ;;
    *)
        echo "Unknown tier '$tier'; expected 'identity' or 'highz'." >&2
        exit 1
        ;;
esac
boxsize=${boxsize:=$default_boxsize}         #Target PHYSICAL boxsize at a_begin, in kpc
gas_particle_mass=${gas_mass:=$default_gas_mass} #Mass of the gas particles in Msun

# Remove and regenerate the ICs (they encode a_begin/rho/boxsize)
if [ -e "$ics_file" ]; then
    rm "$ics_file"
fi
echo "Generating initial conditions for tier=$tier (a_begin=$a_begin)..."
python3 makeIC.py --a_begin $a_begin --boxsize $boxsize --rho $gas_density \
        --mass $gas_particle_mass --star_mass $star_mass \
        -o "$ics_file"

# Get the Grackle cooling table
if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
    echo "Fetching the Cloudy tables required by Grackle..."
    "$scripts_location"/getGrackleCoolingTable.sh
fi

if [ ! -e POPIIsw.h5 ]; then
    echo "Fetching the chemistry tables..."
    "$scripts_location"/getChemistryTable.sh
fi

DIR=snap
if [ -d "$DIR" ]; then
    if [ "$restart" -ne 1 ]; then
        echo "$DIR directory exists. Its content will be removed."
        rm -r $DIR
    fi
else
    mkdir $DIR
fi

if [ "$restart" -eq 1 ]; then
    runtime_param="--restart --verbose=0"
else
    runtime_param="--verbose=0"
fi

printf "Running simulation (tier=%s)...\n" "$tier"
../../../../swift --hydro --stars --external-gravity --feedback --cooling \
                   --cosmology --sync --limiter $runtime_param \
                   --threads=$n_threads \
                   "$params" 2>&1 | tee output.log

if [ -z "$run_name" ]; then
    echo "run_name is empty."
else
    if [ -d "$run_name" ]; then
        echo "$run_name directory exists. Nothing will be moved."
    else
        mkdir -p $run_name
        mv snap $run_name
        mv output.log $run_name
        mv timesteps.txt $run_name
        mv statistics.txt $run_name
        mv unused_parameters.yml $run_name
        mv used_parameters.yml $run_name
    fi
fi
