#!/bin/bash


# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../GEAR_ICs_and_SCRIPTS"
# make run.sh fail if a subcommand fails. pipefail matters here specifically:
# swift's own exit code is piped into `tee output.log`, and without it a
# crashed/errored swift run still lets the pipeline "succeed" (tee's own
# exit code), silently producing a run_name output directory that looks
# complete but only contains a startup-error log.
set -eo pipefail

n_threads=${n_threads:=8}     # Number of threads to use
gas_density=${gas_density:=5} # Gas density in atom/cm^3
gas_particle_mass=${gas_mass:=1.0} # Mass of the gas particles (Msun); see README
star_mass=${star_mass:=29.7}  # Mass of each star (Msun)
star_type=${star_type:="single_star"}
with_cooling=${with_cooling:=1}
L=${boxsize:=0.24}            # boxsize in kpc
n_cells=${n_cells:=3}         # must match Scheduler:max_top_level_cells in params.yml; must stay 3, see README
run_name=${run_name:=""}
restart=${restart:=0}

if [ -e ICs_stromgren_cluster.hdf5 ]
then
    rm ICs_stromgren_cluster.hdf5
fi

if [ ! -e ICs_stromgren_cluster.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC_cluster.py --boxsize $L --rho $gas_density \
                --mass $gas_particle_mass --star_mass $star_mass \
                --n_cells $n_cells --star_type $star_type \
            -o ICs_stromgren_cluster.hdf5
fi

if [ ! -e CloudyData_UVB=HM2012.h5 ]
then
    echo "Fetching the Cloudy tables required by Grackle..."
    "$scripts_location"/getGrackleCoolingTable.sh
fi

if [ ! -e POPIIsw.h5 ]
then
    echo "Fetching the chemistry tables..."
    "$scripts_location"/getChemistryTable.sh
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
../../../../swift --hydro --stars --external-gravity --feedback --cooling \
               --sync --limiter $runtime_param --threads=$n_threads params.yml 2>&1 | tee output.log
else
../../../../swift --hydro --stars --external-gravity --feedback \
                --sync --limiter $runtime_param --threads=$n_threads \
               params.yml 2>&1 | tee output.log
fi

# Plot the final mass-weighted temperature projection: the Checking section
# above's pass criterion (all 7 regions merged, no un-ionized gap facing the
# center star) is a qualitative, visual check, not a pass/fail script.
last_snap_file=$(ls snap/snapshot_*.hdf5 2>/dev/null | sort | tail -n 1)
if [ -n "$last_snap_file" ]; then
    snap_num=$(echo "$last_snap_file" | grep -oE '[0-9]{4}')
    last_snap=$((10#$snap_num))
    python3 plot_temperature_map.py -s snap/snapshot -i $last_snap
else
    echo "No snapshots found in snap/ directory! Skipping plot_temperature_map.py."
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
        if ls *.png >/dev/null 2>&1; then
            mv *.png $run_name
        fi
    fi
fi
