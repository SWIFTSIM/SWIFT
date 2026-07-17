#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

n_threads=${n_threads:=8}          # Number of threads to use
gas_density=${gas_density:=100}    # Diffuse background gas density in atom/cm^3
gas_particle_mass=${gas_mass:=20.0} # Mass of the gas particles (diffuse and clump);
                                    # Smith et al. 2021's own gas resolution.
                                    # See README for other levels.
star_mass=${star_mass:=19.2}       # Mass of each source in M_sun (Q_H ~ 2.5e48/s)
n_stars=${n_stars:=4}               # Number of co-located sources at the box center
star_type=${star_type:="single_star"}
density_factor=${density_factor:=100.0} # Clump density, as a multiple of gas_density
clump_distance_pc=${clump_distance_pc:=20.0} # Clump center distance from the star, in parsec
clump_radius_pc=${clump_radius_pc:=10.0}     # Clump radius, in parsec
with_cooling=${with_cooling:=1}
L=${boxsize:=0.1}                  # boxsize in kpc
n_cells=${n_cells:=3}               # must match Scheduler:max_top_level_cells in params.yml
nside=${nside:=0}                   # GEARFeedback:HII_angular_nside override: 0 (spherical) or N (12*N^2 pixels);
                                     # must fit the build's --with-number-of-hii-angular-pixels (default 12, i.e. N<=1)
rebuild_time_myr=${rebuild_time_myr:=0.01} # GEARFeedback:HII_region_rebuild_time_Myr override
deterministic=${deterministic:=0}   # GEARFeedback:HII_deterministic_boundary_ionization override
initial_metallicity=${initial_metallicity:=0} # GEARChemistry:initial_metallicity override
                                    # (Z/Zsun with scale_initial_metallicity: 1 below)
run_name=${run_name:=""}
restart=${restart:=0}
run_analysis=${run_analysis:=1}         # Run plot_temperature_hii.py on this run's final snapshot
run_matrix_analysis=${run_matrix_analysis:=0} # Also run the full-matrix erosion scripts
                                         # (clump_erosion_matrix.py, clump_hemisphere_erosion.py,
                                         # clump_radial_profile.py); only set this to 1 once every
                                         # cell of the matrix you care about has been run -- they
                                         # glob every gas${gas_mass}_z${z_label}_nside${nside}_
                                         # rebuild${rebuild_time_myr} directory in the current
                                         # directory and fail loudly if part of it is missing.
only_analysis=${only_analysis:=0}       # Skip IC generation and the simulation entirely, only run
                                         # the analysis steps below (run_analysis/run_matrix_analysis)
                                         # against whatever already exists -- e.g.
                                         # `only_analysis=1 run_matrix_analysis=1 ./run.sh` once a
                                         # sweep is done, without paying for a redundant simulation.

if [ "$only_analysis" -eq 0 ]; then

if [ -e ICs_stromgren_clump.hdf5 ]
then
    rm ICs_stromgren_clump.hdf5
fi

if [ ! -e ICs_stromgren_clump.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC_clump.py --boxsize $L --rho $gas_density \
                --mass $gas_particle_mass --star_mass $star_mass \
                --n_stars $n_stars \
                --density_factor $density_factor \
                --clump_distance_pc $clump_distance_pc \
                --clump_radius_pc $clump_radius_pc \
                --n_cells $n_cells --star_type $star_type \
            -o ICs_stromgren_clump.hdf5
fi

if [ ! -e CloudyData_UVB=HM2012_high_density.h5 ]
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
               -P GEARFeedback:HII_region_rebuild_time_Myr:$rebuild_time_myr \
               -P GEARFeedback:HII_deterministic_boundary_ionization:$deterministic \
               -P GEARChemistry:initial_metallicity:$initial_metallicity \
               params.yml 2>&1 | tee output.log
else
../../../swift --hydro --stars --external-gravity --feedback \
                --sync --limiter $runtime_param --threads=$n_threads \
               -P GEARFeedback:HII_angular_nside:$nside \
               -P GEARFeedback:HII_region_rebuild_time_Myr:$rebuild_time_myr \
               -P GEARFeedback:HII_deterministic_boundary_ionization:$deterministic \
               -P GEARChemistry:initial_metallicity:$initial_metallicity \
               params.yml 2>&1 | tee output.log
fi

fi # only_analysis

if [ "$only_analysis" -eq 0 ] && [ -n "$run_name" ]; then
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
elif [ -z "$run_name" ]; then
    echo "run_name is empty."
fi

# Where this run's snapshots actually ended up: $run_name/snap if a
# run_name is set and that move already happened (this call or a
# previous one, relevant for only_analysis=1), otherwise the plain
# top-level snap/ from this call.
output_dir="."
if [ -n "$run_name" ] && [ -d "$run_name/snap" ]; then
    output_dir="$run_name"
fi

if [ "$run_analysis" -eq 1 ]; then
    echo "Running plot_temperature_hii.py on the final snapshot..."
    python3 plot_temperature_hii.py -s "$output_dir/snap/snapshot"
fi

if [ "$run_matrix_analysis" -eq 1 ]; then
    echo "Running full-matrix erosion analysis scripts..."
    python3 clump_erosion_matrix.py
    python3 clump_hemisphere_erosion.py
    python3 clump_radial_profile.py
fi
