#!/bin/bash


# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../GEAR_ICs_and_SCRIPTS"
# make run.sh fail if a subcommand fails. pipefail matters here specifically:
# swift's own exit code is piped into `tee output.log`, and without it a
# crashed/errored swift run still lets the pipeline "succeed" (tee's own
# exit code), silently producing a run_name output directory that looks
# complete but only contains a startup-error log.
set -eo pipefail

# Defaults reproduce Hu et al. (2017, MNRAS 471, 2151) Appendix A's D-type
# expansion convergence test: 4 co-located 19.2 Msun sources (Q_H~2.5e48/s
# each, 1e49/s combined), n_H=100/cm3 uniform pure-hydrogen medium,
# T_o=1e3 K, run to t=8 Myr. gas_mass=4.0 Msun is their own coarsest tested
# resolution (they also tested 0.5, 0.0625, 0.0078125 Msun -- see README
# for the resolution/runtime tradeoff before going finer). No clump --
# reuses StromgrenSphereClumpSmith2021's makeIC_clump.py with
# clump_radius_pc=0, which cleanly degenerates to a uniform box (verified:
# 0 clump particles, no carve-out). Validated against the STARBENCH
# semi-empirical formula (hu_smith_analytic_check.py), the same way Hu et
# al. and Smith et al. (2021) validate their own D-type test -- see README.
n_threads=${n_threads:=8}          # Number of threads to use
gas_density=${gas_density:=100}    # Uniform gas density in atom/cm^3 (Hu et al. 2017)
gas_particle_mass=${gas_mass:=4.0} # Gas particle mass (Msun) -- Hu et al.'s own coarsest resolution
star_mass=${star_mass:=19.2}       # Mass of each source in M_sun (Q_H ~ 2.5e48/s)
n_stars=${n_stars:=4}               # Number of co-located sources at the box center
star_type=${star_type:="single_star"}
density_factor=${density_factor:=1.0} # Unused at clump_radius_pc=0, kept for makeIC_clump.py compatibility
clump_distance_pc=${clump_distance_pc:=20.0} # Unused at clump_radius_pc=0
clump_radius_pc=${clump_radius_pc:=0.0}      # 0 = no clump, uniform box (verified: 0 clump particles)
with_cooling=${with_cooling:=1}
L=${boxsize:=0.1}                  # boxsize in kpc (100 pc, 50 pc half-width -- see README for
                                    # the STARBENCH-curve-based box-size derivation)
time_end=${time_end:=8.182e-3}     # TimeIntegration:time_end override (internal units) -- 8.0 Myr,
                                    # matching Hu et al. (2017) Fig. A1's reported test duration.
dt_max=${dt_max:=1e-5}             # TimeIntegration:dt_max override (internal units). Must stay
                                    # below time_end (engine_config() errors otherwise) -- lower
                                    # this if you shorten time_end below the default.
n_cells=${n_cells:=3}               # must match Scheduler:max_top_level_cells in params.yml
nside=${nside:=0}                   # GEARFeedback:HII_angular_nside override -- 0 by default (Hu et al.'s
                                     # own algorithm has no angular splitting; this test is uniform anyway)
rebuild_time_myr=${rebuild_time_myr:=0.01} # GEARFeedback:HII_rebuild_time_Myr override (only used if adaptive_rebuild_cadence=0)
adaptive_rebuild_cadence=${adaptive_rebuild_cadence:=1} # GEARFeedback:HII_adaptive_rebuild_cadence override -- on by default (T11)
rebuild_safety_factor=${rebuild_safety_factor:=0.2} # GEARFeedback:HII_rebuild_safety_factor override
deterministic=${deterministic:=0}   # GEARFeedback:HII_deterministic_boundary_ionization override
initial_metallicity=${initial_metallicity:=0} # GEARChemistry:initial_metallicity override --
                        # Z=0, literally pure hydrogen matching Hu/Smith's own medium; gives this
                        # code's own T_i~47500 K floor, not the papers' flat 1e4 K. See README.
max_search_radius=${max_search_radius:=0.04} # Stars:HII_max_search_radius override (internal units, 40 pc,
                                    # just inside the 50 pc box half-width)
max_retry_full_buffer=${max_retry_full_buffer:=30} # Stars:HII_max_retry_full_buffer override
max_radius_expansion_tries=${max_radius_expansion_tries:=5} # Stars:HII_max_radius_expansion_tries override
radius_expansion_factor=${radius_expansion_factor:=1.1} # Stars:HII_radius_expansion_factor override
run_name=${run_name:=""}
restart=${restart:=0}

if [ -e ICs_husmith.hdf5 ]
then
    rm ICs_husmith.hdf5
fi

if [ ! -e ICs_husmith.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC_clump.py --boxsize $L --rho $gas_density \
                --mass $gas_particle_mass --star_mass $star_mass \
                --n_stars $n_stars \
                --density_factor $density_factor \
                --clump_distance_pc $clump_distance_pc \
                --clump_radius_pc $clump_radius_pc \
                --n_cells $n_cells --star_type $star_type \
            -o ICs_husmith.hdf5
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
               --sync --limiter $runtime_param --threads=$n_threads \
               -P TimeIntegration:time_end:$time_end \
               -P TimeIntegration:dt_max:$dt_max \
               -P GEARFeedback:HII_angular_nside:$nside \
               -P GEARFeedback:HII_rebuild_time_Myr:$rebuild_time_myr \
               -P GEARFeedback:HII_adaptive_rebuild_cadence:$adaptive_rebuild_cadence \
               -P GEARFeedback:HII_rebuild_safety_factor:$rebuild_safety_factor \
               -P GEARFeedback:HII_deterministic_boundary_ionization:$deterministic \
               -P GEARChemistry:initial_metallicity:$initial_metallicity \
               -P Stars:HII_max_search_radius:$max_search_radius \
               -P Stars:HII_max_retry_full_buffer:$max_retry_full_buffer \
               -P Stars:HII_max_radius_expansion_tries:$max_radius_expansion_tries \
               -P Stars:HII_radius_expansion_factor:$radius_expansion_factor \
               params.yml 2>&1 | tee output.log
else
../../../../swift --hydro --stars --external-gravity --feedback \
                --sync --limiter $runtime_param --threads=$n_threads \
               -P TimeIntegration:time_end:$time_end \
               -P TimeIntegration:dt_max:$dt_max \
               -P GEARFeedback:HII_angular_nside:$nside \
               -P GEARFeedback:HII_rebuild_time_Myr:$rebuild_time_myr \
               -P GEARFeedback:HII_adaptive_rebuild_cadence:$adaptive_rebuild_cadence \
               -P GEARFeedback:HII_rebuild_safety_factor:$rebuild_safety_factor \
               -P GEARFeedback:HII_deterministic_boundary_ionization:$deterministic \
               -P GEARChemistry:initial_metallicity:$initial_metallicity \
               -P Stars:HII_max_search_radius:$max_search_radius \
               -P Stars:HII_max_retry_full_buffer:$max_retry_full_buffer \
               -P Stars:HII_max_radius_expansion_tries:$max_radius_expansion_tries \
               -P Stars:HII_radius_expansion_factor:$radius_expansion_factor \
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
