#!/bin/bash


# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../GEAR_ICs_and_SCRIPTS"
# make run.sh fail if a subcommand fails. pipefail matters here specifically:
# swift's own exit code is piped into `tee output.log`, and without it a
# crashed/errored swift run still lets the pipeline "succeed" (tee's own
# exit code), silently producing a run_name output directory that looks
# complete but only contains a startup-error log.
set -eo pipefail

# Defaults reproduce STARBENCH's late-phase D-type expansion test; see README.
n_threads=${n_threads:=8}  #Number of threads to use
gas_density=${gas_density:=3115} #Gas density in atom/cm^3
gas_particle_mass=${gas_mass:=0.1} #Mass of the gas particles (Msun); see README's Resolution section
star_mass=${star_mass:=26.75} #Star mass (Msun); see README's Star mass section
star_type=${star_type:="single_star"}
with_cooling=${with_cooling:=1}
L=${boxsize:=0.01} #boxsize in kpc
time_end=${time_end:=3.068e-3} #TimeIntegration:time_end override (internal units) -- 3.0 Myr
dt_max=${dt_max:=1e-5} #TimeIntegration:dt_max override (internal units). Must stay
                        #below time_end (engine_config() errors otherwise) -- lower
                        #this if you shorten time_end below the default.
initial_metallicity=${initial_metallicity:=0} #GEARChemistry:initial_metallicity override; see README's Temperature section
rebuild_time_myr=${rebuild_time_myr:=0.01} #GEARFeedback:HII_rebuild_time_Myr override; see README's Cadence section
max_search_radius=${max_search_radius:=0.0049} #Stars:HII_max_search_radius override (internal units, 4.9 pc)
max_retry_full_buffer=${max_retry_full_buffer:=30} #Stars:HII_max_retry_full_buffer override
max_radius_expansion_tries=${max_radius_expansion_tries:=5} #Stars:HII_max_radius_expansion_tries override
radius_expansion_factor=${radius_expansion_factor:=1.1} #Stars:HII_radius_expansion_factor override
run_name=${run_name:=""}
restart=${restart:=0}

# Remove the ICs
if [ -e ICs_starbench.hdf5 ]
then
    rm ICs_starbench.hdf5
fi

#Create the ICs if they do not exist
if [ ! -e ICs_starbench.hdf5 ]
then
    echo "Generating initial conditions to run the example..."
    python3 makeIC.py --boxsize $L --rho $gas_density \
		--mass $gas_particle_mass --star_mass $star_mass \
		--star_type $star_type \
	    -o ICs_starbench.hdf5
fi

# Get the Grackle cooling table
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

# Create output directory
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
		   -P GEARChemistry:initial_metallicity:$initial_metallicity \
		   -P GEARFeedback:HII_rebuild_time_Myr:$rebuild_time_myr \
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
	       -P GEARChemistry:initial_metallicity:$initial_metallicity \
	       -P GEARFeedback:HII_rebuild_time_Myr:$rebuild_time_myr \
		   -P Stars:HII_max_search_radius:$max_search_radius \
		   -P Stars:HII_max_retry_full_buffer:$max_retry_full_buffer \
		   -P Stars:HII_max_radius_expansion_tries:$max_radius_expansion_tries \
		   -P Stars:HII_radius_expansion_factor:$radius_expansion_factor \
		   params.yml 2>&1 | tee output.log
fi

# Compare to the analytical solutions
python3 starbench_analytic_check.py -s "snap/snapshot_*.hdf5"

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
    fi
fi
