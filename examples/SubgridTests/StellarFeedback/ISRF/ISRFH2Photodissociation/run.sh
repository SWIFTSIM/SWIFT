#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

# make run.sh fail if a subcommand fails. pipefail matters here specifically:
# swift's own exit code is piped into `tee output.log`, and without it a
# crashed/errored swift run still lets the pipeline "succeed" (tee's own
# exit code), silently producing a run_name output directory that looks
# complete but only contains a startup-error log.
set -eo pipefail

config=${config:="thin"}  #H2 column regime: thin or thick

# The two regimes differ only in the initial H2 abundance, which sets the
# column Grackle's local self-shielding estimate sees; the box is identical.
case "$config" in
    thin)
	nH2_ratio_default=5e-10
	;;
    thick)
	nH2_ratio_default=3.1e-5
	;;
    *)
	echo "Unknown config '$config'. Use config=thin or config=thick."
	exit 1
	;;
esac

n_threads=${n_threads:=8}  #Number of threads to use
gas_density=${gas_density:=1e3} #Gas density in atom/cm^3
gas_particle_mass=${gas_mass:=0.1} #Mass of the gas particles (Msun)
star_mass=${star_mass:=29.7} #Star mass (Msun)
star_type=${star_type:="single_star"}
level=${level:=5} #Resolution level: N = (2**level)**3 gas particles
nH2_ratio=${nH2_ratio:=$nH2_ratio_default} #GrackleCooling:initial_nH2I_to_nH_ratio override
h2_self_shielding=${h2_self_shielding:=3} #GrackleCooling:H2_self_shielding override (0=off, 3=local Jeans length)
time_end=${time_end:=2.6e-8} #TimeIntegration:time_end override (internal units)
dt_max=${dt_max:=1.6e-10} #TimeIntegration:dt_max override (internal units)
delta_time=${delta_time:=6.5e-10} #Snapshots:delta_time override (internal units)
max_star_dt_myr=${max_star_dt_myr:=1e-7} #Stars:max_timestep_young_Myr override
min_star_dt_myr=${min_star_dt_myr:=1e-9} #Stars:min_star_timestep_Myr override
initial_metallicity=${initial_metallicity:=0} #GEARChemistry:initial_metallicity override (Z/Zsun)
run_name=${run_name:=""}

# Remove the ICs
if [ -e ICs_isrf_h2_photodissociation.hdf5 ]
then
    rm ICs_isrf_h2_photodissociation.hdf5
fi

glass_n=$((2**level))
glass_file="glassCube_${glass_n}.hdf5"
if [ ! -e "$glass_file" ]
then
    echo "Fetching initial glass file (${glass_file})..."
    ./getGlass.sh $glass_n
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

echo "Generating initial conditions to run the example..."
ic_output=$(python3 makeIC.py --level $level --rho $gas_density \
	--mass $gas_particle_mass --star_mass $star_mass \
	--star_type $star_type \
    -o ICs_isrf_h2_photodissociation.hdf5)
echo "$ic_output"

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

printf "Running the %s configuration...\n" "$config"

../../../../../swift --hydro --stars --external-gravity --feedback --cooling \
		   --sync --limiter --verbose=0 --threads=$n_threads \
		   -P TimeIntegration:time_end:$time_end \
		   -P TimeIntegration:dt_max:$dt_max \
		   -P Snapshots:delta_time:$delta_time \
		   -P Statistics:delta_time:$delta_time \
		   -P GEARChemistry:initial_metallicity:$initial_metallicity \
		   -P GrackleCooling:initial_nH2I_to_nH_ratio:$nH2_ratio \
		   -P GrackleCooling:H2_self_shielding:$h2_self_shielding \
		   -P Stars:max_timestep_young_Myr:$max_star_dt_myr \
		   -P Stars:min_star_timestep_Myr:$min_star_dt_myr \
		   params.yml 2>&1 | tee output.log

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
