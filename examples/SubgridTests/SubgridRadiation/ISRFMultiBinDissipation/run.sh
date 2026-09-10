#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../GEAR_ICs_and_SCRIPTS"

# make run.sh fail if a subcommand fails. pipefail matters here specifically:
# swift's own exit code is piped into `tee output.log`, and without it a
# crashed/errored swift run still lets the pipeline "succeed" (tee's own
# exit code), silently producing a run_name output directory that looks
# complete but only contains a startup-error log.
set -eo pipefail

n_threads=${n_threads:=8}
gas_density=${gas_density:=1e3}                # cold-phase n_H, atom/cm^3
gas_particle_mass=${gas_mass:=0.1}              # cold-phase particle mass, Msun
level=${level:=6}
variant=${variant:="twophase"}                  # uniform | twophase
bin_delta=${bin_delta:=2}
bulk_temperature_K=${bulk_temperature_K:=500}   # cold phase
interface_offset_h=${interface_offset_h:=1.0}
v_over_c_hyp=${v_over_c_hyp:=0.25}
star_velocity_km_s=${star_velocity_km_s:=-1}    # <0: derive from v_over_c_hyp
star_mass=${star_mass:=29.7}
star_type=${star_type:="single_star"}
n_stars=${n_stars:=3}
c_hyp_margin=${c_hyp_margin:=0.5}
c_hyp_pin=${c_hyp_pin:=0}                       # LW_FUV_c_hyp_pin_for_debugging, km/s
alpha_max=${alpha_max:=0.25}
alpha_pin=${alpha_pin:=0}
propagation=${propagation:=1}
initial_metallicity=${initial_metallicity:=1}
time_end=${time_end:=5e-4}
dt_max=${dt_max:=1e-4}
dt_min=${dt_min:=1e-14}
delta_time=${delta_time:=5e-5}
run_name=${run_name:=""}

# Remove the ICs
if [ -e ICs_isrf_multibin_dissipation.hdf5 ]
then
    rm ICs_isrf_multibin_dissipation.hdf5
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
	--mass $gas_particle_mass --variant $variant --bin-delta $bin_delta \
	--bulk-temperature-K $bulk_temperature_K \
	--interface-offset-h $interface_offset_h --v-over-c-hyp $v_over_c_hyp \
	--star-velocity-km-s $star_velocity_km_s --c-hyp-margin $c_hyp_margin \
	--cfl 0.1 --star_mass $star_mass --star_type $star_type \
	--n-stars $n_stars -o ICs_isrf_multibin_dissipation.hdf5)
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

printf "Running simulation..."

# --external-gravity with no Potential: block gives zero gravitational
# force. --cooling is on so cooling_init() actually resolves
# chemistry_data (including local_dust_to_gas_ratio); the per-particle
# energy update itself is skipped via
# GrackleCooling:disable_cooling_for_debugging in params.yml, which keeps
# the two-phase temperature (hence timestep-bin) structure frozen.
../../../../swift --hydro --stars --external-gravity --feedback --cooling \
		   --sync --limiter --verbose=0 --threads=$n_threads \
		   -P TimeIntegration:time_end:$time_end \
		   -P TimeIntegration:dt_max:$dt_max \
		   -P TimeIntegration:dt_min:$dt_min \
		   -P Snapshots:delta_time:$delta_time \
		   -P GEARChemistry:initial_metallicity:$initial_metallicity \
		   -P GEARFeedback:LW_FUV_propagation:$propagation \
		   -P GEARFeedback:LW_FUV_c_hyp_margin:$c_hyp_margin \
		   -P GEARFeedback:LW_FUV_c_hyp_pin_for_debugging:$c_hyp_pin \
		   -P GEARFeedback:LW_FUV_dissipation_alpha_max:$alpha_max \
		   -P GEARFeedback:LW_FUV_dissipation_alpha_pin_for_debugging:$alpha_pin \
		   -P SPH:initial_temperature:0 \
		   params.yml 2>&1 | tee output.log

# Per-run metrics (report-only, never exits nonzero); see README.
python3 isrf_multibin_dissipation_check.py --c-hyp-margin $c_hyp_margin --c-hyp-pin $c_hyp_pin

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
	mv multibin_metrics.json $run_name
	mv multibin_ic.json $run_name
    fi
fi
