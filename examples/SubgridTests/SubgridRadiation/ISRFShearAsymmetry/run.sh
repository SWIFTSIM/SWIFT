#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../GEAR_ICs_and_SCRIPTS"

# make run.sh fail if a subcommand fails. pipefail matters here specifically:
# swift's own exit code is piped into `tee output.log`, and without it a
# crashed/errored swift run still lets the pipeline "succeed" (tee's own
# exit code), silently producing a run_name output directory that looks
# complete but only contains a startup-error log.
set -eo pipefail

n_threads=${n_threads:=8}  #Number of threads to use
gas_density=${gas_density:=1e3} #Gas density in atom/cm^3
gas_particle_mass=${gas_mass:=0.1} #Mass of the gas particles (Msun)
level=${level:=5} #Resolution level: N = (2**level)**3 gas particles
variant=${variant:="shear"}
source_geometry=${source_geometry:="blobs"}
v_shear=${v_shear:=1.0}                # km/s, full v_rel between the streams
layer_width_h=${layer_width_h:=4.0}
density_ratio=${density_ratio:=2.0}
bulk_temperature_K=${bulk_temperature_K:=1000}
pulse_amplitude=${pulse_amplitude:=1.0}
pulse_sigma_h=${pulse_sigma_h:=2.0}
c_hyp_margin=${c_hyp_margin:=0.5}
c_hyp_pin=${c_hyp_pin:=0}
alpha_max=${alpha_max:=0.25}
alpha_pin=${alpha_pin:=0}
propagation=${propagation:=1}
initial_metallicity=${initial_metallicity:=0.05}
time_end=${time_end:=1.0e-3}
dt_max=${dt_max:=1e-5}
dt_min=${dt_min:=1e-14}
delta_time=${delta_time:=1.0e-4}
expect_stable=${expect_stable:=1}
run_name=${run_name:=""}

# Remove the ICs
if [ -e ICs_isrf_shear_asymmetry.hdf5 ]
then
    rm ICs_isrf_shear_asymmetry.hdf5
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
	--mass $gas_particle_mass --variant $variant \
	--source-geometry $source_geometry --v-shear-km-s $v_shear \
	--layer-width-h $layer_width_h --density-ratio $density_ratio \
	--bulk-temperature-K $bulk_temperature_K \
	--pulse-amplitude $pulse_amplitude --pulse-sigma-h $pulse_sigma_h \
    -o ICs_isrf_shear_asymmetry.hdf5)
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
# GrackleCooling:disable_cooling_for_debugging in params.yml, which
# isolates the propagation field from gas thermal/dynamical response
# without leaving chemistry_data unresolved (see README).
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

# Per-run metrics (asymmetry, KH-contamination control); see README.
stable_flag="--expect-stable"
if [ "$expect_stable" = "0" ]; then
    stable_flag="--expect-unstable"
fi
python3 isrf_shear_asymmetry_check.py --variant $variant --source-geometry $source_geometry \
    --v-shear $v_shear --c-hyp-pin $c_hyp_pin --c-hyp-margin $c_hyp_margin

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
	mv shear_metrics.json $run_name
    fi
fi
