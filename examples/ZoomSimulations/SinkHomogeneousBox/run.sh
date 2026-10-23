#!/bin/bash

set -e
set -o pipefail

# run_case: local (level 3, laptop test), level5 (production), level8 (Jed)
run_case=${run_case:=local}
case $run_case in
	local)
		level=${level:=3}
		n_threads=${n_threads:=2}
		case_params="-P TimeIntegration:time_end:0.1"
		;;
	level5)
		level=${level:=5}
		n_threads=${n_threads:=8}
		case_params=""
		;;
	level8)
		level=${level:=8}
		n_threads=${n_threads:=72}
		case_params="-P Restarts:delta_hours:2 -P Snapshots:delta_time:0.05"
		;;
	*)
		echo "Unknown run_case $run_case (local, level5 or level8)"
		exit 1
		;;
esac

embed=${embed:=8}                             # Periodic box side in gas box units (power of 2)
dm_ratio=${dm_ratio:=5.36}                    # Dark matter to gas density ratio in the gas box
delta=${delta:=0}                             # Density contrast of the gas box
jeans_length=${jeans_length:=0.250}           # Jeans wavelength in gas box units (sets the temperature)
with_star_formation=${with_star_formation:=0} # 1: GEAR star formation instead of sinks
cosmological=${cosmological:=0}               # 1: comoving ICs and --cosmology (not tested)

scripts_location="../../GEAR_ICs_and_SCRIPTS"

ic_params=""
cosmo_flag=""
if [ "$cosmological" -eq 1 ]; then
	ic_params="--cosmological"
	cosmo_flag="--cosmology"
fi

echo "Generating initial conditions for SinkHomogeneousBox (level $level)..."
python3 makeIC.py --level $level --embed $embed --dm-ratio $dm_ratio \
	--delta $delta --lJ $jeans_length $ic_params -o ICs_zoom_homogeneous_box.hdf5 \
	--zoom-params-file zoom_runtime_params.txt

if [ ! -e CloudyData_UVB=HM2012_high_density.h5 ]; then
	echo "Fetching the Cloudy tables required by Grackle..."
	$scripts_location/getGrackleCoolingTable.sh
fi

if [ ! -e POPII.hdf5 ]; then
	echo "Fetching the chemistry tables..."
	$scripts_location/getChemistryTable.sh --with-winds
fi

if [ "$with_star_formation" -eq 0 ]; then
	runtime_params="--sinks"
	parameter_file="params_sinks.yml"
else
	runtime_params="--star-formation"
	parameter_file="params_star_formation.yml"
fi

rm -rf snap statistics.txt timesteps.txt output.log
mkdir snap

../../../swift --hydro $runtime_params --stars --self-gravity --feedback \
	--cooling --sync --limiter --zoom $cosmo_flag --threads=$n_threads \
	$(cat zoom_runtime_params.txt) $case_params $parameter_file 2>&1 | tee output.log
