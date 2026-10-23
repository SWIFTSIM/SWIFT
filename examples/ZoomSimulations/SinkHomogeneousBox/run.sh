#!/bin/bash

set -e
set -o pipefail

n_threads=${n_threads:=8}                     # Number of threads to use
level=${level:=5}                             # Number of gas particles = 2^(3*level)
embed=${embed:=8}                             # Periodic box side in gas box units (power of 2)
dm_ratio=${dm_ratio:=5.36}                    # Dark matter to gas density ratio in the gas box
delta=${delta:=0}                             # Density contrast of the gas box
with_star_formation=${with_star_formation:=0} # 1: GEAR star formation instead of sinks
cosmological=${cosmological:=0}               # 1: comoving ICs and --cosmology (not tested)

scripts_location="../../GEAR_ICs_and_SCRIPTS"

echo "Generating initial conditions for SinkHomogeneousBox..."
ic_params=""
cosmo_flag=""
if [ "$cosmological" -eq 1 ]; then
	ic_params="--cosmological"
	cosmo_flag="--cosmology"
fi

python3 makeIC.py --level $level --embed $embed --dm-ratio $dm_ratio \
	--delta $delta $ic_params -o ICs_zoom_homogeneous_box.hdf5 \
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
	$(cat zoom_runtime_params.txt) $parameter_file 2>&1 | tee output.log
