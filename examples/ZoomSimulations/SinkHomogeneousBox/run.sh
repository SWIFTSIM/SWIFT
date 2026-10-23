#!/bin/bash

set -e
set -o pipefail

# run_case: level5 (default, production), local (level 3, laptop test), level8 (Jed)
run_case=${run_case:=level5}
case $run_case in
	level5)
		level=${level:=5}
		n_threads=${n_threads:=8}
		case_params=""
		;;
	local)
		level=${level:=3}
		n_threads=${n_threads:=2}
		case_params="-P TimeIntegration:time_end:0.1"
		;;
	level8)
		level=${level:=8}
		n_threads=${n_threads:=72}
		case_params="-P Restarts:delta_hours:2 -P Snapshots:delta_time:0.05"
		;;
	*)
		echo "Unknown run_case $run_case (level5, local or level8)"
		exit 1
		;;
esac

gas_mass=${gas_mass:=50}                      # Gas particle mass in Msun
embed=${embed:=8}                             # Periodic box side in gas box units (power of 2)
dm_ratio=${dm_ratio:=5.36}                    # Dark matter to gas density ratio in the gas box
delta=${delta:=0}                             # Density contrast of the gas box
jeans_length=${jeans_length:=0.250}           # Jeans wavelength in gas box units (sets the temperature)
with_star_formation=${with_star_formation:=0} # 1: GEAR star formation instead of sinks
cosmological=${cosmological:=0}               # 1: comoving ICs and --cosmology (not tested)

scripts_location="../../GEAR_ICs_and_SCRIPTS"

# Lengths tied to the particle spacing scale as (gas_mass / 50)^(1/3).
# The values in the parameter files are the ones for 50 Msun.
length_scale=$(python3 -c "print(($gas_mass / 50.0) ** (1.0 / 3.0))")
h_max=$(python3 -c "print(0.08 * $length_scale)")
softening=$(python3 -c "print(0.005 * $length_scale)")
cut_off_radius=$(python3 -c "print(0.005 * $length_scale)")

ic_params="--mass $gas_mass --h-max $h_max --softening $softening"
cosmo_flag=""
mass_params="-P Gravity:max_physical_baryon_softening:$softening"
mass_params="$mass_params -P Gravity:max_physical_DM_softening:$softening"
if [ "$cosmological" -eq 1 ]; then
	ic_params="$ic_params --cosmological"
	cosmo_flag="--cosmology"
else
	# makeIC.py writes the comoving SPH:h_max in cosmological runs
	mass_params="$mass_params -P SPH:h_max:$h_max"
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
	mass_params="$mass_params -P GEARSink:cut_off_radius:$cut_off_radius"
else
	runtime_params="--star-formation"
	parameter_file="params_star_formation.yml"
fi

rm -rf snap statistics.txt timesteps.txt output.log
mkdir snap

../../../swift --hydro $runtime_params --stars --self-gravity --feedback \
	--cooling --sync --limiter --zoom $cosmo_flag --threads=$n_threads \
	$(cat zoom_runtime_params.txt) $mass_params $case_params $parameter_file 2>&1 | tee output.log
