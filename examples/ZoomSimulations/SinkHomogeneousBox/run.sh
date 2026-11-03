#!/bin/bash

set -e
set -o pipefail

# Every variable can be set on the command line, for example
#   level=3 padding=2 n_threads=2 time_end=0.1 ./run.sh
level=${level:=5}                             # levelmax: 2^(3 level) gas particles
levelmin=${levelmin:=}                        # Coarsest background level (empty: level - 3)
padding=${padding:=8}                         # Cells of each coarser level around the finer one
box_size=${box_size:=}                        # Box side: kpc, or comoving Mpc if cosmological=1 (empty: 8 gas boxes)
n_threads=${n_threads:=8}
gas_mass=${gas_mass:=50}                      # Gas particle mass in Msun
jeans_length=${jeans_length:=0.250}           # Jeans wavelength in gas box units (sets the temperature)
dm_ratio=${dm_ratio:=5.36}                    # Dark matter to gas density ratio in the gas box
delta=${delta:=0}                             # Density contrast of the gas box
with_star_formation=${with_star_formation:=0} # 1: GEAR star formation instead of sinks
cosmological=${cosmological:=0}               # 1: comoving ICs and --cosmology
time_end=${time_end:=0.5}                     # End time in internal units (cosmological=0)
a_end_factor=${a_end_factor:=2}               # a_end / a_begin (cosmological=1)
snapshot_delta=${snapshot_delta:=0.01}        # Time between snapshots (cosmological=0)
snapshot_delta_cosmo=${snapshot_delta_cosmo:=1.02} # Ratio of scale factors between snapshots, > 1 (cosmological=1)
velocity_fraction=${velocity_fraction:=0.1}   # Random velocities in units of the gas sound speed (cosmological=1)
max_bkg_cdim=${max_bkg_cdim:=91}              # Cap on ZoomRegion:bkg_top_level_cells
swift=${swift:=../../../swift}
extra_params=${extra_params:=}                # More -P options for swift

scripts_location="../../GEAR_ICs_and_SCRIPTS"

# Lengths tied to the particle spacing scale as (gas_mass / 50)^(1/3).
# The values in the parameter files are the ones for 50 Msun.
length_scale=$(python3 -c "print(($gas_mass / 50.0) ** (1.0 / 3.0))")
h_max=$(python3 -c "print(0.08 * $length_scale)")
softening=$(python3 -c "print(0.005 * $length_scale)")
cut_off_radius=$(python3 -c "print(0.005 * $length_scale)")

ic_params="--level $level --padding $padding --mass $gas_mass --h-max $h_max"
ic_params="$ic_params --softening $softening --dm-ratio $dm_ratio --delta $delta"
ic_params="$ic_params --lJ $jeans_length --max-bkg-cdim $max_bkg_cdim"
if [ -n "$levelmin" ]; then
	ic_params="$ic_params --levelmin $levelmin"
fi
if [ -n "$box_size" ]; then
	ic_params="$ic_params --box-size $box_size"
fi

cosmo_flag=""
run_params="-P Gravity:max_physical_baryon_softening:$softening"
run_params="$run_params -P Gravity:max_physical_DM_softening:$softening"
if [ "$cosmological" -eq 1 ]; then
	ic_params="$ic_params --cosmological --a-end-factor $a_end_factor"
	ic_params="$ic_params --velocity-fraction $velocity_fraction"
	cosmo_flag="--cosmology"
	run_params="$run_params -P Snapshots:delta_time:$snapshot_delta_cosmo"
else
	# makeIC.py writes the comoving SPH:h_max in cosmological runs
	run_params="$run_params -P SPH:h_max:$h_max"
	run_params="$run_params -P TimeIntegration:time_end:$time_end"
	run_params="$run_params -P Snapshots:delta_time:$snapshot_delta"
fi

echo "Generating initial conditions for SinkHomogeneousBox (level $level)..."
python3 makeIC.py $ic_params -o ICs_zoom_homogeneous_box.hdf5 \
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
	run_params="$run_params -P GEARSink:cut_off_radius:$cut_off_radius"
else
	runtime_params="--star-formation"
fi

rm -rf snap statistics.txt timesteps.txt output.log
mkdir snap

$swift --hydro $runtime_params --stars --self-gravity --feedback \
	--cooling --sync --limiter --zoom $cosmo_flag --threads=$n_threads \
	$(cat zoom_runtime_params.txt) $run_params $extra_params params.yml 2>&1 | tee output.log
