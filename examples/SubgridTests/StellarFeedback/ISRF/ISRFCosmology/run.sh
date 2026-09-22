#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

# pipefail: swift is piped into tee, whose exit code would hide a crash.
set -eo pipefail

config=${config:="free_field"}  #free_field, dust_absorption, photoelectric, photoelectric_dark, injection, injection_dusty or h2_shielded
redshift=${redshift:=0}         #Starting redshift, 0 runs without cosmology

# Physical inputs at the starting redshift, see README.
level_default=5
gas_density_default=1        # atom/cm^3
gas_mass_default=1           # Msun
temperature_default=100      # K
metallicity_default=0        # Z/Zsun, or an absolute mass fraction when scale_metallicity=0
scale_metallicity_default=1  # GEARChemistry:scale_initial_metallicity
u_pe_default=0              # erg/g
u_lw_default=0               # erg/g
nH2_ratio_default=1e-8
h2_self_shielding_default=3
propagation_default=1
c_hyp_pin_default=0          # km/s, 0 = off
disable_cooling_default=0    # GrackleCooling:disable_cooling_for_debugging
star_mass_default=0          # Msun, 0 = no star
duration_default=0.22283119056961848  # internal time (218 Myr, z = 9 to a = 0.125)
snapshots_default=40
steps_default=2230
max_star_dt_myr_default=1e-7  # Myr, young-star step cap
star_age_default=0
case "$config" in
    free_field)
	u_pe_default=6.9955e4
	u_lw_default=6.9955e4
	nH2_ratio_default=2e-4
	h2_self_shielding_default=0  # the check's closed form (A2) is unshielded
	;;
    dust_absorption)
	metallicity_default=1
	u_pe_default=1e5
	u_lw_default=1e5
	c_hyp_pin_default=4
	;;
    photoelectric|photoelectric_dark)
	metallicity_default=1
	temperature_default=10
	c_hyp_pin_default=1
	duration_default=3.07e-6
	snapshots_default=20
	steps_default=200
	if [ "$config" = "photoelectric" ]; then
	    u_pe_default=1.6e9
	    u_lw_default=1.6e9
	fi
	;;
    injection|injection_dusty)
	gas_density_default=1e3
	gas_mass_default=0.1
	temperature_default=50
	propagation_default=0
	star_mass_default=29.7
	duration_default=3e-5
	snapshots_default=3
	steps_default=3
	star_age_default=1e-6
	max_star_dt_myr_default=10  # the star steps at dt_max
	if [ "$config" = "injection_dusty" ]; then
	    # Grackle's own solar metal mass fraction, absolute, so the
	    # dust-to-gas ratio relative to the Milky Way is exactly 1 and both
	    # bands are strongly and differently attenuated (see the README).
	    metallicity_default=0.01295
	    scale_metallicity_default=0
	    # Metal cooling would otherwise drive the gas onto time-bins below
	    # the star's, and the check reads the star's Delta_t from the step
	    # table, which then reports the step's dt instead. The gas state is
	    # held fixed so the snapshot's h and rho are the ones the injection
	    # actually used.
	    disable_cooling_default=1
	fi
	;;
    h2_shielded)
	gas_density_default=1e3
	gas_mass_default=0.1
	temperature_default=50
	nH2_ratio_default=3.1e-5
	h2_self_shielding_default=3
	star_mass_default=29.7
	duration_default=2.6e-8
	star_age_default=1e-6  # a star of age 0 skips its first injection
	snapshots_default=41  # outputs never on a star step boundary (256 steps)
	steps_default=162
	;;
    *)
	echo "Unknown config '$config'."
	exit 1
	;;
esac

n_threads=${n_threads:=8}
level=${level:=$level_default} #N = (2**level)**3 gas particles
gas_density=${gas_density:=$gas_density_default}
gas_mass=${gas_mass:=$gas_mass_default}
temperature=${temperature:=$temperature_default}
metallicity=${metallicity:=$metallicity_default}
scale_metallicity=${scale_metallicity:=$scale_metallicity_default}
disable_cooling=${disable_cooling:=$disable_cooling_default}
u_pe=${u_pe:=$u_pe_default}
u_lw=${u_lw:=$u_lw_default}
nH2_ratio=${nH2_ratio:=$nH2_ratio_default}
h2_self_shielding=${h2_self_shielding:=$h2_self_shielding_default} #0, 2 or 3
propagation=${propagation:=$propagation_default}
c_hyp_pin=${c_hyp_pin:=$c_hyp_pin_default}
star_mass=${star_mass:=$star_mass_default}
star_age=${star_age:=$star_age_default} #Internal time units
duration=${duration:=$duration_default} #Proper time, internal units
snapshots=${snapshots:=$snapshots_default}
steps=${steps:=$steps_default} #Number of dt_max steps over the run
max_star_dt_myr=${max_star_dt_myr:=$max_star_dt_myr_default}
extinction_path=${extinction_path:="constant_kernel_path"}
extinction_path_in_kernel_radii=${extinction_path_in_kernel_radii:=1.0}
extinction_jeans_temperature_cap_K=${extinction_jeans_temperature_cap_K:=40}
run_name=${run_name:=""}
swift=${swift:="../../../../../swift"}

glass_n=$((2**level))
if [ ! -e "glassCube_${glass_n}.hdf5" ]; then
    ./getGlass.sh $glass_n
fi
if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
    "$scripts_location"/getGrackleCoolingTable.sh
fi
if [ ! -e POPIIsw.h5 ]; then
    "$scripts_location"/getChemistryTable.sh
fi

# Stop here on a table the radiation reader cannot use, rather than
# aborting at start-up once the initial conditions are built.
"$scripts_location"/checkRadiationTable.sh POPIIsw.h5 --with-isrf || exit 1

python3 makeIC.py --level $level --rho $gas_density --mass $gas_mass \
    --temperature $temperature --redshift $redshift --u-pe $u_pe \
    --u-lw $u_lw --star-mass $star_mass --star-age $star_age -o ICs_isrf_cosmology.hdf5

eval "$(python3 cosmo_timeline.py --redshift $redshift --duration $duration \
    --snapshots $snapshots --steps $steps)"

time_args=(-P TimeIntegration:dt_max:$dt_max
	   -P Snapshots:delta_time:$delta_time
	   -P Statistics:delta_time:$delta_time)
if [ "$(python3 -c "print(int(float('$redshift') > 0))")" = "1" ]; then
    time_args+=(--cosmology
		-P Cosmology:a_begin:$a_begin
		-P Cosmology:a_end:$a_end
		-P Snapshots:scale_factor_first:$a_begin
		-P Statistics:scale_factor_first:$a_begin)
else
    # The star is born at time 0 and the run starts at its age.
    time_begin=$star_age
    time_end=$(python3 -c "print(repr($star_age + $time_end))")
    time_args+=(-P TimeIntegration:time_begin:$time_begin
		-P TimeIntegration:time_end:$time_end
		-P Snapshots:time_first:$time_begin
		-P Statistics:time_first:$time_begin)
fi

rm -rf snap
mkdir snap

"$swift" --hydro --stars --external-gravity --feedback --cooling \
    --sync --limiter --verbose=0 --threads=$n_threads "${time_args[@]}" \
    -P GEARChemistry:initial_metallicity:$metallicity \
    -P GEARChemistry:scale_initial_metallicity:$scale_metallicity \
    -P GrackleCooling:disable_cooling_for_debugging:$disable_cooling \
    -P GrackleCooling:initial_nH2I_to_nH_ratio:$nH2_ratio \
    -P GrackleCooling:H2_self_shielding:$h2_self_shielding \
    -P GEARFeedback:ISRF_propagation:$propagation \
    -P GEARFeedback:ISRF_extinction_path:$extinction_path \
    -P GEARFeedback:ISRF_extinction_path_in_kernel_radii:$extinction_path_in_kernel_radii \
    -P GEARFeedback:ISRF_extinction_jeans_temperature_cap_K:$extinction_jeans_temperature_cap_K \
    -P GEARFeedback:ISRF_c_hyp_pin_for_debugging:$c_hyp_pin \
    -P Stars:max_timestep_young_Myr:$max_star_dt_myr \
    params.yml 2>&1 | tee output.log

if [ -n "$run_name" ]; then
    mkdir -p "$run_name"
    mv snap output.log timesteps.txt statistics.txt used_parameters.yml \
       unused_parameters.yml "$run_name"
fi
