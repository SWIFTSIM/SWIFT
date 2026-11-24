#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

# swift's exit status is piped into tee: without pipefail a refused or crashed run looks complete.
set -eo pipefail

leg=${leg:="main"}   # main, switch_off, stale_key, restart, cap_fine, cap_coarse, cap_half
n_threads=${n_threads:=8}   # Number of threads
ranks=${ranks:=1}   # MPI ranks, 2 needs an MPI build named swift_mpi
level=${level:=5}   # Resolution level: N = (2**level)**3 gas particles
gas_density=${gas_density:=1e3}   # Gas density in atom/cm^3
gas_mass=${gas_mass:=0.1}   # Gas particle mass in Msun
h2_ratio=${h2_ratio:=3e-4}   # Initial nH2 / nH
run_name=${run_name:=$leg}
ic_file=ICs_isrf_h2_star_lit.hdf5
params=params.yml
extra_ic_args=""
extra_swift_args=""
stale_key_arg=""
mode_flag="--cooling"   # the limited-speed legs run --temperature: the species stay at their start values

case "$leg" in
    main) switch=1 ;;
    switch_off) switch=0 ;;   # the line fields are not written: the check must refuse the run
    stale_key) switch=1 ; stale_key_arg="-P GeajayCooling:h2_transported_lw_rate_shielding:shielded_sum" ; extra_swift_args="-n 1" ;;   # the removed key must be ignored
    restart) switch=1 ; params=params_restart.yml ; extra_swift_args="-n 16" ;;   # limited light speed, stops after step 17
    cap_fine) switch=1 ; params=params_limited.yml ; mode_flag="--temperature" ; extra_swift_args="-n 71" ;;   # cap light speed, the star steps every gas step
    cap_coarse) switch=1 ; params=params_limited.yml ; mode_flag="--temperature" ; extra_swift_args="-n 71 -P Stars:max_timestep_young_Myr:1.84e-6" ;;   # the star steps every 8 gas steps
    cap_half) switch=1 ; params=params_limited.yml ; mode_flag="--temperature" ; extra_swift_args="-n 109 -P GEARFeedback:ISRF_packet_light_speed_fraction_of_c:0.5" ;;   # half the light speed, the same steady state
    *) echo "Unknown leg '$leg'." ; exit 1 ;;
esac

glass_file="glassCube_$((2**level)).hdf5"
if [ ! -e "$glass_file" ]; then
    echo "Fetching the glass file ($glass_file)..."
    ./getGlass.sh $((2**level))
fi
if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
    "$scripts_location"/getGrackleCoolingTable.sh
fi
"$scripts_location"/getRadiationTable.sh PopII_parsec_spectral.hdf5 || exit 1
"$scripts_location"/checkRadiationTable.sh PopII_parsec_spectral.hdf5 --with-isrf || exit 1

rm -f "$ic_file"
python3 makeIC.py --level $level --rho $gas_density --mass $gas_mass --h2_ratio $h2_ratio $extra_ic_args -o "$ic_file"

rm -rf snap restart
mkdir snap
swift_exe=../../../../../swift
[ "$ranks" -gt 1 ] && launcher="mpirun -np $ranks" && swift_exe=../../../../../swift_mpi || launcher=""

# Parameter overrides: the switch of the leg, and for stale_key the removed library key.
set +e
$launcher $swift_exe --hydro --stars --external-gravity --feedback $mode_flag \
    --sync --limiter --verbose=0 --threads=$n_threads \
    -P GEARFeedback:ISRF_H2_line_absorption:$switch \
    $stale_key_arg \
    $extra_swift_args $params 2>&1 | tee output.log
swift_status=${PIPESTATUS[0]}
set -e

if [ "$leg" = "stale_key" ]; then
    python3 ../ISRFH2StarLit/checkH2Attenuation.py stale_key --log output.log --exit-status $swift_status
    exit $?
fi
[ $swift_status -eq 0 ] || { echo "swift stopped with status $swift_status"; exit $swift_status; }

if [ "$leg" = "restart" ]; then
    echo "Stopped at the restart point. Resume with the same binary and the same options, with -r in place of -n 16."
    exit 0
fi

mkdir -p "$run_name"
cp used_parameters.yml output.log timesteps.txt fixture.json "$run_name"/
mv snap "$run_name"/
case "$leg" in
    cap_fine|cap_half)
        python3 ../ISRFH2StarLit/checkH2LimitedSpeed.py limited --run "$run_name"
        exit $? ;;
    cap_coarse)
        python3 ../ISRFH2StarLit/checkH2LimitedSpeed.py limited --run "$run_name" || exit $?
        if [ ! -d cap_fine ]; then
            echo "The ripple check needs the cap_fine leg and found no cap_fine directory."
            echo "Run leg=cap_fine ./run.sh, then: python3 checkH2LimitedSpeed.py ripple --coarse $run_name --fine cap_fine"
            exit 4
        fi
        python3 ../ISRFH2StarLit/checkH2LimitedSpeed.py ripple --coarse "$run_name" --fine cap_fine
        exit $? ;;
esac
python3 ../ISRFH2StarLit/checkH2Attenuation.py starlit \
    --first "$run_name"/snap/snapshot_0001.hdf5 \
    --last "$(ls "$run_name"/snap/snapshot_*.hdf5 | tail -n 1)" \
    --fixture "$run_name"/fixture.json --params "$run_name"/used_parameters.yml --startup-log "$run_name"/output.log \
    --timesteps "$run_name"/timesteps.txt
