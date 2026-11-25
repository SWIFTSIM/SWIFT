#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

# swift's exit status is piped into tee: without pipefail a refused or crashed run looks complete.
set -eo pipefail

leg=${leg:="main"}   # main, merge, switch_off, stale_key
n_threads=${n_threads:=8}   # Number of threads
ranks=${ranks:=1}   # MPI ranks, 2 needs an MPI build named swift_mpi
level=${level:=5}   # Resolution level: N = (2**level)**3 gas particles
gas_density=${gas_density:=1e3}   # Gas density in atom/cm^3
gas_mass=${gas_mass:=0.1}   # Gas particle mass in Msun
if [ "$leg" = merge ]; then
    default_h2_ratio=3e-6; default_n_steps=7
else
    default_h2_ratio=3e-4; default_n_steps=6
fi
h2_ratio=${h2_ratio:=$default_h2_ratio}   # Initial nH2 / nH
n_steps=${n_steps:=$default_n_steps}   # Steps of the run, keep below 8; swift -n N runs N+1 steps
in_line_pc=${in_line_pc:=0.2}   # merge leg: distance d in pc from star A to the receiver, star B sits 2 d behind star A
star_mass=${star_mass:=10}   # merge leg: mass of star A in Msun
star_mass_b=${star_mass_b:=50}   # merge leg: mass of star B in Msun
dt_max=${dt_max:=2.5e-14}   # merge leg: step in code time units, small so that the H2 stays frozen while the light arrives
run_name=${run_name:=$leg}
ic_file=ICs_isrf_h2_two_sources.hdf5
params=params.yml
extra_ic_args="${contrast:+--contrast $contrast}"
extra_swift_args="-n $((n_steps - 1))"
stale_key_arg=""

case "$leg" in
    main) switch=1 ;;
    merge)
        switch=1
        extra_ic_args="--contrast 1 --in_line_pc $in_line_pc --star_mass $star_mass --star_mass_b $star_mass_b"
        first_snapshot_time=$(awk -v d="$dt_max" 'BEGIN { printf "%.9g", d * 1.00001 }')   # just after a step end, so a snapshot holds whole steps
        extra_swift_args="$extra_swift_args -P TimeIntegration:dt_max:$dt_max -P Statistics:delta_time:$dt_max -P Snapshots:delta_time:$dt_max -P Snapshots:time_first:$first_snapshot_time"
        ;;
    switch_off) switch=0 ;;   # the line fields are not written: the check must refuse the run
    stale_key) switch=1 ; stale_key_arg="-P GeajayCooling:h2_transported_lw_rate_shielding:shielded_sum" ; extra_swift_args="-n 1" ;;   # the removed key must be ignored
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

swift_exe=../../../../../swift
[ "$ranks" -gt 1 ] && launcher="mpirun -np $ranks" && swift_exe=../../../../../swift_mpi || launcher=""

# Run swift on the star set $1 ("a", "b" or "ab") and keep the run in the directory $2.
# Parameter overrides: the switch of the leg, and for stale_key the removed library key.
run_stars() {
    rm -f "$ic_file"
    python3 makeIC.py --level $level --rho $gas_density --mass $gas_mass --h2_ratio $h2_ratio --stars $1 $extra_ic_args -o "$ic_file"
    rm -rf snap restart
    mkdir snap
    set +e
    $launcher $swift_exe --hydro --stars --external-gravity --feedback --cooling \
        --sync --limiter --verbose=0 --threads=$n_threads \
        -P GEARFeedback:ISRF_H2_line_absorption:$switch \
        -P GeajayCooling:initial_h2i_to_h_ratio:$h2_ratio \
        $stale_key_arg \
        $extra_swift_args $params 2>&1 | tee output.log
    swift_status=${PIPESTATUS[0]}
    set -e
    [ "$leg" = "stale_key" ] && return
    [ $swift_status -eq 0 ] || { echo "swift stopped with status $swift_status"; exit $swift_status; }
    rm -rf "$2"
    mkdir -p "$2"
    cp used_parameters.yml output.log timesteps.txt fixture.json "$2"/
    mv snap "$2"/
}

if [ "$leg" = "stale_key" ]; then
    run_stars ab
    python3 ../ISRFH2StarLit/checkH2Attenuation.py stale_key --log output.log --exit-status $swift_status
    exit $?
fi

# The check compares the two-star run with the sum of one-star runs of the same binary,
# and star A runs twice to measure the reproducibility of a run.
if [ "$leg" = "main" ] || [ "$leg" = "merge" ]; then
    run_stars a only_a
    run_stars b only_b
fi
[ "$leg" = "main" ] && run_stars a only_a2
run_stars ab "$run_name"

slots=$(sed -n 's/^#define RADIATION_PACKET_SLOTS //p' "$(dirname "$swift_exe")/config.h")
[ -n "$slots" ] || { echo "RADIATION_PACKET_SLOTS not found in $(dirname "$swift_exe")/config.h"; exit 1; }
if [ "$leg" = "merge" ]; then
    python3 ../ISRFH2StarLit/checkH2Attenuation.py merge \
        --first "$run_name"/snap/snapshot_0001.hdf5 \
        --last "$(ls "$run_name"/snap/snapshot_*.hdf5 | tail -n 1)" \
        --fixture "$run_name"/fixture.json --params "$run_name"/used_parameters.yml --startup-log "$run_name"/output.log \
        --timesteps "$run_name"/timesteps.txt \
        --only-a only_a --only-b only_b --slots "$slots"
    exit $?
fi
python3 ../ISRFH2StarLit/checkH2Attenuation.py twosource \
    --first "$run_name"/snap/snapshot_0001.hdf5 \
    --last "$(ls "$run_name"/snap/snapshot_*.hdf5 | tail -n 1)" \
    --fixture "$run_name"/fixture.json --params "$run_name"/used_parameters.yml --startup-log "$run_name"/output.log \
    --timesteps "$run_name"/timesteps.txt \
    --only-a only_a --only-b only_b --repeat-a only_a2 --slots "$slots"
