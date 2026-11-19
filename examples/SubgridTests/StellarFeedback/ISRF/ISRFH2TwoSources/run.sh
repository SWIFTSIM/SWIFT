#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

# swift's exit status is piped into tee: without pipefail a refused or crashed run looks complete.
set -eo pipefail

leg=${leg:="main"}   # main, switch_off, stale_key
n_threads=${n_threads:=8}   # Number of threads
ranks=${ranks:=1}   # MPI ranks, 2 needs an MPI build named swift_mpi
level=${level:=5}   # Resolution level: N = (2**level)**3 gas particles
gas_density=${gas_density:=1e3}   # Gas density in atom/cm^3
gas_mass=${gas_mass:=0.1}   # Gas particle mass in Msun
h2_ratio=${h2_ratio:=3e-4}   # Initial nH2 / nH
run_name=${run_name:=$leg}
ic_file=ICs_isrf_h2_two_sources.hdf5
params=params.yml
extra_ic_args="${contrast:+--contrast $contrast}"
extra_swift_args=""
stale_key_arg=""

case "$leg" in
    main) switch=1 ;;
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
if [ "$leg" = "main" ]; then
    run_stars a only_a
    run_stars b only_b
    run_stars a only_a2
fi
run_stars ab "$run_name"

slots=$(sed -n 's/^#define RADIATION_PACKET_SLOTS //p' "$(dirname "$swift_exe")/config.h")
[ -n "$slots" ] || { echo "RADIATION_PACKET_SLOTS not found in $(dirname "$swift_exe")/config.h"; exit 1; }
python3 ../ISRFH2StarLit/checkH2Attenuation.py twosource \
    --first "$run_name"/snap/snapshot_0001.hdf5 \
    --last "$(ls "$run_name"/snap/snapshot_*.hdf5 | tail -n 1)" \
    --fixture "$run_name"/fixture.json --params "$run_name"/used_parameters.yml --startup-log "$run_name"/output.log \
    --only-a only_a --only-b only_b --repeat-a only_a2 --slots "$slots"
