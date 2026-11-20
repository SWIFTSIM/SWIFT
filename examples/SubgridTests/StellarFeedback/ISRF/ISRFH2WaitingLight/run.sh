#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

# swift's exit status is piped into tee: without pipefail a refused or crashed run looks complete.
set -eo pipefail

n_threads=${n_threads:=8}   # Number of threads
level=${level:=5}   # Resolution level: N = (2**level)**3 gas particles
gas_density=${gas_density:=1e3}   # Gas density in atom/cm^3
gas_mass=${gas_mass:=0.1}   # Gas particle mass in Msun
h2_ratio=${h2_ratio:=3e-5}   # Initial nH2 / nH
light_speed=${light_speed:=fixed}   # Light speed of the limited arm: fixed or cap
ic_file=ICs_isrf_h2_waiting_light.hdf5

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
python3 makeIC.py --level $level --rho $gas_density --mass $gas_mass --h2_ratio $h2_ratio -o "$ic_file"

swift_exe=../../../../../swift

# Two arms from the same ICs: the control arm has unlimited light speed (no light waits) and the
# same time step as the limited arm. Each arm writes the same snapshot times.
run_arm() {
    local arm=$1; shift
    rm -rf snap restart
    mkdir snap
    set +e
    $swift_exe --hydro --stars --external-gravity --feedback --cooling \
        --sync --limiter --verbose=0 --threads=$n_threads \
        -P GEARFeedback:ISRF_H2_line_absorption:1 \
        -P GeajayCooling:initial_h2i_to_h_ratio:$h2_ratio \
        "$@" params.yml 2>&1 | tee output_$arm.log
    local status=${PIPESTATUS[0]}
    set -e
    [ $status -eq 0 ] || { echo "swift stopped with status $status"; exit $status; }
    mkdir -p "$arm"
    cp used_parameters.yml timesteps.txt fixture.json "$arm"/
    cp output_$arm.log "$arm"/output.log
    mv snap "$arm"/
}

run_arm control -P GEARFeedback:ISRF_packet_light_speed:unlimited -P TimeIntegration:dt_max:2.34375e-10
run_arm limited -P GEARFeedback:ISRF_packet_light_speed:$light_speed

python3 checkH2WaitingLight.py --limited limited --control control
