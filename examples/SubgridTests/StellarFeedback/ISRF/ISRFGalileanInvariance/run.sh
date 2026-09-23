#!/bin/bash

# Shared GEAR example scripts (tables, plotting)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

set -eo pipefail

config=${config:="galilean"}           # galilean or lag (see README)
n_threads=${n_threads:=8}
gas_density=${gas_density:=1e3}        # atom/cm^3
gas_particle_mass=${gas_mass:=0.1}     # Msun
star_mass=${star_mass:=29.7}           # Msun
initial_metallicity=${initial_metallicity:=1}

if [ "$config" = "galilean" ]; then
    level=${level:=5}
    time_end=${time_end:=5e-5}
    dt_max=${dt_max:=2.5e-6}
    delta_time=${delta_time:=1e-5}
    c_hyp_pin=${c_hyp_pin:=0}          # 0: shipped closure c_hyp = C_hyp*h/dt
    alpha_max=${alpha_max:=0.5}
    alpha_floor=${alpha_floor:=0.5}
    boost_factors=${boost_factors:="0.1 1 10"}  # bulk velocity in units of the measured c_hyp
    boost_kms=${boost_kms:="30"}                # bulk velocities in km/s
elif [ "$config" = "lag" ]; then
    level=${level:=6}                  # same h as level 5, twice the box
    time_end=${time_end:=5e-4}
    dt_max=${dt_max:=2.5e-6}
    delta_time=${delta_time:=5e-5}
    c_hyp_pin=${c_hyp_pin:=10}         # km/s; a uniform c_hyp gives a uniform tau
    alpha_max=${alpha_max:=0}
    alpha_floor=${alpha_floor:=0}
    v_rel_factors=${v_rel_factors:="0.25 0.5 1"}  # star velocity in units of c_hyp_pin
else
    echo "Unknown config=$config (galilean or lag)" >&2
    exit 1
fi

glass_n=$((2**level))
if [ ! -e "glassCube_${glass_n}.hdf5" ]; then
    echo "Fetching initial glass file (glassCube_${glass_n}.hdf5)..."
    ./getGlass.sh $glass_n
fi
if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
    "$scripts_location"/getGrackleCoolingTable.sh
fi
"$scripts_location"/getRadiationTable.sh PopII_parsec_spectral.hdf5 || exit 1

# Stop here on a table the radiation reader cannot use, rather than
# aborting at start-up once the initial conditions are built.
"$scripts_location"/checkRadiationTable.sh PopII_parsec_spectral.hdf5 --with-isrf || exit 1

# run_one <run_dir> <bulk vx> <star vx> [position shift]
run_one() {
    local dir=$1 vbulk=$2 vstar=$3 shift=${4:-"[0.0,0.0,0.0]"}
    rm -rf "$dir" snap
    mkdir -p snap
    python3 makeIC.py --level $level --rho $gas_density --mass $gas_particle_mass \
        --star_mass $star_mass --bulk-velocity $vbulk 0 0 --star-velocity $vstar 0 0 \
        -o ICs_isrf_galilean_invariance.hdf5
    ../../../../../swift --hydro --stars --external-gravity --feedback --cooling \
        --sync --limiter --verbose=0 --threads=$n_threads \
        -P TimeIntegration:time_end:$time_end \
        -P TimeIntegration:dt_max:$dt_max \
        -P Snapshots:delta_time:$delta_time \
        -P GEARChemistry:initial_metallicity:$initial_metallicity \
        -P GEARFeedback:ISRF_c_hyp_pin_for_debugging:$c_hyp_pin \
        -P GEARFeedback:ISRF_dissipation_alpha_max:$alpha_max \
        -P GEARFeedback:ISRF_dissipation_alpha_floor:$alpha_floor \
        -P "InitialConditions:shift:$shift" \
        params.yml 2>&1 | tee output.log
    grep -q "main: done. Bye." output.log
    mkdir -p "$dir"
    mv snap output.log timesteps.txt statistics.txt used_parameters.yml \
        unused_parameters.yml ICs_isrf_galilean_invariance.hdf5 "$dir"
}

if [ "$config" = "galilean" ]; then
    run_one runs/rest 0 0
    run_one runs/rest_repeat 0 0
    # Same physics, particles assigned to different cells: the round-off floor of the cell layout.
    shift=$(python3 -c "import h5py; L=h5py.File('runs/rest/ICs_isrf_galilean_invariance.hdf5')['Header'].attrs['BoxSize'][0]; print(f'[{0.37*L},{0.13*L},{0.06*L}]')")
    run_one runs/rest_shifted 0 0 "$shift"
    c_hyp=$(python3 galilean_compare.py --measure-c-hyp runs/rest)
    echo "Measured median c_hyp in the rest run: $c_hyp km/s"
    boosted=()
    for k in $boost_factors; do
        v=$(python3 -c "print($k * $c_hyp)")
        run_one runs/boost_${k}c $v 0
        boosted+=(runs/boost_${k}c)
    done
    for v in $boost_kms; do
        run_one runs/boost_${v}kms $v 0
        boosted+=(runs/boost_${v}kms)
    done
    python3 galilean_compare.py --rest runs/rest --repeat runs/rest_repeat --shifted runs/rest_shifted \
        --boosted "${boosted[@]}"
else
    status=0
    for k in $v_rel_factors; do
        v=$(python3 -c "print($k * $c_hyp_pin)")
        run_one runs/lag_${k}c 0 $v
        gate_flag="--gate"
        # Numeric, not string, comparison: alpha_max=0.0 must match "0".
        if ! python3 -c "import sys
sys.exit(0 if float(sys.argv[1]) == 0.0 and float(sys.argv[2]) == 0.0 else 1)" \
            "$alpha_max" "$alpha_floor"; then
            gate_flag="--report-only"
        fi
        python3 lag_check.py --run runs/lag_${k}c $gate_flag || status=1
    done
    exit $status
fi
