#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

set -eo pipefail

n_threads=${n_threads:=8}                   # Number of threads
swift=${swift:="../../../../../swift"}      # SWIFT binary
run_check=${run_check:=1}                   # Run the check after the runs
level=${level:=5}                           # N = (2**level)**3 gas particles
gas_density=${gas_density:=1e3}             # atom/cm^3
gas_mass=${gas_mass:=0.1}                   # Cold phase particle mass, Msun
T_cold=${T_cold:=500}                       # Cold phase temperature, K (hot: 16 T_cold)
star_mass=${star_mass:=29.7}                # Msun
c_hyp_pin=${c_hyp_pin:=5}                   # km/s
dt_fine=${dt_fine:=1.3e-6}                  # dt_max of the single-bin run; floored to the 1.25e-6 bin
runs=${runs:="conservation_hierarchy conservation_single_bin"}

swift=$(realpath "$swift")

glass_n=$((2**level))
if [ ! -e glassCube_${glass_n}.hdf5 ]; then
    ./getGlass.sh $glass_n
fi
if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
    "$scripts_location"/getGrackleCoolingTable.sh
fi
if [ ! -e POPIIsw.h5 ]; then
    "$scripts_location"/getChemistryTable.sh
fi

for run in $runs; do
    # dt_max 2e-5 leaves both phases on their CFL bins; dt_fine forces one bin.
    case $run in
        conservation_hierarchy)    dt=2e-5 ;;
        conservation_single_bin)   dt=$dt_fine ;;
        *) echo "Unknown run $run"; exit 1 ;;
    esac

    echo "=== $run: c_hyp_pin=$c_hyp_pin dt_max=$dt"
    rm -rf "$run"
    mkdir -p "$run/snap"
    python3 makeIC.py --level $level --rho $gas_density --mass $gas_mass \
        --T_cold $T_cold --star_mass $star_mass -o "$run/ICs_isrf_hierarchy.hdf5"

    (cd "$run" && "$swift" --hydro --stars --external-gravity --feedback \
        --cooling --sync --limiter --verbose=0 --threads=$n_threads \
        -P GrackleCooling:cloudy_table:../CloudyData_UVB=HM2012.h5 \
        -P GEARFeedback:yields_table:../POPIIsw.h5 \
        -P GEARFeedback:yields_table_first_stars:../POPIIsw.h5 \
        -P TimeIntegration:time_end:1.28e-3 \
        -P TimeIntegration:dt_max:$dt \
        -P Snapshots:delta_time:2e-5 \
        -P Statistics:delta_time:2e-5 \
        -P GEARChemistry:initial_metallicity:0 \
        -P GEARFeedback:ISRF_c_hyp_pin_for_debugging:$c_hyp_pin \
        ../params.yml 2>&1 | tee output.log)
    grep -q "main: done. Bye." "$run/output.log"
done

if [ "$run_check" = 1 ]; then
    python3 isrf_timestep_hierarchy_injection_check.py
fi
