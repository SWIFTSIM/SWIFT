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
# The cadence sweep: runs="cadence_base cadence_half cadence_quarter".

swift=$(realpath "$swift")

glass_n=$((2**level))
if [ ! -e glassCube_${glass_n}.hdf5 ]; then
    ./getGlass.sh $glass_n
fi
if [ ! -e CloudyData_UVB=HM2012.h5 ]; then
    "$scripts_location"/getGrackleCoolingTable.sh
fi
"$scripts_location"/getRadiationTable.sh PopII_parsec_spectral.hdf5 || exit 1

# Stop here on a table the radiation reader cannot use, rather than
# aborting at start-up once the initial conditions are built.
"$scripts_location"/checkRadiationTable.sh PopII_parsec_spectral.hdf5 --with-isrf || exit 1

for run in $runs; do
    # dt_max 2e-5 leaves both phases on their CFL bins; dt_fine forces one
    # bin. The cadence sweep halves dt_fine twice, so every gas particle and
    # the star share one bin in all three and the bin itself is the only
    # thing that moves.
    case $run in
        conservation_hierarchy)    dt=2e-5 ;;
        conservation_single_bin)   dt=$dt_fine ;;
        cadence_base)              dt=$dt_fine ;;
        cadence_half)              dt=$(python3 -c "print(0.5*$dt_fine)") ;;
        cadence_quarter)           dt=$(python3 -c "print(0.25*$dt_fine)") ;;
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
        -P GEARFeedback:yields_table:../PopII_parsec_spectral.hdf5 \
        -P GEARFeedback:yields_table_first_stars:../PopII_parsec_spectral.hdf5 \
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
    case " $runs " in
        *" conservation_"*) python3 isrf_timestep_hierarchy_injection_check.py ;;
    esac
    case " $runs " in
        *" cadence_"*) python3 isrf_injection_cadence_independence_check.py ;;
    esac
fi
