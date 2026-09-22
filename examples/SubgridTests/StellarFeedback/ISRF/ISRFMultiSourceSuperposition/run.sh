#!/bin/bash

# Shared GEAR example scripts (tables)
scripts_location="../../../../GEAR_ICs_and_SCRIPTS"

set -eo pipefail

n_threads=${n_threads:=8}                   # Number of threads
swift=${swift:="../../../../../swift"}      # SWIFT binary
run_check=${run_check:=1}                   # Run the check after the runs
level=${level:=5}                           # N = (2**level)**3 gas particles
gas_density=${gas_density:=1e3}             # atom/cm^3
gas_mass=${gas_mass:=0.1}                   # Msun
star_mass=${star_mass:=29.7}                # Msun
n_side=${n_side:=8}                         # Lattice stars per side
Z_pair=${Z_pair:=0.2}                       # Z/Zsun of the A, B, AB runs
Z_lattice=${Z_lattice:=0.1}                 # Z/Zsun of the lattice runs
time_end=${time_end:=1.28e-3}               # Propagation runs
dt_max=${dt_max:=2.6e-6}                    # Propagation runs; floored to the 2.5e-6 bin
delta_time=${delta_time:=8e-5}              # Propagation runs
sep_injection=${sep_injection:=0.0625}      # A to B distance over L, propagation off
sep_propagation=${sep_propagation:=0.25}    # A to B distance over L, propagation on
runs=${runs:="injection_A injection_B injection_AB A B AB lattice lattice_single"}

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
    # injection_*: propagation off, four steps; others: propagation on.
    case $run in
        injection_*) sources=${run#injection_}; Z=$Z_pair; prop=0
                     t_end=1e-5; dt=2.6e-6; dsnap=1e-5; sep=$sep_injection ;;
        lattice*)    sources=$run; Z=$Z_lattice; prop=1
                     t_end=$time_end; dt=$dt_max; dsnap=$delta_time; sep=$sep_propagation ;;
        *)           sources=$run; Z=$Z_pair; prop=1
                     t_end=$time_end; dt=$dt_max; dsnap=$delta_time; sep=$sep_propagation ;;
    esac

    echo "=== $run: sources=$sources Z=$Z propagation=$prop"
    rm -rf "$run"
    mkdir -p "$run/snap"
    python3 makeIC.py --level $level --rho $gas_density --mass $gas_mass \
        --star_mass $star_mass --sources $sources --n_side $n_side --separation $sep \
        -o "$run/ICs_isrf_multi_source.hdf5"

    (cd "$run" && "$swift" --hydro --stars --external-gravity --feedback \
        --cooling --sync --limiter --verbose=0 --threads=$n_threads \
        -P InitialConditions:file_name:ICs_isrf_multi_source.hdf5 \
        -P GrackleCooling:cloudy_table:../CloudyData_UVB=HM2012.h5 \
        -P GEARFeedback:yields_table:../POPIIsw.h5 \
        -P GEARFeedback:yields_table_first_stars:../POPIIsw.h5 \
        -P TimeIntegration:time_end:$t_end \
        -P TimeIntegration:dt_max:$dt \
        -P Snapshots:delta_time:$dsnap \
        -P Statistics:delta_time:$dsnap \
        -P GEARChemistry:initial_metallicity:$Z \
        -P GEARFeedback:ISRF_propagation:$prop \
        ../params.yml 2>&1 | tee output.log)
    grep -q "main: done. Bye." "$run/output.log"
done

if [ "$run_check" = 1 ]; then
    python3 isrf_multi_source_superposition_check.py
fi
