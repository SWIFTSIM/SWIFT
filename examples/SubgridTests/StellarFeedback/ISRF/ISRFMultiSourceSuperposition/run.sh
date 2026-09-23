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

# Set to 1 to run the lattice below the kernel-support bar on purpose. The
# pre-flight then warns instead of refusing, and every lattice run's
# output.log opens with that warning. The upper bound is not overridable.
study_below_bar=${ISRF_SUPERPOSITION_STUDY_BELOW_KERNEL_SUPPORT_BAR:=0}

swift=$(realpath "$swift")

# The lattice run holds its gates only while the star injection kernel
# reaches the neighbouring sources without covering the lattice. The
# support is gamma eta times the gas interparticle spacing, the source
# spacing is the box over n_side. Both bounds come from the check itself.
kernel_support_warning=""
case " $runs " in
    *" lattice "*)
        eta=$(sed -n 's/^ *resolution_eta: *\([0-9.eE+-]*\).*/\1/p' params.yml)
        if [ -z "$eta" ]; then
            echo "Cannot read SPH:resolution_eta from params.yml" >&2
            exit 1
        fi
        warning_file=$(mktemp)
        status=0
        python3 - "$eta" "$n_side" "$level" "$study_below_bar" \
                "$warning_file" <<'EOF' || status=$?
import sys

sys.path.insert(0, ".")
from isrf_multi_source_superposition_check import (
    GAMMA_3D,
    KERNEL_SUPPORT_OVER_SPACING_BAR as bar,
    KERNEL_SUPPORT_OVER_SPACING_COVERED as covered,
)

eta, n_side, level = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
study_below_bar, warning_file = sys.argv[4] == "1", sys.argv[5]
ratio = GAMMA_3D * eta * n_side / 2 ** level
print(
    f"Star kernel support over source spacing: {ratio:.3f} "
    f"(>= {bar:.2f}, < {covered:.3f})"
)
if not ratio < covered:
    sys.exit(
        f"Refusing to run: at level {level} with n_side {n_side} the star "
        f"injection kernel spans {ratio:.3f} of the source spacing, at or "
        f"above the {covered:.3f} that covers a cubic lattice. Every "
        "particle then sits inside a kernel and the isotropic-closure gate "
        "has no gas to measure. Lower n_side, or raise the resolution "
        "level."
    )
if not ratio >= bar:
    refusal = (
        f"at level {level} with n_side {n_side} the star injection kernel "
        f"spans {ratio:.3f} of the source spacing, under {bar:.2f}, a guard "
        "under the 0.598 this example is measured good at. At 0.299 the "
        "superposed field keeps its mean while its spatial variance exceeds "
        "the continuum lattice sum by an order of magnitude, with no fix "
        "available; the ratio at which that turns over has not been "
        "measured. Raise n_side, or lower the resolution level, until the "
        "kernel reaches the neighbouring sources."
    )
    if not study_below_bar:
        sys.exit(f"Refusing to run: {refusal}")
    # run.sh repeats this in every lattice run's own log.
    with open(warning_file, "w") as f:
        f.write(
            "WARNING: ISRF_SUPERPOSITION_STUDY_BELOW_KERNEL_SUPPORT_BAR is "
            "set, so this run is OUTSIDE the regime this example is "
            f"validated in: {refusal} Its gates do not hold here, and the "
            "check's L3 gate fails on the same bound.\n"
        )
EOF
        kernel_support_warning=$(cat "$warning_file")
        rm -f "$warning_file"
        [ $status -eq 0 ] || exit 1
        if [ -n "$kernel_support_warning" ]; then
            printf '%s\n' "$kernel_support_warning" >&2
        fi
        ;;
esac

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

    (cd "$run" && {
        case $run in
            lattice*)
                if [ -n "$kernel_support_warning" ]; then
                    printf '%s\n' "$kernel_support_warning"
                fi ;;
        esac
        "$swift" --hydro --stars --external-gravity --feedback \
        --cooling --sync --limiter --verbose=0 --threads=$n_threads \
        -P InitialConditions:file_name:ICs_isrf_multi_source.hdf5 \
        -P GrackleCooling:cloudy_table:../CloudyData_UVB=HM2012.h5 \
        -P GEARFeedback:yields_table:../PopII_parsec_spectral.hdf5 \
        -P GEARFeedback:yields_table_first_stars:../PopII_parsec_spectral.hdf5 \
        -P TimeIntegration:time_end:$t_end \
        -P TimeIntegration:dt_max:$dt \
        -P Snapshots:delta_time:$dsnap \
        -P Statistics:delta_time:$dsnap \
        -P GEARChemistry:initial_metallicity:$Z \
        -P GEARFeedback:ISRF_propagation:$prop \
        ../params.yml
     } 2>&1 | tee output.log)
    grep -q "main: done. Bye." "$run/output.log"
done

if [ "$run_check" = 1 ]; then
    python3 isrf_multi_source_superposition_check.py
fi
