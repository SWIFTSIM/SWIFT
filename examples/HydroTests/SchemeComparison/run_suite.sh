#!/bin/bash
#
# Build several hydro flavours of SWIFT and run the standard suite of hydro
# tests with each of them, with identical initial conditions, then produce
# comparison plots and a summary table.
#
# Usage: ./run_suite.sh [options]
#   -s DIR    SWIFT source tree (default: the tree containing this script)
#   -o DIR    output directory (default: ./suite_output)
#   -c FILE   scheme list, "<name> <configure flags> [| <run-time -P overrides>]"
#             per line (default: schemes.txt next to this script)
#   -t LIST   comma-separated tests to run (default: all, see TESTS below)
#   -j N      number of threads for the runs and the builds (default: 8)
#   -k NAME   kernel passed to --with-kernel (default: wendland-C4)
#   -b        (re)build the binaries (done automatically when missing)
#   -n        do not run, only build and/or plot
#   -p        do not plot
#   -h        this help
#
# Layout of the output directory:
#   builds/<scheme>_{2d,3d}/   source copies + binaries
#   ics/<test>/                initial conditions (shared by all schemes)
#   runs/<scheme>/<test>/      snapshots, logs, statistics, example plots
#   plots/                     comparison figures and summary.md
#
# Every test directory contains the example's own plotSolution.py output and
# the SWIFT output.log, so that the standard views are available too.

set -u

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SRC_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
OUT_DIR=$PWD/suite_output
SCHEMES_FILE=$SCRIPT_DIR/schemes.txt
THREADS=8
KERNEL=wendland-C4
DO_BUILD=0
DO_RUN=1
DO_PLOT=1

# ---------------------------------------------------------------------------
# Knobs for the individual tests (resolution / duration). Change them here or
# export them before calling the script.
# ---------------------------------------------------------------------------
: "${KH_L2:=128}"               # KH: particles per edge in the low-density region
: "${EVRARD_NPARTS:=100000}"    # Evrard: number of particles
: "${KEPLERIAN_T_END:=50}"      # Keplerian ring: end time
: "${ZELDOVICH_Z_END:=0.9}"     # Zeldovich: final redshift
: "${ZELDOVICH_PERTURB:=0.1}"   # Zeldovich_pert: transverse random displacement in units of the spacing
: "${SEDOV_GLASS:=glassCube_64}"
: "${NOH_GLASS:=glassCube_64}"

GLASS_URL=https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs
REF_URL=https://virgodb.cosma.dur.ac.uk/swift-webstorage/ReferenceSolutions

ALL_TESTS="gresho square zeldovich sod keplerian kh noh evrard sedov zeldovich_pert"
TESTS=$ALL_TESTS

while getopts "s:o:c:t:j:k:bnph" opt; do
  case $opt in
    s) SRC_DIR=$(cd "$OPTARG" && pwd) ;;
    o) OUT_DIR=$OPTARG ;;
    c) SCHEMES_FILE=$OPTARG ;;
    t) TESTS=$(echo "$OPTARG" | tr ',' ' ') ;;
    j) THREADS=$OPTARG ;;
    k) KERNEL=$OPTARG ;;
    b) DO_BUILD=1 ;;
    n) DO_RUN=0 ;;
    p) DO_PLOT=0 ;;
    h) sed -n 2,30p "$0"; exit 0 ;;
    *) exit 1 ;;
  esac
done

mkdir -p "$OUT_DIR" || exit 1
OUT_DIR=$(cd "$OUT_DIR" && pwd)
EXAMPLES=$SRC_DIR/examples
# The softened point mass is what the Keplerian ring ICs assume (see its README)
COMMON_FLAGS="--with-kernel=$KERNEL --with-ext-potential=point-mass-softened --disable-doxygen-doc --disable-mpi"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# Scheme names, configure flags and (optional, after a '|') run-time
# parameter overrides passed to every run of the scheme, e.g.
#   magma_eta15 --with-hydro=magma2 | -P SPH:resolution_eta:1.5
SCHEME_NAMES=()
declare -A SCHEME_FLAGS
declare -A SCHEME_RUNPARAMS
while read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  runparams=""
  if [[ "$line" == *"|"* ]]; then runparams="${line#*|}"; line="${line%%|*}"; fi
  read -r name flags <<< "$line"
  [[ -z "$name" ]] && continue
  SCHEME_NAMES+=("$name")
  SCHEME_FLAGS[$name]=$flags
  SCHEME_RUNPARAMS[$name]=$runparams
done < "$SCHEMES_FILE"
[[ ${#SCHEME_NAMES[@]} -eq 0 ]] && { echo "No scheme in $SCHEMES_FILE"; exit 1; }

# ---------------------------------------------------------------------------
# Test definitions
# ---------------------------------------------------------------------------
# test_<name> sets: DIM, EXAMPLE (relative to examples/), YML, IC (file name),
# FLAGS (SWIFT run-time flags), PARAMS (-P overrides), PLOT_SNAP (snapshot
# passed to the example's plotSolution.py; empty: none), and defines make_ic().
test_sod() {
  DIM=3; EXAMPLE=HydroTests/SodShock_3D; YML=sodShock.yml; IC=sodShock.hdf5
  FLAGS="--hydro"; PARAMS=""; PLOT_SNAP=1
  make_ic() { get_glass glassCube_64 glassCube_32; python3 makeIC.py; }
}
test_sedov() {
  DIM=3; EXAMPLE=HydroTests/SedovBlast_3D; YML=sedov.yml; IC=sedov.hdf5
  FLAGS="--hydro --limiter"; PARAMS=""; PLOT_SNAP=5
  make_ic() {
    get_glass "$SEDOV_GLASS"
    sed -i "s/glassCube_64.hdf5/$SEDOV_GLASS.hdf5/" makeIC.py
    python3 makeIC.py
  }
}
test_noh() {
  DIM=3; EXAMPLE=HydroTests/Noh_3D; YML=noh.yml; IC=noh.hdf5
  FLAGS="--hydro"; PARAMS=""; PLOT_SNAP=12
  make_ic() {
    get_glass "$NOH_GLASS"
    sed -i "s/glassCube_64.hdf5/$NOH_GLASS.hdf5/" makeIC.py
    python3 makeIC.py
  }
}
test_evrard() {
  DIM=3; EXAMPLE=HydroTests/EvrardCollapse_3D; YML=evrard.yml; IC=evrard.hdf5
  FLAGS="--hydro --self-gravity"; PARAMS=""; PLOT_SNAP=8
  make_ic() {
    python3 makeIC.py -n "$EVRARD_NPARTS"
    [[ -e evrardCollapse3D_exact.txt ]] || wget -q $REF_URL/evrardCollapse3D_exact.txt
  }
}
test_kh() {
  DIM=2; EXAMPLE=HydroTests/KelvinHelmholtz_2D; YML=kelvinHelmholtz.yml; IC=kelvinHelmholtz.hdf5
  FLAGS="--hydro"; PARAMS="-P Snapshots:delta_time:0.5"; PLOT_SNAP=9
  make_ic() { sed -i "s/^L2 = .*/L2 = $KH_L2/" makeIC.py; python3 makeIC.py; }
}
test_square() {
  DIM=2; EXAMPLE=HydroTests/SquareTest_2D; YML=square.yml; IC=square.hdf5
  FLAGS="--hydro"; PARAMS="-P Snapshots:delta_time:0.5"; PLOT_SNAP=8
  make_ic() { python3 makeIC.py; }
}
test_gresho() {
  DIM=2; EXAMPLE=HydroTests/GreshoVortex_2D; YML=gresho.yml; IC=greshoVortex.hdf5
  FLAGS="--hydro"; PARAMS=""; PLOT_SNAP=10
  make_ic() { get_glass glassPlane_128; python3 makeIC.py; }
}
test_keplerian() {
  # A razor-thin disc run with the 3D code (see the example's README)
  DIM=3; EXAMPLE=HydroTests/KeplerianRing; YML=keplerian_ring.yml; IC=initial_conditions.hdf5
  FLAGS="--hydro --external-gravity"
  PARAMS="-P TimeIntegration:time_end:$KEPLERIAN_T_END -P Snapshots:delta_time:1"
  PLOT_SNAP=""
  make_ic() { python3 makeIC.py; }
}
test_zeldovich() {
  DIM=3; EXAMPLE=Cosmology/ZeldovichPancake_3D; YML=zeldovichPancake.yml; IC=zeldovichPancake.hdf5
  FLAGS="--hydro --self-gravity --cosmology"
  PARAMS="-P Cosmology:a_end:$(python3 -c "print(1./(1.+$ZELDOVICH_Z_END))")"
  PLOT_SNAP=""
  make_ic() { python3 makeIC.py; }
}

# Same as zeldovich, but the lattice symmetry is broken by random transverse
# displacements of the particles (the 1D flow is unchanged).
test_zeldovich_pert() {
  test_zeldovich
  make_ic() {
    python3 - <<PYEOF
s = open("makeIC.py").read()
s = s.replace("# Unit conversion", """# Break the lattice symmetry (transverse random displacements)
_rng = numpy.random.default_rng(42)
coords[:, 1:] += $ZELDOVICH_PERTURB * delta_x * (_rng.random((numPart, 2)) - 0.5)
coords[:, 1:] = numpy.mod(coords[:, 1:], boxSize)

# Unit conversion""")
s = s.replace("import h5py", "import h5py\nimport numpy")
open("makeIC.py", "w").write(s)
PYEOF
    python3 makeIC.py
  }
}

get_glass() {
  for g in "$@"; do
    [[ -e $g.hdf5 ]] || wget -q "$GLASS_URL/$g.hdf5" || { echo "Cannot fetch $g"; return 1; }
  done
}

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
build_one() {
  local name=$1 dim=$2
  local dir=$OUT_DIR/builds/${name}_${dim}d
  if [[ $DO_BUILD -eq 0 && -x $dir/swift ]]; then return 0; fi
  log "Building $name (${dim}D) in $dir"
  mkdir -p "$dir"
  # Copy the sources only: no build products of the source tree (in particular
  # not its swift binary, which would mask a failed build here).
  rsync -a --delete --exclude='.git' --exclude='*.hdf5' --exclude='*.o' \
    --exclude='*.lo' --exclude='*.la' --exclude='*.a' --exclude='.libs' \
    --exclude='/swift' --exclude='/swift_mpi' --exclude='/swift_fof' \
    --exclude='/swift_fof_mpi' --exclude='*.png' --exclude='*.pdf' \
    --exclude='restart' --exclude='suite_output' --exclude='results_*' \
    "$SRC_DIR/" "$dir/" || return 1
  local dimflag=""
  [[ $dim -eq 2 ]] && dimflag="--with-hydro-dimension=2"
  if ! ( cd "$dir" && ./autogen.sh && ./configure $COMMON_FLAGS $dimflag ${SCHEME_FLAGS[$name]} \
      && make -j"$THREADS" ) > "$dir/build.log" 2>&1; then
    echo "Build of ${name}_${dim}d failed, see $dir/build.log"; return 1
  fi
  [[ -x $dir/swift ]] || { echo "Build of ${name}_${dim}d produced no binary"; return 1; }
}

# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------
prepare_ic() {
  local test=$1
  test_$test
  local dir=$OUT_DIR/ics/$test
  if [[ -e $dir/$IC ]]; then return 0; fi
  log "Generating initial conditions for $test"
  mkdir -p "$dir"
  cp -r "$EXAMPLES/$EXAMPLE/." "$dir/"
  # Reuse glass files already downloaded for other tests
  for g in "$OUT_DIR"/ics/*/glass*.hdf5; do
    [[ -e $g && ! -e $dir/$(basename "$g") ]] && cp "$g" "$dir/"
  done
  ( cd "$dir" && make_ic ) > "$dir/makeIC.log" 2>&1 || { echo "IC generation for $test failed, see $dir/makeIC.log"; return 1; }
  [[ -e $dir/$IC ]] || { echo "IC file $IC missing for $test"; return 1; }
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
run_one() {
  local scheme=$1 test=$2
  test_$test
  local icdir=$OUT_DIR/ics/$test
  local dir=$OUT_DIR/runs/$scheme/$test
  local swift=$OUT_DIR/builds/${scheme}_${DIM}d/swift
  if [[ -e $dir/DONE ]]; then log "$scheme/$test already done"; return 0; fi
  log "Running $test with $scheme"
  rm -rf "$dir"; mkdir -p "$dir"
  # Everything the example ships with (plot scripts, reference solutions, ...)
  for f in "$icdir"/*; do
    case $(basename "$f") in *.hdf5) ln -s "$f" "$dir/" ;; *) cp -r "$f" "$dir/" ;; esac
  done
  local start=$(date +%s)
  ( cd "$dir" && "$swift" $FLAGS --threads="$THREADS" $PARAMS ${SCHEME_RUNPARAMS[$scheme]} "$YML" > output.log 2>&1 )
  local status=$?
  local end=$(date +%s)
  echo "$((end - start))" > "$dir/walltime_s"
  if [[ $status -ne 0 ]]; then echo "$scheme/$test FAILED (exit $status), see $dir/output.log"; return 1; fi
  # Standard plot of the example (best effort)
  if [[ -n "$PLOT_SNAP" && -e $dir/plotSolution.py ]]; then
    ( cd "$dir" && sed -i "s|\.\./\.\./\.\./tools/stylesheets|$SRC_DIR/tools/stylesheets|" plotSolution.py \
        && python3 plotSolution.py "$PLOT_SNAP" ) > "$dir/plot.log" 2>&1
  fi
  touch "$dir/DONE"
  log "$scheme/$test finished in $((end - start)) s"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
for s in "${SCHEME_NAMES[@]}"; do
  for d in 2 3; do build_one "$s" "$d" || exit 1; done
done

if [[ $DO_RUN -eq 1 ]]; then
  for t in $TESTS; do prepare_ic "$t" || exit 1; done
  for t in $TESTS; do
    for s in "${SCHEME_NAMES[@]}"; do run_one "$s" "$t"; done
  done
fi

if [[ $DO_PLOT -eq 1 ]]; then
  log "Plotting"
  python3 "$SCRIPT_DIR/plot_comparison.py" "$OUT_DIR" $TESTS
fi
log "All done. Results in $OUT_DIR"
