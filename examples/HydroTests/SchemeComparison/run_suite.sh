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
#   -i        prepare the initial conditions only (no build, run or plot)
#   -p        do not plot
#   -h        this help
#
# Schemes with identical configure flags share one build: the build
# directory of a scheme is a symbolic link to builds/by_flags/<hash>_<dim>,
# so that a parameter study with many run-time variants of one scheme
# compiles the code once per flavour.
#
# Layout of the output directory:
#   builds/<scheme>_{1d,2d,2drt,3d,3dnfw}/   source copies + binaries
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
ICS_ONLY=0

# ---------------------------------------------------------------------------
# Knobs for the individual tests (resolution / duration). Change them here or
# export them before calling the script.
# ---------------------------------------------------------------------------
: "${KH_L2:=128}"               # KH: particles per edge in the low-density region
: "${EVRARD_NPARTS:=100000}"    # Evrard: number of particles
: "${KEPLERIAN_T_END:=50}"      # Keplerian ring: end time
: "${KEPLERIAN2D_T_END:=30}"    # Keplerian ring, 2D code: end time (the 2D run segfaults later on particles without neighbours)
: "${KEPLERIAN_N:=64}"          # Keplerian disc (Hopkins): effective resolution across the disc (makeIC.py -n; Hopkins: 256)
: "${KEPLERIAN_DT_MIN:=1e-9}"   # Keplerian disc/ring: minimal time step (particles plunging onto the point mass need it)
: "${KEPLERIAN_CD10_T_END:=40}" # Keplerian ring (Cullen & Dehnen): end time (paper: 40 = 6.4 orbits at the ring centre)
: "${KEPLERIAN_CD10_N:=9745}"   # Keplerian ring (Cullen & Dehnen): number of particles (paper: 9745)
: "${ZELDOVICH_Z_END:=0.9}"     # Zeldovich: final redshift
: "${ZELDOVICH_PERTURB:=0.1}"   # Zeldovich_pert: transverse random displacement in units of the spacing
: "${ZELDOVICH_N:=32}"          # Zeldovich (lattice): particles per dimension
: "${ZELDOVICH_GLASS:=glassCube_32}" # Zeldovich_glass: glass file for the Lagrangian positions
: "${ZELDOVICH_Z_START:=10}"    # Zeldovich_glass: starting redshift (limits the growth of the glass noise)
: "${SEDOV_GLASS:=glassCube_64}"
: "${NOH_GLASS:=glassCube_64}"
: "${BLOB_N:=64}"               # Blob: makeIC.py's num_on_side (BCC lattice: N^3/4 particles per unit cube, N even)
: "${BLOB_T_END:=1.75}"         # Blob: end time (tau_KH = 0.348 for the example's set-up: wind speed 1, Mach 2.7)
: "${NFW_T_END:=1.0}"           # NFW hydrostatic halo: end time (internal units, ~Gyr)
: "${TG_T_END:=5.0}"            # Taylor-Green vortex: end time (eddy turnover 1.6)
: "${SHEARING_N:=4}"            # Shearing shock: glass cubes along x
: "${SHEARING_S:=10}"           # Shearing shock: shear amplitude in units of the stream speed (paper: 100, which heats the box through the diffused shear by 10-100x in every SPH scheme)
: "${SHEARING_T_END:=0.6}"      # Shearing shock: end time (shocks and edge rarefactions meet at t ~ 1 for 4 cubes)
: "${BLAST_N:=4}"               # Strong blast wave: glass cubes along x
: "${BLAST_T_END:=0.02}"        # Strong blast wave: end time (the waves from the periodic boundary stay out of the middle half of the box for 4 cubes)
: "${RT_N:=128}"                # Rayleigh-Taylor: particles along x (the box is 1 x 1.5)
: "${RT_T_END:=4.0}"            # Rayleigh-Taylor: end time (Hopkins 2015, Rosswog 2020: 4)
: "${KHGR_T_END:=2.0}"          # KH growth rate: end time (the seeded mode saturates at t ~ 1.2)
: "${KHGR_SEED:=0.01}"          # KH growth rate: amplitude of the seed velocity (the example's 5e-4 is below the shear noise of SPH)
: "${KERNEL_1D:=wendland-C2}"   # Kernel of the 1D build (the Wendland C4 and C6 kernels have no 1D coefficients in SWIFT)
: "${SOUNDWAVE_GLASS:=glassCube_32}" # Sound wave: 3D glass (32: 32 particles per wavelength for one wave per box)
: "${SOUNDWAVE_AMPLITUDE:=1e-3}" # Sound wave: relative density amplitude (linear for ~100 periods)
: "${SOUNDWAVE_NWAVE:=1}"       # Sound wave: wavelengths across the box
: "${SOUNDWAVE_T_END:=5}"       # Sound wave (travelling): end time in periods of the n = 1 wave
: "${WAVEDAMPING_T_END:=10}"    # Wave damping (standing wave): end time
: "${THERMAL_CONTRAST:=4}"      # Thermal discontinuity: density contrast of the cold slab
: "${THERMAL_T_END:=3}"         # Thermal discontinuity: end time (sound crossing of the slab: 1.5)
: "${SOD_WEAK_P_RIGHT:=0.5}"    # Sod series: right-state pressure of the weak shock (standard 0.1)
: "${SOD_STRONG_P_RIGHT:=0.01}" # Sod series: right-state pressure of the strong shock
: "${KH8_CONTRAST:=8}"         # High-contrast KH: density contrast (the example's IC is a sharp interface of equal-mass particles)
: "${KH8_L2:=128}"              # High-contrast KH: particles per edge in the low-density region
: "${KH8_T_END:=6}"             # High-contrast KH: end time (the KH time of the seeded mode is 1.6 at contrast 8, 1.1 at contrast 2)
: "${GRESHO_BOOST:=3}"          # Boosted Gresho vortex: bulk velocity along x (sound speed 2.9-3.1; integer: the vortex is back at the centre at t = 1)

GLASS_URL=https://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs
REF_URL=https://virgodb.cosma.dur.ac.uk/swift-webstorage/ReferenceSolutions

# The NFW halo (nfw) is not in the default list: expensive and so far not
# discriminating between schemes; run it with -t nfw.
ALL_TESTS="gresho square square_mass zeldovich zeldovich_glass sod keplerian keplerian2d keplerian_cd10 kh noh evrard sedov zeldovich_pert blob taylorgreen shearing blast rt khgr ibw soundwave wavedamping thermal gresho_boost kh8 sod_weak sod_strong"
TESTS=$ALL_TESTS

while getopts "s:o:c:t:j:k:bniph" opt; do
  case $opt in
    s) SRC_DIR=$(cd "$OPTARG" && pwd) ;;
    o) OUT_DIR=$OPTARG ;;
    c) SCHEMES_FILE=$OPTARG ;;
    t) TESTS=$(echo "$OPTARG" | tr ',' ' ') ;;
    j) THREADS=$OPTARG ;;
    k) KERNEL=$OPTARG ;;
    b) DO_BUILD=1 ;;
    n) DO_RUN=0 ;;
    i) ICS_ONLY=1 ;;
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
# Tests needing a build with different common flags set BUILD (default: the
# "<dim>d" build of the scheme), see BUILD_FLAGS below.
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
# The square with the paper's set-up (Rosswog 2020, Saitoh & Makino 2013):
# equal-spacing lattice with 4:1 particle MASSES instead of 4:1 spacing. The
# example sets u from the nominal densities, so the pressure is not exactly
# uniform across the kernel-smoothed contact at t = 0; the blip relaxes within
# a sound crossing of the interface and does not change the t = 4 result
# (checked against u set from the SPH density: L1 0.017 vs 0.024 for MI2).
test_square_mass() {
  DIM=2; EXAMPLE=HydroTests/SquareTest_2D; YML=square.yml; IC=square.hdf5
  FLAGS="--hydro"; PARAMS="-P Snapshots:delta_time:0.5"; PLOT_SNAP=8
  make_ic() { python3 makeICDifferentMasses.py; }
}
test_gresho() {
  DIM=2; EXAMPLE=HydroTests/GreshoVortex_2D; YML=gresho.yml; IC=greshoVortex.hdf5
  FLAGS="--hydro"; PARAMS=""; PLOT_SNAP=10
  make_ic() { get_glass glassPlane_128; python3 makeIC.py; }
}
test_keplerian() {
  # Hopkins (2015) Keplerian disc as a razor-thin disc run with the 3D code (see the example's README)
  DIM=3; EXAMPLE=HydroTests/KeplerianRing; YML=keplerian_ring.yml; IC=initial_conditions.hdf5
  FLAGS="--hydro --external-gravity"
  PARAMS="-P TimeIntegration:time_end:$KEPLERIAN_T_END -P Snapshots:delta_time:1 -P TimeIntegration:dt_min:$KEPLERIAN_DT_MIN"
  PLOT_SNAP=""
  make_ic() { python3 makeIC.py -n "$KEPLERIAN_N"; }
}
# The same ring with the 2D code: the particles are put in the z = 0 plane and
# the point mass with them. In the 3D planar set-up every MAGMA particle has
# a singular correction matrix (no z extent), which exercises the regularised
# inverse; in 2D the 2x2 matrices are regular, so this is the reference that
# tests the matrix-inversion gradients themselves on a differentially rotating
# flow. Late-time (t > 10) ring metrics scatter between runs of the same code.
test_keplerian2d() {
  DIM=2; EXAMPLE=HydroTests/KeplerianRing; YML=keplerian_ring.yml; IC=initial_conditions.hdf5
  FLAGS="--hydro --external-gravity"
  PARAMS="-P TimeIntegration:time_end:$KEPLERIAN2D_T_END -P Snapshots:delta_time:1 -P PointMassPotential:position:[5.,5.,0.] -P TimeIntegration:dt_min:$KEPLERIAN_DT_MIN"
  PLOT_SNAP=""
  make_ic() {
    python3 makeIC.py -n "$KEPLERIAN_N"
    python3 -c "
import h5py
with h5py.File('initial_conditions.hdf5', 'r+') as f:
    pos = f['/PartType0/Coordinates']; x = pos[:]; x[:, 2] = 0.0; pos[...] = x
    vel = f['/PartType0/Velocities']; v = vel[:]; v[:, 2] = 0.0; vel[...] = v
    f['/Header'].attrs['Dimension'] = 2
"
  }
}
# The Keplerian ring of Cullen & Dehnen (2010), Section 4.3: Gaussian ring,
# GM = 1000, Mach 1000, 9745 equal-mass particles placed with the method of
# Cartwright et al. (2009), run with the 2D code. The test used by Cullen &
# Dehnen to build their viscosity switch (their Fig. 8 = particle positions).
test_keplerian_cd10() {
  DIM=2; EXAMPLE=HydroTests/KeplerianRing_CullenDehnen; YML=keplerian_ring.yml; IC=initial_conditions.hdf5
  FLAGS="--hydro --external-gravity"
  PARAMS="-P TimeIntegration:time_end:$KEPLERIAN_CD10_T_END -P Snapshots:delta_time:1 -P TimeIntegration:dt_min:$KEPLERIAN_DT_MIN"
  PLOT_SNAP=""
  make_ic() { python3 makeIC.py -n "$KEPLERIAN_CD10_N"; }
}

test_zeldovich() {
  DIM=3; EXAMPLE=Cosmology/ZeldovichPancake_3D; YML=zeldovichPancake.yml; IC=zeldovichPancake.hdf5
  FLAGS="--hydro --self-gravity --cosmology"
  PARAMS="-P Cosmology:a_end:$(python3 -c "print(1./(1.+$ZELDOVICH_Z_END))")"
  PLOT_SNAP=""
  make_ic() { sed -i "s/^numPart_1D = .*/numPart_1D = $ZELDOVICH_N  # Number of particles along each dimension/" makeIC.py; python3 makeIC.py; }
}

# Same problem on a glass (makeIC_glass.py): no lattice symmetry, hence no
# column regime, at the price of particle-scale density noise. The cold gas is
# Jeans-unstable at all resolved scales, so the run starts at ZELDOVICH_Z_START
# (default z = 10) rather than z = 100 to limit the linear growth of that noise
# (factor 6 instead of 50); the thermal state is the adiabatic one of the z =
# 100 start, so the analytic solution is unchanged.
test_zeldovich_glass() {
  test_zeldovich
  local a_begin
  a_begin=$(python3 -c "print(1./(1.+$ZELDOVICH_Z_START))")
  PARAMS="$PARAMS -P Cosmology:a_begin:$a_begin -P Snapshots:scale_factor_first:$a_begin"
  make_ic() {
    get_glass "$ZELDOVICH_GLASS"
    python3 makeIC_glass.py --glass "$ZELDOVICH_GLASS.hdf5" --z_start "$ZELDOVICH_Z_START" --z_end "$ZELDOVICH_Z_END"
  }
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

# Blob test (Agertz et al. 2007 set-up as in the SPHENIX paper): a dense
# (rho = 10) sphere of radius 0.1 in a Mach 2.7 wind, 4 x 1 x 1 box.
test_blob() {
  DIM=3; EXAMPLE=HydroTests/BlobTest_3D; YML=blob.yml; IC=blob.hdf5
  FLAGS="--hydro"
  PARAMS="-P TimeIntegration:time_end:$BLOB_T_END -P Snapshots:delta_time:0.348"
  PLOT_SNAP=""
  make_ic() {
    sed -i "s/num_on_side=64/num_on_side=$BLOB_N/" makeIC.py
    python3 makeIC.py
    # Recent swiftsimio writers use the plural dataset name; SWIFT reads the singular one
    python3 -c "
import h5py
with h5py.File('blob.hdf5', 'r+') as f:
    g = f['PartType0']
    if 'SmoothingLength' not in g and 'SmoothingLengths' in g:
        g.move('SmoothingLengths', 'SmoothingLength')
"
  }
}
# Gas in hydrostatic equilibrium in a fixed NFW potential (plus self-gravity);
# needs a build with the NFW external potential.
test_nfw() {
  DIM=3; BUILD=3dnfw; EXAMPLE=HydroTests/NFW_Hydrostatic; YML=NFW_Hydrostatic.yml; IC=nfw.hdf5
  FLAGS="--hydro --external-gravity --self-gravity"
  # Current SWIFT requires Gravity:epsilon_fmm > 0 even with the geometric MAC of
  # the example's file; use the standard adaptive MAC instead.
  PARAMS="-P TimeIntegration:time_end:$NFW_T_END -P Gravity:MAC:adaptive -P Gravity:epsilon_fmm:0.001"
  PLOT_SNAP=""
  make_ic() { [[ -e nfw.hdf5 ]] || wget -q $GLASS_URL/NFW_Hydrostatic/nfw.hdf5; }
}

# Taylor-Green vortex (Phantom paper, section 5.4.1) at Mach 0.1 on a 2D
# glass: a stationary Euler solution, so the kinetic-energy decay measures the
# effective numerical viscosity of the scheme in a smooth subsonic flow.
test_taylorgreen() {
  DIM=2; EXAMPLE=HydroTests/TaylorGreenVortex_2D; YML=taylorGreen.yml; IC=taylorGreen.hdf5
  FLAGS="--hydro"; PARAMS="-P TimeIntegration:time_end:$TG_T_END"; PLOT_SNAP=""
  make_ic() { get_glass glassPlane_128; python3 makeIC.py; }
}
# Shearing shock (Cullen & Dehnen 2010, section 5.3): Mach-1 colliding streams
# with a sinusoidal shear of SHEARING_S times the stream speed superposed
# (paper: 100; default 10, see the README of the example). Tests that
# a shear-aware viscosity switch still captures the shock.
test_shearing() {
  DIM=3; EXAMPLE=HydroTests/ShearingShock_3D; YML=shearingShock.yml; IC=shearingShock.hdf5
  FLAGS="--hydro"; PARAMS="-P TimeIntegration:time_end:$SHEARING_T_END -P Snapshots:delta_time:$(python3 -c "print($SHEARING_T_END/2)")"; PLOT_SNAP=""
  make_ic() { get_glass glassCube_32; python3 makeIC.py -n "$SHEARING_N" -s "$SHEARING_S"; }
}
# Strong blast wave (Rosswog 2020, section 3.4.3): planar Riemann problem with
# a pressure jump of 1e4 at uniform density (shock Mach ~ 100).
test_blast() {
  DIM=3; EXAMPLE=HydroTests/StrongBlastWave_3D; YML=blastWave.yml; IC=blastWave.hdf5
  FLAGS="--hydro --limiter"; PARAMS="-P TimeIntegration:time_end:$BLAST_T_END -P Snapshots:delta_time:$(python3 -c "print($BLAST_T_END/2)")"; PLOT_SNAP=""
  make_ic() { get_glass glassCube_32; python3 makeIC.py -n "$BLAST_N"; }
}
# Rayleigh-Taylor instability (2D, Abel 2011 / Saitoh & Makino 2013 set-up of
# the example, gamma = 7/5, constant gravity, fixed boundary layers): needs the
# 2drt build (constant potential, boundary-particle forcing, gamma = 7/5).
test_rt() {
  DIM=2; BUILD=2drt; EXAMPLE=HydroTests/Rayleigh-Taylor_2D; YML=rayleigh_taylor.yml; IC=rayleigh_taylor.hdf5
  FLAGS="--hydro --external-gravity"
  PARAMS="-P TimeIntegration:time_end:$RT_T_END -P Snapshots:delta_time:0.1"
  # The example prints the highest boundary-particle ID when the ICs are made
  [[ -e $OUT_DIR/ics/rt/boundary_id.txt ]] && PARAMS="$PARAMS -P BoundaryParticles:boundary_particle_max_id:$(cat "$OUT_DIR/ics/rt/boundary_id.txt")"
  PLOT_SNAP=""
  make_ic() {
    sed -i "s/^N = \[128, 192\]/N = [$RT_N, $((RT_N * 3 / 2))]/" makeIC.py
    python3 makeIC.py | tee makeIC.out
    grep -o "boundary_particle_max_id: [0-9]*" makeIC.out | awk '{print $2}' > boundary_id.txt
    [[ -s boundary_id.txt ]] || { echo "No boundary-particle ID reported by makeIC.py"; return 1; }
  }
}
# Kelvin-Helmholtz growth rate (McNally et al. 2012 style seed on a uniform
# density shear layer with linear ramps): the growth of the seeded mode is
# compared with the analytic rate of the piecewise-linear shear layer.
test_khgr() {
  DIM=2; EXAMPLE=HydroTests/KelvinHelmholtzGrowthRate_2D; YML=kelvinHelmholtzGrowthRate.yml; IC=kelvinHelmholtzGrowthRate.hdf5
  FLAGS="--hydro"; PARAMS="-P TimeIntegration:time_end:$KHGR_T_END"; PLOT_SNAP=""
  make_ic() { get_glass glassPlane_128; sed -i "s/^B = .*/B = $KHGR_SEED  # Amplitude of the seed velocity/" makeIC.py; python3 makeIC.py; }
}
# Interacting blast waves (Woodward & Colella 1984), 1D, gamma = 7/5: needs the
# 1d build. The reference is a high-resolution grid solution.
test_ibw() {
  DIM=1; BUILD=1d; EXAMPLE=HydroTests/InteractingBlastWaves_1D; YML=interactingBlastWaves.yml; IC=interactingBlastWaves.hdf5
  FLAGS="--hydro --limiter"; PARAMS=""; PLOT_SNAP=4
  make_ic() {
    python3 makeIC.py
    [[ -e interactingBlastWaves1D_exact.txt ]] || wget -q $REF_URL/interactingBlastWaves1D_exact.txt
  }
}

test_sod_weak() {
  # Sod shock with a weaker shock (right pressure 0.5): the contact discontinuity dominates
  DIM=3; EXAMPLE=HydroTests/SodShock_3D; YML=sodShock.yml; IC=sodShock.hdf5
  FLAGS="--hydro"; PARAMS=""; PLOT_SNAP=""
  make_ic() { get_glass glassCube_64 glassCube_32; python3 makeIC.py --P-right $SOD_WEAK_P_RIGHT; }
}
test_sod_strong() {
  # Sod shock with a stronger shock (right pressure 0.01)
  DIM=3; EXAMPLE=HydroTests/SodShock_3D; YML=sodShock.yml; IC=sodShock.hdf5
  FLAGS="--hydro"; PARAMS=""; PLOT_SNAP=""
  make_ic() { get_glass glassCube_64 glassCube_32; python3 makeIC.py --P-right $SOD_STRONG_P_RIGHT; }
}
test_kh8() {
  # Kelvin-Helmholtz at a high density contrast (sharp interface, equal-mass particles)
  DIM=2; EXAMPLE=HydroTests/KelvinHelmholtz_2D; YML=kelvinHelmholtz.yml; IC=kelvinHelmholtz.hdf5
  FLAGS="--hydro"; PARAMS="-P TimeIntegration:time_end:$KH8_T_END -P Snapshots:delta_time:0.5"; PLOT_SNAP=""
  make_ic() { sed -i "s/^L2 = .*/L2 = $KH8_L2/; s/^rho1 = .*/rho1 = $KH8_CONTRAST  # Central density/" makeIC.py; python3 makeIC.py; }
}
test_soundwave() {
  # Travelling linear sound wave: phase speed (dispersion), damping, transverse noise
  DIM=3; EXAMPLE=HydroTests/SoundWave_3D; YML=soundWave.yml; IC=soundWave.hdf5
  FLAGS="--hydro"; PARAMS="-P TimeIntegration:time_end:$SOUNDWAVE_T_END"; PLOT_SNAP=""
  make_ic() { get_glass $SOUNDWAVE_GLASS; python3 makeIC.py --glass $SOUNDWAVE_GLASS.hdf5 --amplitude $SOUNDWAVE_AMPLITUDE --nwave $SOUNDWAVE_NWAVE; }
}
test_wavedamping() {
  # Standing sound wave over many periods: damping rate of the mode energy, frequency
  DIM=3; EXAMPLE=HydroTests/SoundWave_3D; YML=soundWave.yml; IC=soundWave.hdf5
  FLAGS="--hydro"; PARAMS="-P TimeIntegration:time_end:$WAVEDAMPING_T_END"; PLOT_SNAP=""
  make_ic() { get_glass $SOUNDWAVE_GLASS; python3 makeIC.py --glass $SOUNDWAVE_GLASS.hdf5 --amplitude $SOUNDWAVE_AMPLITUDE --nwave $SOUNDWAVE_NWAVE --standing; }
}
test_thermal() {
  # Thermal discontinuity in pressure equilibrium: artificial conduction at a contact, pressure blip
  DIM=2; EXAMPLE=HydroTests/ThermalDiffusion_2D; YML=thermalDiffusion.yml; IC=thermalDiffusion.hdf5
  FLAGS="--hydro"; PARAMS="-P TimeIntegration:time_end:$THERMAL_T_END"; PLOT_SNAP=""
  make_ic() { get_glass glassPlane_128; python3 makeIC.py --contrast $THERMAL_CONTRAST; }
}
test_gresho_boost() {
  # Gresho vortex advected at ~Mach 1 through the box: Galilean invariance (compare with gresho)
  DIM=2; EXAMPLE=HydroTests/GreshoVortex_2D; YML=gresho.yml; IC=greshoVortex.hdf5
  FLAGS="--hydro"; PARAMS=""; PLOT_SNAP=""
  make_ic() { get_glass glassPlane_128; python3 makeIC.py --boost $GRESHO_BOOST 0; }
}
get_glass() {
  for g in "$@"; do
    [[ -e $g.hdf5 ]] || wget -q "$GLASS_URL/$g.hdf5" || { echo "Cannot fetch $g"; return 1; }
  done
}

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
# Build flavours: suffix -> configure flags added to COMMON_FLAGS (the
# external potential of COMMON_FLAGS is replaced for the NFW build).
declare -A BUILD_FLAGS=(
  [1d]="--with-hydro-dimension=1 --with-adiabatic-index=7/5 --with-kernel=$KERNEL_1D"
  [2d]="--with-hydro-dimension=2"
  [2drt]="--with-hydro-dimension=2 --with-ext-potential=constant --with-forcing=boundary-particles --with-adiabatic-index=7/5"
  [3d]=""
  [3dnfw]="--with-ext-potential=nfw"
)

build_one() {
  local name=$1 build=$2
  local link=$OUT_DIR/builds/${name}_${build}
  local common=$COMMON_FLAGS
  [[ ${BUILD_FLAGS[$build]} == *--with-ext-potential* ]] && common=${common/--with-ext-potential=point-mass-softened/}
  [[ ${BUILD_FLAGS[$build]} == *--with-kernel* ]] && common=${common/--with-kernel=$KERNEL/}
  # One build per distinct set of configure flags, shared through a link
  local key
  key=$(echo "$common ${BUILD_FLAGS[$build]} ${SCHEME_FLAGS[$name]}" | md5sum | cut -c1-10)
  local dir=$OUT_DIR/builds/by_flags/${key}_${build}
  mkdir -p "$OUT_DIR/builds/by_flags"
  if [[ ! -L $link || $(readlink "$link") != "by_flags/${key}_${build}" ]]; then rm -rf "$link"; ln -s "by_flags/${key}_${build}" "$link"; fi
  if [[ $DO_BUILD -eq 0 && -x $dir/swift ]]; then return 0; fi
  log "Building $name ($build) in $dir"
  mkdir -p "$dir"
  # Copy the sources only: no build products of the source tree (in particular
  # not its swift binary, which would mask a failed build here).
  rsync -a --delete --exclude='.git' --exclude='*.hdf5' --exclude='*.o' \
    --exclude='*.lo' --exclude='*.la' --exclude='*.a' --exclude='.libs' \
    --exclude='/swift' --exclude='/swift_mpi' --exclude='/swift_fof' \
    --exclude='/swift_fof_mpi' --exclude='*.png' --exclude='*.pdf' \
    --exclude='restart' --exclude='suite_output' --exclude='results_*' \
    "$SRC_DIR/" "$dir/" || return 1
  if ! ( cd "$dir" && ./autogen.sh && ./configure $common ${BUILD_FLAGS[$build]} ${SCHEME_FLAGS[$name]} \
      && make -j"$THREADS" ) > "$dir/build.log" 2>&1; then
    echo "Build of ${name}_${build} failed, see $dir/build.log"; return 1
  fi
  [[ -x $dir/swift ]] || { echo "Build of ${name}_${build} produced no binary"; return 1; }
}

# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------
prepare_ic() {
  local test=$1
  test_$test
  local dir=$OUT_DIR/ics/$test
  if [[ -e $dir/$IC ]]; then return 0; fi
  mkdir -p "$dir"
  # Concurrent suite instances (parameter studies) may ask for the same ICs
  # at the same time: let one generate them while the others wait.
  if command -v flock > /dev/null; then
    (
      flock 9 || exit 1
      [[ -e $dir/$IC ]] && exit 0
      prepare_ic_locked "$test" "$dir"
    ) 9> "$dir/.lock"
    return $?
  fi
  prepare_ic_locked "$test" "$dir"
}

prepare_ic_locked() {
  local test=$1 dir=$2
  log "Generating initial conditions for $test"
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
  BUILD=""
  test_$test
  [[ -z $BUILD ]] && BUILD=${DIM}d
  local icdir=$OUT_DIR/ics/$test
  local dir=$OUT_DIR/runs/$scheme/$test
  local swift=$OUT_DIR/builds/${scheme}_${BUILD}/swift
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
if [[ $ICS_ONLY -eq 1 ]]; then
  for t in $TESTS; do prepare_ic "$t" || exit 1; done
  log "Initial conditions ready in $OUT_DIR/ics"
  exit 0
fi

# Only the builds needed by the selected tests
NEEDED_BUILDS=""
for t in $TESTS; do
  BUILD=""; test_$t; [[ -z $BUILD ]] && BUILD=${DIM}d
  [[ " $NEEDED_BUILDS " == *" $BUILD "* ]] || NEEDED_BUILDS="$NEEDED_BUILDS $BUILD"
done
for s in "${SCHEME_NAMES[@]}"; do
  for b in $NEEDED_BUILDS; do build_one "$s" "$b" || exit 1; done
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
