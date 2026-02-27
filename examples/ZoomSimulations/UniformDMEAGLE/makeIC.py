"""Generate DM-only ICs for the ``UniformDMEAGLE`` zoom example.

This script constructs a deterministic two-component DM-only IC file:

- ``PartType1``: high-resolution particles in a central zoom cube.
- ``PartType2``: low-resolution background particles outside that cube.

Geometry and hierarchy:
- background grid with ``16`` top-level cells per axis;
- zoom region spanning the central ``2x2x2`` background cells
  (12.5% of the box width per axis);
- zoom top-level depth fixed to ``2`` (``8`` zoom cells per axis);
- ``region_pad_factor = 1.1``.

Mass model:
- cosmology follows EAGLE reference values;
- high-resolution parent masses are chosen so that, after
  ``generate_gas_in_ics``, the remaining DM mass is ~``9.7e6 Msun``;
- background particles are ~``10x`` more massive than high-resolution parent
  particles (i.e. ``10x`` lower mass resolution), while preserving mean matter
  density in the background volume.

Dependencies:
- ``numpy`` for geometry and sampling;
- ``h5py`` for direct HDF5 edits;
- ``unyt`` + ``swiftsimio`` for SWIFT-compliant unit-aware IC writing.
"""

import sys
from typing import Any, cast

try:
    import h5py
except ImportError:
    print("ERROR: Missing dependency 'h5py'.")
    print("Install it with: python3 -m pip install h5py")
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("ERROR: Missing dependency 'numpy'.")
    print("Install it with: python3 -m pip install numpy")
    sys.exit(1)

try:
    from unyt import Mpc, Msun, km, s, unyt_array
except ImportError:
    print("ERROR: Missing dependency 'unyt'.")
    print("Install it with: python3 -m pip install unyt")
    sys.exit(1)

try:
    from swiftsimio import Writer
    from swiftsimio.units import cosmo_units
except ImportError:
    print("ERROR: Missing dependency 'swiftsimio'.")
    print("Install it with: python3 -m pip install swiftsimio")
    sys.exit(1)


# Cosmology and target mass scale.
H_PARAM = 0.6777
OMEGA_CDM = 0.2587481
OMEGA_B = 0.0482519
OMEGA_M = OMEGA_CDM + OMEGA_B
G_NEWTON_INTERNAL = 43.01069
RHO_M = OMEGA_M * 3.0 * (100.0 * H_PARAM) ** 2 / (8.0 * np.pi * G_NEWTON_INTERNAL)

TARGET_DM_POST_SPLIT_MSUN = 9.7e6
TARGET_DM_POST_SPLIT = TARGET_DM_POST_SPLIT_MSUN / 1.0e10
FBARYON = OMEGA_B / OMEGA_M
TARGET_PARENT_MASS = TARGET_DM_POST_SPLIT / (1.0 - FBARYON)

# Zoom geometry.
REGION_PAD_FACTOR = 1.1
BKG_CELLS = 16
ZOOM_BKG_CELLS_SPANNED = 2
ZOOM_TOP_LEVEL_DEPTH = 2
ZOOM_CELLS_PER_AXIS = ZOOM_BKG_CELLS_SPANNED * (2**ZOOM_TOP_LEVEL_DEPTH)

# Particle layout.
HI_PARTICLES_AXIS = 1
HI_PARTICLES_PER_ZOOM_CELL = HI_PARTICLES_AXIS**3
BACKGROUND_MASS_RATIO_TARGET = 10.0
RNG_SEED = 20260227

OUTPUT_FILE = "zoom_uniform_dm_eagle.hdf5"
A_BEGIN = 0.1


def solve_box_size() -> float:
    """Solve the box size required to hit the target high-res mass.

    Returns:
        Box size in internal length units (Mpc for this setup).
    """
    n_hi = ZOOM_CELLS_PER_AXIS**3 * HI_PARTICLES_PER_ZOOM_CELL
    hi_fraction = (ZOOM_BKG_CELLS_SPANNED / BKG_CELLS / REGION_PAD_FACTOR) ** 3
    return float((TARGET_PARENT_MASS * n_hi / (RHO_M * hi_fraction)) ** (1.0 / 3.0))


BOX_SIZE = solve_box_size()
TARGET_ZOOM_EXTENT = ZOOM_BKG_CELLS_SPANNED / BKG_CELLS * BOX_SIZE
HI_REGION_SIZE = TARGET_ZOOM_EXTENT / REGION_PAD_FACTOR


def occupancy_3d(positions: np.ndarray, ncell: int, lo: float, hi: float) -> np.ndarray:
    """Return per-cell occupancy in an ``ncell^3`` Cartesian grid.

    Args:
        positions: ``(N, 3)`` particle coordinates.
        ncell: Number of cells per axis.
        lo: Lower coordinate bound of the grid.
        hi: Upper coordinate bound of the grid.

    Returns:
        Flat array of length ``ncell^3`` with integer counts per cell.
    """
    width = hi - lo
    frac = (positions - lo) / width
    idx = np.floor(frac * ncell).astype(np.int64)
    idx = np.clip(idx, 0, ncell - 1)
    linear = idx[:, 0] + ncell * (idx[:, 1] + ncell * idx[:, 2])
    return np.bincount(linear, minlength=ncell**3)


def sample_stratified_in_cell(
    rng: np.random.Generator, lo: np.ndarray, hi: np.ndarray, n_axis: int
) -> np.ndarray:
    """Sample ``n_axis^3`` points with one jittered point per stratum.

    Args:
        rng: NumPy random generator.
        lo: Lower coordinate bounds of the cell.
        hi: Upper coordinate bounds of the cell.
        n_axis: Number of strata per axis.

    Returns:
        ``(n_axis^3, 3)`` array of positions in ``[lo, hi)``.
    """
    width = hi - lo
    offset = (np.arange(n_axis, dtype=np.float64) + 0.5) / n_axis
    gx, gy, gz = np.meshgrid(offset, offset, offset, indexing="ij")
    centers = np.column_stack((gx.ravel(), gy.ravel(), gz.ravel()))
    jitter = rng.uniform(-0.45, 0.45, size=centers.shape) / n_axis
    u = np.clip(centers + jitter, 0.0, 1.0 - np.finfo(np.float64).eps)
    return lo + u * width


def sample_uniform_in_cell(
    rng: np.random.Generator, lo: np.ndarray, hi: np.ndarray, n_part: int
) -> np.ndarray:
    """Sample ``n_part`` positions uniformly in one Cartesian cell.

    Args:
        rng: NumPy random generator.
        lo: Lower coordinate bounds of the cell.
        hi: Upper coordinate bounds of the cell.
        n_part: Number of particles to sample.

    Returns:
        ``(n_part, 3)`` array of sampled positions.
    """
    return rng.uniform(lo, hi, size=(n_part, 3))


def generate_high_res_particles(rng: np.random.Generator, centre: float) -> np.ndarray:
    """Populate all zoom cells with fixed high-resolution occupancy.

    Args:
        rng: NumPy random generator.
        centre: Zoom region centre coordinate.

    Returns:
        ``(N, 3)`` array of high-resolution ``PartType1`` positions.
    """
    half = 0.5 * HI_REGION_SIZE
    hi_lo = centre - half
    zoom_w = HI_REGION_SIZE / ZOOM_CELLS_PER_AXIS

    chunks = []
    for ix in range(ZOOM_CELLS_PER_AXIS):
        for iy in range(ZOOM_CELLS_PER_AXIS):
            for iz in range(ZOOM_CELLS_PER_AXIS):
                lo = hi_lo + zoom_w * np.array([ix, iy, iz], dtype=np.float64)
                hi = lo + zoom_w
                chunks.append(sample_stratified_in_cell(rng, lo, hi, HI_PARTICLES_AXIS))
    return np.vstack(chunks)


def background_particles_per_top_cell() -> int:
    """Particles per background top cell for target lower mass resolution.

    We enforce a target background-to-highres particle mass ratio and solve for
    the total number of background particles required by mass conservation.

    Returns:
        Integer number of ``PartType2`` particles per non-central top-level
        background cell.
    """
    total_background_mass = RHO_M * (BOX_SIZE**3 - HI_REGION_SIZE**3)
    target_background_mass = BACKGROUND_MASS_RATIO_TARGET * TARGET_PARENT_MASS

    n_background_target = max(
        1, int(round(total_background_mass / target_background_mass))
    )
    n_background_cells = BKG_CELLS**3 - ZOOM_BKG_CELLS_SPANNED**3
    return max(1, int(round(n_background_target / n_background_cells)))


def generate_background_particles(rng: np.random.Generator) -> np.ndarray:
    """Populate all non-central top-level background cells.

    Args:
        rng: NumPy random generator.

    Returns:
        ``(N, 3)`` array of background ``PartType2`` positions.
    """
    top_w = BOX_SIZE / BKG_CELLS
    central = {BKG_CELLS // 2 - 1, BKG_CELLS // 2}
    n_per_cell = background_particles_per_top_cell()

    chunks = []
    for ix in range(BKG_CELLS):
        for iy in range(BKG_CELLS):
            for iz in range(BKG_CELLS):
                if ix in central and iy in central and iz in central:
                    continue
                lo = top_w * np.array([ix, iy, iz], dtype=np.float64)
                hi = lo + top_w
                chunks.append(sample_uniform_in_cell(rng, lo, hi, n_per_cell))
    return np.vstack(chunks)


def validate_layout(hi_pos: np.ndarray, bkg_pos: np.ndarray, centre: float) -> None:
    """Validate occupancy and zoom/background spatial separation.

    Args:
        hi_pos: High-resolution ``PartType1`` positions.
        bkg_pos: Background ``PartType2`` positions.
        centre: Zoom region centre coordinate.
    """
    all_pos = np.vstack((hi_pos, bkg_pos))
    bkg_occ = occupancy_3d(all_pos, BKG_CELLS, lo=0.0, hi=BOX_SIZE)
    if np.any(bkg_occ == 0):
        raise RuntimeError("Found empty top-level background cells.")

    half = 0.5 * HI_REGION_SIZE
    lo = centre - half
    hi = centre + half
    zoom_occ = occupancy_3d(hi_pos, ZOOM_CELLS_PER_AXIS, lo=lo, hi=hi)
    if np.any(zoom_occ != HI_PARTICLES_PER_ZOOM_CELL):
        raise RuntimeError("High-resolution zoom-cell occupancy is not uniform.")

    if np.any((hi_pos < lo) | (hi_pos >= hi)):
        raise RuntimeError("High-resolution particles found outside zoom cube.")

    bkg_in_hi = np.all((bkg_pos >= lo) & (bkg_pos < hi), axis=1)
    if np.any(bkg_in_hi):
        raise RuntimeError("Background particles found inside zoom cube.")


def write_ic_file(hi_pos: np.ndarray, bkg_pos: np.ndarray) -> None:
    """Write SWIFT-compatible DM-only ICs with PartType1 + PartType2.

    Args:
        hi_pos: High-resolution ``PartType1`` positions.
        bkg_pos: Background ``PartType2`` positions.
    """
    n_hi = hi_pos.shape[0]
    n_bkg = bkg_pos.shape[0]

    hi_mass = np.full(n_hi, TARGET_PARENT_MASS, dtype=np.float64)
    bkg_mass_total = RHO_M * (BOX_SIZE**3 - HI_REGION_SIZE**3)
    bkg_mass = np.full(n_bkg, bkg_mass_total / n_bkg, dtype=np.float64)

    hi_vel = np.zeros_like(hi_pos)
    bkg_vel = np.zeros_like(bkg_pos)

    # Write the high-resolution DM component with swiftsimio.
    box_size = unyt_array([BOX_SIZE, BOX_SIZE, BOX_SIZE], Mpc)
    writer = Writer(cosmo_units, cast(Any, box_size), dimension=3)
    writer_any: Any = writer
    writer_any.dark_matter.coordinates = cast(Any, hi_pos) * Mpc
    writer_any.dark_matter.velocities = cast(Any, hi_vel) * km / s
    writer_any.dark_matter.masses = cast(Any, hi_mass) * 1e10 * Msun
    writer.write(OUTPUT_FILE)

    # Append the background DM component and update header metadata.
    with h5py.File(OUTPUT_FILE, "r+") as handle:
        part2 = handle.create_group("PartType2")
        part2.create_dataset("Coordinates", data=bkg_pos, compression="gzip")
        part2.create_dataset("Velocities", data=bkg_vel, compression="gzip")
        part2.create_dataset("Masses", data=bkg_mass, compression="gzip")
        part2.create_dataset(
            "ParticleIDs",
            data=np.arange(n_hi + 1, n_hi + n_bkg + 1, dtype=np.int64),
            compression="gzip",
        )

        header = handle["Header"].attrs
        num_this: Any = header["NumPart_ThisFile"]
        num_tot: Any = header["NumPart_Total"]
        mass_table: Any = header["MassTable"]

        num_this[2] = n_bkg
        num_tot[2] = n_bkg
        mass_table[2] = float(bkg_mass[0])

        header["NumPart_ThisFile"] = num_this
        header["NumPart_Total"] = num_tot
        header["MassTable"] = mass_table
        header["Time"] = A_BEGIN

        if "RuntimePars" not in handle:
            runtime = handle.create_group("RuntimePars")
            runtime.attrs["PeriodicBoundariesOn"] = 1

    dm_post_split = TARGET_PARENT_MASS * (1.0 - FBARYON) * 1e10
    gas_post_split = TARGET_PARENT_MASS * FBARYON * 1e10
    print(f"Wrote {OUTPUT_FILE}")
    print(f"Box size [Mpc]:                           {BOX_SIZE:.6f}")
    print(f"High-resolution particles (PartType1):    {n_hi}")
    print(f"Background particles (PartType2):         {n_bkg}")
    print(
        f"Background particles/top-level cell:      {background_particles_per_top_cell()}"
    )
    print(f"Target post-split DM mass [Msun]:         {dm_post_split:.3e}")
    print(f"Target post-split gas mass [Msun]:        {gas_post_split:.3e}")
    print(
        f"Zoom-region width fraction per axis:      {TARGET_ZOOM_EXTENT / BOX_SIZE:.3f}"
    )
    print(
        f"Background/high-res mass ratio target:    {BACKGROUND_MASS_RATIO_TARGET:.1f}"
    )
    print(f"Background/high-res mass ratio actual:    {bkg_mass[0] / hi_mass[0]:.2f}")


def main() -> None:
    """Generate, validate, and write the ``UniformDMEAGLE`` IC file."""
    # Offset the centre by a tiny deterministic amount to avoid exact lattice
    # symmetries at depth-1 multipole checks while staying in the same central
    # background-cell block.
    centre = 0.5003 * BOX_SIZE
    rng = np.random.default_rng(RNG_SEED)

    # Build particle positions from the target zoom hierarchy.
    hi_pos = generate_high_res_particles(rng, centre)
    bkg_pos = generate_background_particles(rng)

    # Validate the layout before writing the output file.
    validate_layout(hi_pos, bkg_pos, centre)
    write_ic_file(hi_pos, bkg_pos)


if __name__ == "__main__":
    main()
