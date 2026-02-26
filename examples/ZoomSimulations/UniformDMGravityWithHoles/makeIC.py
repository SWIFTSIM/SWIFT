"""Generate deterministic zoom-region dark-matter ICs with empty zoom cells.

This script builds a synthetic DMO zoom simulation example with particles
distributed uniformly into the top-level cells of the background and zoom
hierarchy, while deliberately removing a deterministic set of zoom cells to
exercise empty-cell edge cases.

Geometry and hierarchy:
- periodic unit box of side 1.0;
- background grid with 8 top-level cells per axis;
- central zoom region spanning 2x2x2 background cells;
- zoom top-level depth fixed to 2, giving 8 zoom cells per axis;
- user pad factor fixed to 1.1.

Particle construction:
- PartType1 (high-resolution DM): every zoom cell receives 4^3 particles except
  a deterministic set of two 2x2x2 hole blocks (16 empty zoom cells total);
- PartType2 (background DM): every non-central background top-level cell is
  populated with exactly 4^3 = 64 particles, also stratified within-cell.

Example requires:
- ``h5py`` for direct HDF5 manipulation to append PartType2.
- ``numpy`` for array handling and random sampling.
- ``unyt`` for unit handling and conversion.
- ``swiftsimio`` for writing SWIFT compliant IC files with appropriate
  metadata.

Validation performed before writing:
- no empty background top-level cells in the combined distribution;
- empty zoom cells match the requested hole pattern exactly;
- all non-hole zoom cells retain uniform occupancy;
- no PartType2 particles inside the high-resolution cube;
- no PartType1 particles outside the high-resolution cube.

Output:
- ``zoom_uniform_dm_gravity_holes.hdf5`` with SWIFT-compliant Header/Units
  metadata, PartType1 written through ``swiftsimio.Writer``, and PartType2
  appended via ``h5py``.
"""

import sys
from typing import Any

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
    from unyt import Mpc, Msun, km, s
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


# Set up constants describing the set up, cosmology, and particle layout.
BOX_SIZE = 1.0
A_BEGIN = 0.1
H0_INTERNAL = 100.0
G_NEWTON_INTERNAL = 43.01069
OMEGA_CDM = 1.0
RHO_DM = OMEGA_CDM * 3.0 * H0_INTERNAL**2 / (8.0 * np.pi * G_NEWTON_INTERNAL)
REGION_PAD_FACTOR = 1.1
RNG_SEED = 20260227

# We target a zoom extent that spans 2x2x2 background cells for bkg cdim=8.
BKG_CELLS = 8
ZOOM_BKG_CELLS_SPANNED = 2
ZOOM_TOP_LEVEL_DEPTH = 2
TARGET_ZOOM_EXTENT = ZOOM_BKG_CELLS_SPANNED / BKG_CELLS
HI_REGION_SIZE = TARGET_ZOOM_EXTENT / REGION_PAD_FACTOR

# Populate the expected cell hierarchy directly.
ZOOM_CELLS_PER_AXIS = ZOOM_BKG_CELLS_SPANNED * (2**ZOOM_TOP_LEVEL_DEPTH)
HI_PARTICLES_AXIS = 4
BKG_PARTICLES_AXIS = 4
HI_PARTICLES_PER_ZOOM_CELL = HI_PARTICLES_AXIS**3
BKG_PARTICLES_PER_TOP_CELL = BKG_PARTICLES_AXIS**3

# Output file name for the generated ICs, this needs to match the param file.
OUTPUT_FILE = "zoom_uniform_dm_gravity_holes.hdf5"

if BKG_CELLS != 8 or ZOOM_TOP_LEVEL_DEPTH != 2 or ZOOM_CELLS_PER_AXIS != 8:
    raise RuntimeError(
        "This test is defined for bkg_top_level_cells=8, zoom_top_level_depth=2, "
        "and 8 zoom cells per axis."
    )


def occupancy_3d(positions, ncell, lo=0.0, hi=1.0):
    """Return per-cell occupancy in an ``ncell^3`` Cartesian grid.

    Args:
        positions: (N, 3) array of particle positions in the unit box.
        ncell: number of cells per axis in the grid.
        lo: lower bound of the grid (default 0.0).
        hi: upper bound of the grid (default 1.0).

    Returns:
        A flat array of length ``ncell^3`` containing particle counts per cell.
    """
    width = hi - lo
    frac = (positions - lo) / width
    idx = np.floor(frac * ncell).astype(np.int64)
    idx = np.clip(idx, 0, ncell - 1)
    linear = idx[:, 0] + ncell * (idx[:, 1] + ncell * idx[:, 2])
    return np.bincount(linear, minlength=ncell**3)


def sample_stratified_in_cell(rng, lo, hi, n_axis):
    """Sample ``n_axis^3`` points stratified within one Cartesian cell.

    Args:
        rng: NumPy random number generator.
        lo: 3-vector of lower cell bounds.
        hi: 3-vector of upper cell bounds.
        n_axis: Number of strata per axis.

    Returns:
        A ``(n_axis^3, 3)`` array of positions inside ``[lo, hi)``.
    """
    width = hi - lo

    # Build sub-cell centres in normalised [0, 1)^3 coordinates.
    offset = (np.arange(n_axis, dtype=np.float64) + 0.5) / n_axis
    gx, gy, gz = np.meshgrid(offset, offset, offset, indexing="ij")
    centers = np.column_stack((gx.ravel(), gy.ravel(), gz.ravel()))

    # Jitter each centre within its own stratum, then map to physical coords.
    jitter = rng.uniform(-0.45, 0.45, size=centers.shape) / n_axis
    u = np.clip(centers + jitter, 0.0, 1.0 - np.finfo(np.float64).eps)
    return lo + u * width


def zoom_linear_index(ix, iy, iz):
    """Map a zoom-cell ``(i, j, k)`` index to a flattened linear index."""
    n = ZOOM_CELLS_PER_AXIS
    return ix + n * (iy + n * iz)


def hole_linear_indices():
    """Return the deterministic set of zoom-cell indices to leave empty."""
    hole_cells = set()

    # Hole A: 2x2x2 block near one outskirts corner of the zoom region.
    for ix in (5, 6):
        for iy in (1, 2):
            for iz in (4, 5):
                hole_cells.add(zoom_linear_index(ix, iy, iz))

    # Hole B: another 2x2x2 block on a different outskirts side.
    for ix in (1, 2):
        for iy in (5, 6):
            for iz in (2, 3):
                hole_cells.add(zoom_linear_index(ix, iy, iz))

    return hole_cells


def generate_high_res_particles(rng, centre):
    """Generate high-resolution particles, excluding deterministic hole cells.

    Args:
        rng: NumPy random number generator.
        centre: Zoom-region centre in box coordinates.

    Returns:
        Tuple ``(positions, removed_cells, hole_ids)`` where:
        - ``positions`` is the PartType1 position array;
        - ``removed_cells`` is the number of empty zoom cells;
        - ``hole_ids`` is the set of flattened hole cell indices.
    """
    half = 0.5 * HI_REGION_SIZE
    hi_lo = centre - half
    zoom_w = HI_REGION_SIZE / ZOOM_CELLS_PER_AXIS
    hole_ids = hole_linear_indices()

    chunks = []
    removed_cells = 0
    for ix in range(ZOOM_CELLS_PER_AXIS):
        for iy in range(ZOOM_CELLS_PER_AXIS):
            for iz in range(ZOOM_CELLS_PER_AXIS):
                lid = zoom_linear_index(ix, iy, iz)
                if lid in hole_ids:
                    removed_cells += 1
                    continue

                lo = hi_lo + zoom_w * np.array([ix, iy, iz], dtype=np.float64)
                hi = lo + zoom_w
                chunks.append(sample_stratified_in_cell(rng, lo, hi, HI_PARTICLES_AXIS))

    return np.vstack(chunks), removed_cells, hole_ids


def generate_background_particles(rng):
    """Generate low-resolution background particles.

    Fills all top-level background cells except the central ``2x2x2`` block that
    is occupied by the high-resolution region.

    Args:
        rng: NumPy random number generator.

    Returns:
        A ``(N, 3)`` array of PartType2 positions.
    """
    top_w = BOX_SIZE / BKG_CELLS
    central = {BKG_CELLS // 2 - 1, BKG_CELLS // 2}

    chunks = []
    for ix in range(BKG_CELLS):
        for iy in range(BKG_CELLS):
            for iz in range(BKG_CELLS):
                if ix in central and iy in central and iz in central:
                    continue

                lo = top_w * np.array([ix, iy, iz], dtype=np.float64)
                hi = lo + top_w
                chunks.append(
                    sample_stratified_in_cell(rng, lo, hi, BKG_PARTICLES_AXIS)
                )

    return np.vstack(chunks)


def validate_layout(hi_pos, bkg_pos, centre, hole_ids):
    """Validate occupancy and geometric consistency for the holes setup.

    Args:
        hi_pos: PartType1 positions.
        bkg_pos: PartType2 positions.
        centre: Zoom-region centre in box coordinates.
        hole_ids: Set of flattened zoom-cell indices that must be empty.
    """
    all_pos = np.vstack((hi_pos, bkg_pos))

    # Every background top-level cell (8^3) should contain particles.
    bkg_occ = occupancy_3d(all_pos, BKG_CELLS)
    empty_bkg = int(np.count_nonzero(bkg_occ == 0))
    if empty_bkg != 0:
        raise RuntimeError(
            f"Found {empty_bkg} empty top-level background cells (expected 0)."
        )

    half = 0.5 * HI_REGION_SIZE
    hi_lo = centre - half
    hi_hi = centre + half
    zoom_occ = occupancy_3d(hi_pos, ZOOM_CELLS_PER_AXIS, lo=hi_lo, hi=hi_hi)

    # Ensure exactly the requested hole cells are empty.
    empty_zoom_ids = set(np.where(zoom_occ == 0)[0].tolist())
    if empty_zoom_ids != hole_ids:
        raise RuntimeError(
            "Hole-cell pattern does not match requested empty zoom cells."
        )

    # All non-hole cells should have full occupancy.
    for lid, count in enumerate(zoom_occ):
        if lid in hole_ids:
            continue
        if count != HI_PARTICLES_PER_ZOOM_CELL:
            raise RuntimeError("Non-hole zoom cell occupancy is not uniform.")

    if np.any((hi_pos < hi_lo) | (hi_pos >= hi_hi)):
        raise RuntimeError("High-res particles found outside zoom cube.")
    bkg_in_hi = np.all((bkg_pos >= hi_lo) & (bkg_pos < hi_hi), axis=1)
    if np.any(bkg_in_hi):
        raise RuntimeError("Background particles found inside zoom cube.")

    print(
        f"Top-level background cell occupancy min/max: {bkg_occ.min()}/{bkg_occ.max()}"
    )
    print(
        f"High-res zoom sub-cell occupancy min/max:    {zoom_occ.min()}/{zoom_occ.max()}"
    )
    print(f"Number of empty zoom cells from holes:       {len(hole_ids)}")


def write_ic_file(hi_pos, bkg_pos):
    """Write SWIFT-compatible ICs containing PartType1 and PartType2.

    Args:
        hi_pos: PartType1 positions.
        bkg_pos: PartType2 positions.
    """
    hi_count = hi_pos.shape[0]
    bkg_count = bkg_pos.shape[0]

    hi_mass = np.full(hi_count, RHO_DM * (HI_REGION_SIZE**3) / hi_count)
    bkg_mass_total = RHO_DM * (BOX_SIZE**3 - HI_REGION_SIZE**3)
    bkg_mass = np.full(bkg_count, bkg_mass_total / bkg_count)

    hi_vel = np.zeros_like(hi_pos)
    bkg_vel = np.zeros_like(bkg_pos)

    # First write the high-resolution component through swiftsimio.
    writer = Writer(
        cosmo_units,
        np.array([BOX_SIZE, BOX_SIZE, BOX_SIZE]) * Mpc,
        dimension=3,
    )
    writer_any: Any = writer
    writer_any.dark_matter.coordinates = hi_pos * Mpc
    writer_any.dark_matter.velocities = hi_vel * km / s
    writer_any.dark_matter.masses = hi_mass * 1e10 * Msun
    writer.write(OUTPUT_FILE)

    # Then append the background component and finalise header metadata.
    with h5py.File(OUTPUT_FILE, "r+") as handle:
        part2 = handle.create_group("PartType2")
        part2.create_dataset("Coordinates", data=bkg_pos, compression="gzip")
        part2.create_dataset("Velocities", data=bkg_vel, compression="gzip")
        part2.create_dataset("Masses", data=bkg_mass, compression="gzip")
        part2.create_dataset(
            "ParticleIDs",
            data=np.arange(hi_count + 1, hi_count + bkg_count + 1, dtype=np.int64),
            compression="gzip",
        )

        header = handle["Header"].attrs
        num_this: Any = header["NumPart_ThisFile"]
        num_tot: Any = header["NumPart_Total"]
        mass_table: Any = header["MassTable"]

        num_this[2] = bkg_count
        num_tot[2] = bkg_count
        mass_table[2] = float(bkg_mass[0])

        header["NumPart_ThisFile"] = num_this
        header["NumPart_Total"] = num_tot
        header["MassTable"] = mass_table
        header["Time"] = A_BEGIN

        if "RuntimePars" not in handle:
            runtime = handle.create_group("RuntimePars")
            runtime.attrs["PeriodicBoundariesOn"] = 1

    print(f"Wrote {OUTPUT_FILE}")
    print(f"High-resolution particles (PartType1): {hi_count}")
    print(f"Background particles (PartType2):      {bkg_count}")
    print(f"High-res region side-length:           {HI_REGION_SIZE:.6f}")
    print(f"Target zoom extent after pad:          {TARGET_ZOOM_EXTENT:.6f}")


def main():
    """Generate, validate, and write the holes IC file."""
    centre = 0.5
    rng = np.random.default_rng(RNG_SEED)

    # Build particle positions from the target cell hierarchy.
    hi_pos, removed_cells, hole_ids = generate_high_res_particles(rng, centre)
    bkg_pos = generate_background_particles(rng)

    # Validate geometry/occupancy assumptions and write output.
    validate_layout(hi_pos, bkg_pos, centre, hole_ids)
    print(
        f"Removed high-resolution particles in holes: {removed_cells * HI_PARTICLES_PER_ZOOM_CELL}"
    )
    write_ic_file(hi_pos, bkg_pos)


if __name__ == "__main__":
    main()
