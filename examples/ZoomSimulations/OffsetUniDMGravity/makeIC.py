"""Generate deterministic offset zoom-region dark-matter ICs.

This script mirrors ``UniformDMGravity`` but shifts the high-resolution region
centre by +1/4 of the box along each axis. The offset region stays inside the
box (i.e. does not cross periodic boundaries in the ICs).
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


BOX_SIZE = 1.0
A_BEGIN = 0.1
H0_INTERNAL = 100.0
G_NEWTON_INTERNAL = 43.01069
OMEGA_CDM = 1.0
RHO_DM = OMEGA_CDM * 3.0 * H0_INTERNAL**2 / (8.0 * np.pi * G_NEWTON_INTERNAL)
REGION_PAD_FACTOR = 1.1
RNG_SEED = 20260228

BKG_CELLS = 8
ZOOM_BKG_CELLS_SPANNED = 2
ZOOM_TOP_LEVEL_DEPTH = 2
TARGET_ZOOM_EXTENT = ZOOM_BKG_CELLS_SPANNED / BKG_CELLS
HI_REGION_SIZE = TARGET_ZOOM_EXTENT / REGION_PAD_FACTOR

ZOOM_CELLS_PER_AXIS = ZOOM_BKG_CELLS_SPANNED * (2**ZOOM_TOP_LEVEL_DEPTH)
HI_PARTICLES_AXIS = 4
BKG_PARTICLES_AXIS = 4
HI_PARTICLES_PER_ZOOM_CELL = HI_PARTICLES_AXIS**3

OUTPUT_FILE = "offset_uni_dm_gravity.hdf5"


def occupancy_3d(positions, ncell, lo=0.0, hi=1.0):
    """Return per-cell occupancy in an ``ncell^3`` Cartesian grid."""
    width = hi - lo
    frac = (positions - lo) / width
    idx = np.floor(frac * ncell).astype(np.int64)
    idx = np.clip(idx, 0, ncell - 1)
    linear = idx[:, 0] + ncell * (idx[:, 1] + ncell * idx[:, 2])
    return np.bincount(linear, minlength=ncell**3)


def sample_stratified_in_cell(rng, lo, hi, n_axis):
    """Sample ``n_axis^3`` points stratified within one Cartesian cell."""
    width = hi - lo
    offset = (np.arange(n_axis, dtype=np.float64) + 0.5) / n_axis
    gx, gy, gz = np.meshgrid(offset, offset, offset, indexing="ij")
    centers = np.column_stack((gx.ravel(), gy.ravel(), gz.ravel()))
    jitter = rng.uniform(-0.45, 0.45, size=centers.shape) / n_axis
    u = np.clip(centers + jitter, 0.0, 1.0 - np.finfo(np.float64).eps)
    return lo + u * width


def generate_high_res_particles(rng, centre):
    """Generate high-resolution particles in all zoom cells."""
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


def zoom_top_cell_indices(centre):
    """Return top-level background cell indices covered by the zoom cube."""
    half = 0.5 * HI_REGION_SIZE
    lo = centre - half
    hi = centre + half
    width = BOX_SIZE / BKG_CELLS

    i_lo = int(np.floor(lo / width))
    i_hi = int(np.floor((hi - 1e-12) / width))
    return set(range(i_lo, i_hi + 1))


def generate_background_particles(rng, centre):
    """Generate low-resolution background particles outside the zoom top cells."""
    top_w = BOX_SIZE / BKG_CELLS
    zoom_ids = zoom_top_cell_indices(centre)

    chunks = []
    for ix in range(BKG_CELLS):
        for iy in range(BKG_CELLS):
            for iz in range(BKG_CELLS):
                if ix in zoom_ids and iy in zoom_ids and iz in zoom_ids:
                    continue
                lo = top_w * np.array([ix, iy, iz], dtype=np.float64)
                hi = lo + top_w
                chunks.append(
                    sample_stratified_in_cell(rng, lo, hi, BKG_PARTICLES_AXIS)
                )
    return np.vstack(chunks)


def validate_layout(hi_pos, bkg_pos, centre):
    """Validate occupancy and geometric consistency."""
    all_pos = np.vstack((hi_pos, bkg_pos))
    bkg_occ = occupancy_3d(all_pos, BKG_CELLS)
    if np.any(bkg_occ == 0):
        raise RuntimeError("Found empty top-level background cells.")

    half = 0.5 * HI_REGION_SIZE
    hi_lo = centre - half
    hi_hi = centre + half
    zoom_occ = occupancy_3d(hi_pos, ZOOM_CELLS_PER_AXIS, lo=hi_lo, hi=hi_hi)
    if np.any(zoom_occ != HI_PARTICLES_PER_ZOOM_CELL):
        raise RuntimeError("High-res zoom sub-cell occupancy is not uniform.")

    if np.any((hi_pos < hi_lo) | (hi_pos >= hi_hi)):
        raise RuntimeError("High-res particles found outside zoom cube.")

    bkg_in_hi = np.all((bkg_pos >= hi_lo) & (bkg_pos < hi_hi), axis=1)
    if np.any(bkg_in_hi):
        raise RuntimeError("Background particles found inside high-res zoom cube.")


def write_ic_file(hi_pos, bkg_pos):
    """Write SWIFT-compatible ICs with PartType1 and PartType2."""
    hi_count = hi_pos.shape[0]
    bkg_count = bkg_pos.shape[0]

    hi_mass = np.full(hi_count, RHO_DM * (HI_REGION_SIZE**3) / hi_count)
    bkg_mass_total = RHO_DM * (BOX_SIZE**3 - HI_REGION_SIZE**3)
    bkg_mass = np.full(bkg_count, bkg_mass_total / bkg_count)

    hi_vel = np.zeros_like(hi_pos)
    bkg_vel = np.zeros_like(bkg_pos)

    writer = Writer(
        cosmo_units, np.array([BOX_SIZE, BOX_SIZE, BOX_SIZE]) * Mpc, dimension=3
    )
    writer_any: Any = writer
    writer_any.dark_matter.coordinates = hi_pos * Mpc
    writer_any.dark_matter.velocities = hi_vel * km / s
    writer_any.dark_matter.masses = hi_mass * 1e10 * Msun
    writer.write(OUTPUT_FILE)

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


def main():
    """Generate, validate, and write the IC file."""
    centre = 0.75
    rng = np.random.default_rng(RNG_SEED)

    hi_pos = generate_high_res_particles(rng, centre)
    bkg_pos = generate_background_particles(rng, centre)

    validate_layout(hi_pos, bkg_pos, centre)
    write_ic_file(hi_pos, bkg_pos)


if __name__ == "__main__":
    main()
