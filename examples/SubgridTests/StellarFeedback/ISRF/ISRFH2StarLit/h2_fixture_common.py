"""Shared helpers of the H2 line-absorption fixtures: gas box, one or two stars, geometry file."""

import json
import os

import h5py
import numpy as np

UNIT_MASS_CGS = 1.988409870698051e43  # 10^10 solar masses
UNIT_LENGTH_CGS = 3.0856775814913673e21  # kpc
UNIT_VELOCITY_CGS = 1e5  # km/s
MSUN_CGS = 1.988409870698051e33
M_PROTON_CGS = 1.67262192369e-24
KPC_CGS = UNIT_LENGTH_CGS
UNIT_TIME_CGS = UNIT_LENGTH_CGS / UNIT_VELOCITY_CGS


def box_positions(level, mass_msun, rho_hpcm3, glass_dir="."):
    """Return the box size and particle positions in kpc, and the particle number.

    Uses the pre-relaxed glass cube when it is present, a random placement
    otherwise (a random placement is only for tests of the scripts).
    """
    n_part = (2**level) ** 3
    rho = rho_hpcm3 * M_PROTON_CGS
    box = (n_part * mass_msun * MSUN_CGS / rho) ** (1.0 / 3.0) / KPC_CGS
    glass = os.path.join(glass_dir, f"glassCube_{2**level}.hdf5")
    if os.path.exists(glass):
        with h5py.File(glass, "r") as f:
            pos = f["/PartType0/Coordinates"][:, :]
        eps = 1e-6
        pos = (pos - pos.min()) / (pos.max() - pos.min() + eps) * box
        pos = (pos + 0.5 * box / 2**level) % box
    else:
        print(f"No {glass}: random placement, for script tests only.")
        pos = np.random.default_rng(1).random([n_part, 3]) * box
    return box, pos, n_part


def check_away_from_faces(pos, box, label, margin=1e-6):
    """Raise if a coordinate lies within margin * box of a box face.

    A particle on a face can drift out of the box by less than one ulp of the box size
    and be assigned to a cell that does not exist.
    """
    d = np.minimum(pos, box - pos)
    bad = int(np.sum(d < margin * box))
    if bad:
        raise RuntimeError(
            f"{label}: {bad} coordinates lie within {margin:g} box sizes of a box face"
        )


def write_ic(
    filename,
    box,
    pos,
    mass,
    rho,
    vel,
    star_pos,
    star_mass_msun,
    star_type=0,
    star_ids=None,
):
    """Write gas and star particles in the SWIFT IC layout (internal units).

    The star IDs default to n_gas + 0, 1, ... . The ID seeds the random frame of the
    emission pixels, so a star keeps its ID when it is written alone.
    """
    n_gas = len(pos)
    n_star = len(star_pos)
    star_ids = n_gas + np.arange(n_star) if star_ids is None else star_ids
    check_away_from_faces(pos, box, "gas")
    check_away_from_faces(star_pos, box, "star")
    m_star = np.full(n_star, star_mass_msun * MSUN_CGS / UNIT_MASS_CGS)
    h_gas = 1.2348 * (box**3 / n_gas) ** (1.0 / 3.0) * np.ones(n_gas)
    with h5py.File(filename, "w") as f:
        g = f.create_group("/Header")
        g.attrs["BoxSize"] = [box, box, box]
        g.attrs["NumPart_Total"] = [n_gas, 0, 0, 0, n_star, 0]
        g.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
        g.attrs["NumPart_ThisFile"] = [n_gas, 0, 0, 0, n_star, 0]
        g.attrs["Time"] = 0.0
        g.attrs["NumFileOutputsPerSnapshot"] = 1
        g.attrs["MassTable"] = [0.0] * 6
        g.attrs["Flag_Entropy_ICs"] = [0] * 6
        g.attrs["Dimension"] = 3
        u = f.create_group("/Units")
        u.attrs["Unit length in cgs (U_L)"] = UNIT_LENGTH_CGS
        u.attrs["Unit mass in cgs (U_M)"] = UNIT_MASS_CGS
        u.attrs["Unit time in cgs (U_t)"] = UNIT_TIME_CGS
        u.attrs["Unit current in cgs (U_I)"] = 1
        u.attrs["Unit temperature in cgs (U_T)"] = 1
        p = f.create_group("/PartType0")
        p.create_dataset("Coordinates", data=pos, dtype="d")
        p.create_dataset("Velocities", data=vel, dtype="f")
        p.create_dataset("Masses", data=mass, dtype="f")
        p.create_dataset("SmoothingLength", data=h_gas, dtype="f")
        p.create_dataset("InternalEnergy", data=np.ones(n_gas), dtype="f")
        p.create_dataset("ParticleIDs", data=np.arange(n_gas), dtype="L")
        p.create_dataset("Densities", data=rho, dtype="f")
        s = f.create_group("/PartType4")
        s.create_dataset("Coordinates", data=star_pos, dtype="d")
        s.create_dataset("Velocities", data=np.zeros([n_star, 3]), dtype="f")
        s.create_dataset("Masses", data=m_star, dtype="f")
        s.create_dataset("ParticleIDs", data=star_ids, dtype="L")
        s.create_dataset(
            "SmoothingLength", data=np.full(n_star, 3 * h_gas[0]), dtype="f"
        )
        s.create_dataset("BirthMass", data=m_star, dtype="f")
        s.create_dataset("BirthTime", data=np.zeros(n_star), dtype="f")
        s.create_dataset(
            "StellarParticleType", data=np.full(n_star, star_type), dtype="i"
        )


def write_fixture(filename, **entries):
    """Write the geometry file that the check script reads."""
    with open(filename, "w") as f:
        json.dump(entries, f, indent=1)
