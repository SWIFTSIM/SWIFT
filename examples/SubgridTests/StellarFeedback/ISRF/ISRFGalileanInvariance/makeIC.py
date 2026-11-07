#!/usr/bin/env python3
"""Initial conditions: a periodic gas box with one star at its centre, at rest or boosted along x.

The gas sits on a perfect cubic lattice (--ic lattice) or on a glass (--ic glass). Gas and star
share the bulk velocity, so the boosted box is the rest box seen from a moving frame.
"""

import argparse

import h5py
import numpy as np

UNIT_LENGTH = 3.0856775814913673e21  # 1 kpc in cm
UNIT_MASS = 1.988409870698051e43  # 1e10 Msun in g
UNIT_VELOCITY = 1e5  # km/s in cm/s
MSUN = 1.988409870698051e33  # g
M_PROTON = 1.67262192369e-24  # g


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--glass",
        required=True,
        help="Glass cube file (getGlass.sh); also sets the smoothing length",
    )
    p.add_argument(
        "--ic",
        choices=["lattice", "glass"],
        default="lattice",
        help="Gas particle layout",
    )
    p.add_argument(
        "--level", type=int, default=5, help="Gas particle number is (2**level)**3"
    )
    p.add_argument("--rho", type=float, default=1e3, help="Gas density in atom/cm^3")
    p.add_argument("--mass", type=float, default=0.1, help="Gas particle mass in Msun")
    p.add_argument("--star_mass", type=float, default=29.7, help="Star mass in Msun")
    p.add_argument(
        "--velocity",
        type=float,
        default=0.0,
        help="Bulk velocity along x of gas and star in km/s",
    )
    p.add_argument("-o", dest="out", default="ICs.hdf5", help="Output file")
    opt = p.parse_args()

    n1 = 2**opt.level
    n = n1**3
    m = opt.mass * MSUN / UNIT_MASS
    rho = opt.rho * M_PROTON / UNIT_MASS * UNIT_LENGTH**3
    box = (n * m / rho) ** (1.0 / 3.0)

    with h5py.File(opt.glass, "r") as g:
        glass_pos = g["PartType0/Coordinates"][:]
        glass_h = g["PartType0/SmoothingLength"][:]
    if len(glass_pos) != n:
        raise RuntimeError(f"{opt.glass} has {len(glass_pos)} particles, expected {n}.")
    if opt.ic == "glass":
        pos = (
            (glass_pos - glass_pos.min())
            / (glass_pos.max() - glass_pos.min() + 1e-6)
            * box
        )
        h = glass_h * box
    else:
        g1 = (np.arange(n1) + 0.5) * box / n1
        pos = np.stack(np.meshgrid(g1, g1, g1, indexing="ij"), axis=-1).reshape(-1, 3)
        h = np.full(n, np.median(glass_h) * box)

    v_bulk = np.array([opt.velocity, 0.0, 0.0])
    ms = opt.star_mass * MSUN / UNIT_MASS
    star_pos = np.full((1, 3), 0.5 * box)
    star_h = 3.0 * box / n1

    kw = dict(compression="gzip", compression_opts=4, shuffle=True)
    with h5py.File(opt.out, "w") as f:
        hd = f.create_group("Header")
        hd.attrs["BoxSize"] = [box, box, box]
        hd.attrs["NumPart_Total"] = [n, 0, 0, 0, 1, 0]
        hd.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
        hd.attrs["NumPart_ThisFile"] = [n, 0, 0, 0, 1, 0]
        hd.attrs["Time"] = 0.0
        hd.attrs["NumFileOutputsPerSnapshot"] = 1
        hd.attrs["MassTable"] = [0.0] * 6
        hd.attrs["Flag_Entropy_ICs"] = [0] * 6
        hd.attrs["Dimension"] = 3
        un = f.create_group("Units")
        un.attrs["Unit length in cgs (U_L)"] = UNIT_LENGTH
        un.attrs["Unit mass in cgs (U_M)"] = UNIT_MASS
        un.attrs["Unit time in cgs (U_t)"] = UNIT_LENGTH / UNIT_VELOCITY
        un.attrs["Unit current in cgs (U_I)"] = 1.0
        un.attrs["Unit temperature in cgs (U_T)"] = 1.0
        g = f.create_group("PartType0")
        g.create_dataset("Coordinates", data=pos, dtype="d", **kw)
        g.create_dataset("Velocities", data=np.tile(v_bulk, (n, 1)), dtype="f", **kw)
        g.create_dataset("Masses", data=np.full(n, m), dtype="f", **kw)
        g.create_dataset("SmoothingLength", data=h, dtype="f", **kw)
        g.create_dataset("InternalEnergy", data=np.ones(n), dtype="f", **kw)
        g.create_dataset("ParticleIDs", data=np.arange(n), dtype="L", **kw)
        g.create_dataset("Densities", data=np.full(n, rho), dtype="f", **kw)
        s = f.create_group("PartType4")
        s.create_dataset("Coordinates", data=star_pos, dtype="d")
        s.create_dataset("Velocities", data=v_bulk.reshape(1, 3), dtype="f")
        s.create_dataset("Masses", data=[ms], dtype="f")
        s.create_dataset("ParticleIDs", data=[n + 1], dtype="L")
        s.create_dataset("SmoothingLength", data=[star_h], dtype="f")
        s.create_dataset("BirthMass", data=[ms], dtype="f")
        s.create_dataset("BirthTime", data=[0.0], dtype="f")
        s.create_dataset(
            "StellarParticleType", data=[0], dtype="i"
        )  # discrete single star
    print(
        f"{opt.out}: {n} gas particles, box {box:.6e} kpc, bulk velocity {opt.velocity} km/s, layout {opt.ic}"
    )


if __name__ == "__main__":
    main()
