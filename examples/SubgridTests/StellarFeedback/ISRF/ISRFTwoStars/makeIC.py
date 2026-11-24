#!/usr/bin/env python3
"""Initial conditions: two identical stars on the x axis, centred in a periodic glass cube of uniform gas."""

import argparse
import json

import h5py
import numpy as np

UNIT_LENGTH = 3.0856775814913673e21  # 1 kpc in cm
UNIT_MASS = 1.988409870698051e43  # 1e10 Msun in g
UNIT_VELOCITY = 1e5  # km/s in cm/s
MSUN = 1.988409870698051e33  # g
M_PROTON = 1.67262192369e-24  # g


def glass_box(glass_file, rho_per_cm3, mass_msun):
    """Scale a glass cube to the requested density and gas mass.

    Return positions, smoothing lengths, box side, gas particle mass and
    density, all in internal units.
    """
    with h5py.File(glass_file, "r") as g:
        pos = g["PartType0/Coordinates"][:]
        h = g["PartType0/SmoothingLength"][:]
    n = len(pos)
    m = mass_msun * MSUN / UNIT_MASS
    rho = rho_per_cm3 * M_PROTON / UNIT_MASS * UNIT_LENGTH**3
    box = (n * m / rho) ** (1.0 / 3.0)
    pos = (pos - pos.min()) / (pos.max() - pos.min() + 1e-6) * box
    return pos, h * box, box, m, rho


def write_ic(path, pos, h, box, m, rho, star_pos=None, star_mass_msun=0.0, metals=None):
    """Write gas, and optionally single-star particles, in the SWIFT IC format."""
    n = len(pos)
    ns = 0 if star_pos is None else len(star_pos)
    kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
    with h5py.File(path, "w") as f:
        hd = f.create_group("Header")
        hd.attrs["BoxSize"] = [box, box, box]
        hd.attrs["NumPart_Total"] = [n, 0, 0, 0, ns, 0]
        hd.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
        hd.attrs["NumPart_ThisFile"] = [n, 0, 0, 0, ns, 0]
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
        g.create_dataset("Velocities", data=np.zeros((n, 3)), dtype="f", **kw)
        g.create_dataset("Masses", data=np.full(n, m), dtype="f", **kw)
        g.create_dataset("SmoothingLength", data=h, dtype="f", **kw)
        g.create_dataset("InternalEnergy", data=np.ones(n), dtype="f", **kw)
        g.create_dataset("ParticleIDs", data=np.arange(n), dtype="L", **kw)
        g.create_dataset("Densities", data=np.full(n, rho), dtype="f", **kw)
        if metals is not None:
            g.create_dataset("MetalMassFraction", data=metals, dtype="d", **kw)
        if ns:
            ms = star_mass_msun * MSUN / UNIT_MASS
            s = f.create_group("PartType4")
            s.create_dataset("Coordinates", data=star_pos, dtype="d")
            s.create_dataset("Velocities", data=np.zeros((ns, 3)), dtype="f")
            s.create_dataset("Masses", data=np.full(ns, ms), dtype="f")
            s.create_dataset("ParticleIDs", data=np.arange(n, n + ns), dtype="L")
            s.create_dataset("SmoothingLength", data=np.full(ns, 3.0 * box / n ** (1.0 / 3.0)), dtype="f")
            s.create_dataset("BirthMass", data=np.full(ns, ms), dtype="f")
            s.create_dataset("BirthTime", data=np.zeros(ns), dtype="f")
            s.create_dataset("StellarParticleType", data=np.zeros(ns, dtype=np.int32), dtype="i")  # discrete single star


def parse(extra):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--glass", required=True, help="Glass cube file (getGlass.sh)")
    p.add_argument("--rho", type=float, default=1e3, help="Gas density in atom/cm^3")
    p.add_argument("--mass", type=float, default=0.1, help="Gas particle mass in Msun")
    p.add_argument("-o", dest="out", default="ICs.hdf5", help="Output file")
    extra(p)
    return p.parse_args()


def main():
    def extra(p):
        p.add_argument("--star_mass", type=float, default=29.7, help="Star mass in Msun")
        p.add_argument("--distance_pc", type=float, default=3.0, help="Distance between the stars in pc")

    opt = parse(extra)
    pos, h, box, m, rho = glass_box(opt.glass, opt.rho, opt.mass)
    half = 0.5 * 1e-3 * opt.distance_pc  # kpc
    star = np.array([[0.5 * box - half, 0.5 * box, 0.5 * box], [0.5 * box + half, 0.5 * box, 0.5 * box]])
    write_ic(opt.out, pos, h, box, m, rho, star, opt.star_mass)
    print(json.dumps(dict(n_gas=len(pos), box_kpc=box, h_mean_kpc=float(h.mean()), distance_h=2 * half / float(h.mean()))))


if __name__ == "__main__":
    main()
