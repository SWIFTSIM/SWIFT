#!/usr/bin/env python3
"""Initial conditions: a periodic glass box of gas with one star at its centre.

Every input is physical at the starting scale factor a = 1 / (1 + z). With cosmology SWIFT reads
positions and densities as comoving, so the file holds the box size L / a, the density rho a^3, the
physical internal energy and, for the star, the birth scale factor in BirthTime. Without cosmology
the file is a plain IC and the run starts at the age of the star (run.sh does so).
"""

import argparse

import h5py
import numpy as np
import yaml
from cosmologyTimeline import hubble_rate

UNIT_LENGTH = 3.0856775814913673e21  # 1 kpc in cm
UNIT_MASS = 1.988409870698051e43  # 1e10 Msun in g
UNIT_VELOCITY = 1e5  # km/s in cm/s
MSUN = 1.988409870698051e33  # g
M_PROTON = 1.67262192369e-24  # g
K_BOLTZMANN = 1.380649e-16  # erg/K
GAMMA = 5.0 / 3.0
HYDROGEN_MASS_FRACTION = 0.76


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--glass", required=True, help="Glass cube file (getGlass.sh)")
    p.add_argument(
        "--level", type=int, default=5, help="Gas particle number is (2**level)**3"
    )
    p.add_argument(
        "--rho", type=float, default=1.0, help="Physical gas density in atom/cm^3"
    )
    p.add_argument("--mass", type=float, default=1.0, help="Gas particle mass in Msun")
    p.add_argument(
        "--temperature", type=float, default=100.0, help="Gas temperature in K"
    )
    p.add_argument(
        "--cosmology", type=int, default=1, help="1: expanding box, 0: static box"
    )
    p.add_argument("--redshift", type=float, default=0.0, help="Starting redshift")
    p.add_argument("--star_mass", type=float, default=3.0, help="Star mass in Msun")
    p.add_argument(
        "--star_age",
        type=float,
        default=1e-6,
        help="Age of the star at the start, internal time",
    )
    p.add_argument(
        "--params", default="params.yml", help="Parameter file, for the cosmology"
    )
    p.add_argument("-o", dest="out", default="ICs.hdf5", help="Output file")
    opt = p.parse_args()

    a = 1.0 / (1.0 + opt.redshift) if opt.cosmology else 1.0
    n1 = 2**opt.level
    n = n1**3
    m_g = opt.mass * MSUN
    rho_cgs = opt.rho * M_PROTON
    box_physical = (n * m_g / rho_cgs) ** (1.0 / 3.0)
    mu = 4.0 / (1.0 + 3.0 * HYDROGEN_MASS_FRACTION)  # neutral gas
    u_cgs = K_BOLTZMANN * opt.temperature / ((GAMMA - 1.0) * mu * M_PROTON)

    box = box_physical / a / UNIT_LENGTH
    rho = rho_cgs * a**3 * UNIT_LENGTH**3 / UNIT_MASS
    u = u_cgs / UNIT_VELOCITY**2
    m = m_g / UNIT_MASS
    ms = opt.star_mass * MSUN / UNIT_MASS

    with h5py.File(opt.glass, "r") as g:
        glass_pos = g["PartType0/Coordinates"][:]
        glass_h = g["PartType0/SmoothingLength"][:]
    if len(glass_pos) != n:
        raise RuntimeError(f"{opt.glass} has {len(glass_pos)} particles, expected {n}.")
    pos = (
        (glass_pos - glass_pos.min()) / (glass_pos.max() - glass_pos.min() + 1e-6) * box
    )
    h = glass_h * box

    birth = 0.0
    if opt.cosmology:
        with open(opt.params) as stream:
            cosmology = yaml.safe_load(stream)["Cosmology"]
        # H t << 1 over the age, so ln a is linear in t. SWIFT stores the birth scale factor in
        # float32: round it down, or the star is born after the start.
        a_birth = a * np.exp(-hubble_rate(a, cosmology) * opt.star_age)
        birth = np.float32(a_birth)
        if float(birth) > a_birth:
            birth = np.nextafter(birth, np.float32(0.0))

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
        g.create_dataset("Velocities", data=np.zeros((n, 3)), dtype="f", **kw)
        g.create_dataset("Masses", data=np.full(n, m), dtype="f", **kw)
        g.create_dataset("SmoothingLength", data=h, dtype="f", **kw)
        g.create_dataset("InternalEnergy", data=np.full(n, u), dtype="f", **kw)
        g.create_dataset("ParticleIDs", data=np.arange(n), dtype="L", **kw)
        g.create_dataset("Density", data=np.full(n, rho), dtype="f", **kw)
        s = f.create_group("PartType4")
        s.create_dataset("Coordinates", data=np.full((1, 3), 0.5 * box), dtype="d")
        s.create_dataset("Velocities", data=np.zeros((1, 3)), dtype="f")
        s.create_dataset("Masses", data=[ms], dtype="f")
        s.create_dataset("ParticleIDs", data=[n + 1], dtype="L")
        s.create_dataset("SmoothingLength", data=[3.0 * box / n1], dtype="f")
        s.create_dataset("BirthMass", data=[ms], dtype="f")
        s.create_dataset("BirthTime", data=[birth], dtype="f")
        s.create_dataset(
            "StellarParticleType", data=[0], dtype="i"
        )  # discrete single star
    print(
        f"{opt.out}: {n} gas particles, scale factor {a:.10g}, physical box {box * a * 1e3:.6g} pc, "
        f"star birth {float(birth):.10g}"
    )


if __name__ == "__main__":
    main()
