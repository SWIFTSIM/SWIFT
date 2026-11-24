#!/usr/bin/env python3
"""Initial conditions: uniform gas in a periodic glass cube. One pulse of light is injected by a test parameter, written to beam_args.txt."""

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


def write_ic(path, pos, h, box, m, rho):
    """Write the gas particles in the SWIFT IC format."""
    n = len(pos)
    kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
    with h5py.File(path, "w") as f:
        hd = f.create_group("Header")
        hd.attrs["BoxSize"] = [box, box, box]
        hd.attrs["NumPart_Total"] = [n, 0, 0, 0, 0, 0]
        hd.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
        hd.attrs["NumPart_ThisFile"] = [n, 0, 0, 0, 0, 0]
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


def parse(extra):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--glass", required=True, help="Glass cube file (getGlass.sh)")
    p.add_argument("--rho", type=float, default=1e3, help="Gas density in atom/cm^3")
    p.add_argument("--mass", type=float, default=0.1, help="Gas particle mass in Msun")
    p.add_argument("-o", dest="out", default="ICs.hdf5", help="Output file")
    extra(p)
    return p.parse_args()


# Direction of the beam, at least 15 degrees from the borders of the
# reception slots (12 slots, default frame), so that the beam fills one slot.
DIRECTION = (-0.3977, -0.6757, -0.6207)
START_DISTANCE_H = 12.0  # pulse centre to the box centre at the start
PULSE_H = 10.0  # side of the cubic pulse
ENERGY_PE = 1e-6  # energy per particle of the pulse, PE band, internal units
ENERGY_LW = 0.25e-6  # LW band


def main():
    opt = parse(
        lambda p: p.add_argument(
            "--energy_scale",
            type=float,
            default=1.0,
            help="Factor on the energy of the pulse",
        )
    )
    pos, h, box, m, rho = glass_box(opt.glass, opt.rho, opt.mass)
    h0 = float(h.mean())
    centre = np.full(3, 0.5 * box)
    write_ic(opt.out, pos, h, box, m, rho)
    d = np.asarray(DIRECTION) / np.linalg.norm(DIRECTION)
    c = centre - START_DISTANCE_H * h0 * d
    lo = (c - 0.5 * PULSE_H * h0).astype(np.float32).astype(np.float64)
    hi = (c + 0.5 * PULSE_H * h0).astype(np.float32).astype(np.float64)
    value = ",".join(
        "%.9g" % v
        for v in [
            *lo,
            *hi,
            *d,
            ENERGY_PE * opt.energy_scale,
            ENERGY_LW * opt.energy_scale,
        ]
    )
    open("beam_args.txt", "w").write(
        f"-P\nGEARFeedback:ISRF_packet_test_beam_1:[{value}]\n"
    )
    geometry = dict(
        box_kpc=box,
        h0_kpc=h0,
        centre_kpc=centre.tolist(),
        direction=d.tolist(),
    )
    json.dump(geometry, open("geometry.json", "w"), indent=1)
    print(json.dumps(geometry))


if __name__ == "__main__":
    main()
