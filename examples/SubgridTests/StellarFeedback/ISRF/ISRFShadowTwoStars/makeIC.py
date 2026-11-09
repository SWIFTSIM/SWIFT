#!/usr/bin/env python3
"""Initial conditions: one or two stars and an absorbing clump in a non-periodic cube of glass. The clump is the gas of higher metallicity."""

import argparse
import json

import h5py
import numpy as np

import checkShadow

UNIT_LENGTH = 3.0856775814913673e21  # 1 kpc in cm
UNIT_MASS = 1.988409870698051e43  # 1e10 Msun in g
UNIT_VELOCITY = 1e5  # km/s in cm/s
MSUN = 1.988409870698051e33  # g
M_PROTON = 1.67262192369e-24  # g

# Constants of src/feedback/GEAR/radiation.h (dust cross-section of the PE band,
# mean molecular weight of the H nuclei, Grackle solar metal fraction).
SIGMA_D_PE_CGS = 9e-22
MU_H = 1.4
Z_SOLAR_GRACKLE = 0.01295
AMBIENT = 1e-2  # ambient metals in units of the table's solar abundances

# Lengths in units of h0, the mean smoothing length of the gas.
CLUMP_RADIUS_H = 12.0
STAR_DISTANCE_H = 24.0  # star to clump centre: shadow half angle of 30 degrees
ZONE_MAX_H = 20.0  # outer radius of the zones, from the clump centre
MARGIN_H = 4.0  # zones and stars to the faces of the box
EDGE_H = 3.0  # shadow edge that the zones leave out
SKIRT_H = 3.9  # metal skirt of the smoothed clump
TAU = 10.0  # PE optical depth across the clump diameter


def tiled_cube(glass_file, rho_per_cm3, mass_msun, side_h):
    """Tile the periodic glass cube and keep the particles of a cube of side_h smoothing lengths.

    Return positions, smoothing lengths, h0, side, particle mass and density in internal units.
    """
    with h5py.File(glass_file, "r") as g:
        pos = g["PartType0/Coordinates"][:] % 1.0
        h = g["PartType0/SmoothingLength"][:]
    m = mass_msun * MSUN / UNIT_MASS
    rho = rho_per_cm3 * M_PROTON / UNIT_MASS * UNIT_LENGTH**3
    tile = (len(pos) * m / rho) ** (1.0 / 3.0)
    h0 = float(h.mean()) * tile
    side = side_h * h0
    parts, hs = [], []
    n_tiles = int(np.ceil(side / tile))
    for i in range(n_tiles):
        for j in range(n_tiles):
            for k in range(n_tiles):
                p = (pos + np.array([i, j, k])) * tile
                keep = np.all(p < side, axis=1)
                parts.append(p[keep])
                hs.append(h[keep] * tile)
    return np.concatenate(parts), np.concatenate(hs), h0, side, m, rho


def write_ic(path, pos, h, box, m, rho, u, star_pos, star_mass_msun, metals):
    """Write gas and single-star particles in the SWIFT IC format."""
    n, ns = len(pos), len(star_pos)
    kw = dict(compression="gzip", compression_opts=4, shuffle=True, chunks=True)
    with h5py.File(path, "w") as f:
        hd = f.create_group("Header")
        hd.attrs["BoxSize"] = [box, box, box]
        hd.attrs["NumPart_Total"] = [n, 0, 0, 0, ns, 0]
        hd.attrs["NumPart_Total_HighWord"] = [0] * 6
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
        g.create_dataset("InternalEnergy", data=np.full(n, u), dtype="f", **kw)
        g.create_dataset("ParticleIDs", data=np.arange(n), dtype="L", **kw)
        g.create_dataset("Densities", data=np.full(n, rho), dtype="f", **kw)
        g.create_dataset("MetalMassFraction", data=metals, dtype="d", **kw)
        ms = star_mass_msun * MSUN / UNIT_MASS
        s = f.create_group("PartType4")
        s.create_dataset("Coordinates", data=star_pos, dtype="d")
        s.create_dataset("Velocities", data=np.zeros((ns, 3)), dtype="f")
        s.create_dataset("Masses", data=np.full(ns, ms), dtype="f")
        s.create_dataset("ParticleIDs", data=np.arange(n, n + ns), dtype="L")
        s.create_dataset(
            "SmoothingLength", data=np.full(ns, 3.0 * float(np.median(h))), dtype="f"
        )
        s.create_dataset("BirthMass", data=np.full(ns, ms), dtype="f")
        s.create_dataset("BirthTime", data=np.zeros(ns), dtype="f")
        s.create_dataset(
            "StellarParticleType", data=np.zeros(ns, dtype=np.int32), dtype="i"
        )  # discrete single star


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--glass", required=True, help="Glass cube file (getGlass.sh)")
    p.add_argument(
        "--table", required=True, help="Radiation table, for the solar abundances"
    )
    p.add_argument(
        "--stars",
        type=int,
        choices=(1, 2),
        default=1,
        help="Number of stars: 1 on the x axis of the clump, 2 on its x and y axes",
    )
    p.add_argument(
        "--clump",
        type=int,
        default=1,
        help="1: with the clump, 0: ambient metals everywhere",
    )
    p.add_argument("--rho", type=float, default=1e3, help="Gas density in atom/cm^3")
    p.add_argument("--mass", type=float, default=0.8, help="Gas particle mass in Msun")
    p.add_argument("--star_mass", type=float, default=29.7, help="Star mass in Msun")
    p.add_argument(
        "--internal_energy",
        type=float,
        default=0.1,
        help="Gas internal energy in (km/s)^2",
    )
    p.add_argument("-o", dest="out", default="ICs_isrf_shadow.hdf5", help="Output file")
    o = p.parse_args()

    centre_h = ZONE_MAX_H + MARGIN_H
    side_h = centre_h + STAR_DISTANCE_H + MARGIN_H
    pos, h, h0, box, m, rho = tiled_cube(o.glass, o.rho, o.mass, side_h)
    centre = np.array([centre_h, centre_h, 0.5 * side_h]) * h0
    stars = np.array([centre + STAR_DISTANCE_H * h0 * e for e in np.eye(3)[: o.stars]])
    with h5py.File(o.table, "r") as t:
        solar = np.asarray(t["Data"].attrs["SolarMassAbundances"], dtype=np.float64)
    rho_cgs = rho * UNIT_MASS / UNIT_LENGTH**3
    kappa_per_z = (
        SIGMA_D_PE_CGS / Z_SOLAR_GRACKLE / (MU_H * M_PROTON) * rho_cgs * UNIT_LENGTH
    )
    z_ambient = AMBIENT * solar[-1]
    z_clump = TAU / (kappa_per_z * 2.0 * CLUMP_RADIUS_H * h0)
    inside = np.linalg.norm(pos - centre, axis=1) < CLUMP_RADIUS_H * h0
    scale = np.where(inside & (o.clump == 1), z_clump / z_ambient, 1.0)
    metals = (AMBIENT * solar)[None, :] * scale[:, None]
    write_ic(o.out, pos, h, box, m, rho, o.internal_energy, stars, o.star_mass, metals)

    geo = dict(
        box_kpc=box,
        box_h0=side_h,
        h0_kpc=h0,
        centre_kpc=centre.tolist(),
        stars_kpc=stars.tolist(),
        clump_radius_h=CLUMP_RADIUS_H,
        tau_diameter=TAU,
        zone_radius_max_h=ZONE_MAX_H,
        edge_h=EDGE_H,
        skirt_h=SKIRT_H,
        rho0_internal=rho,
        kappa_ambient_per_kpc=kappa_per_z * z_ambient,
        kappa_ambient_per_h0=kappa_per_z * z_ambient * h0,
        z_ambient=z_ambient,
        z_clump=z_clump,
        n_gas=len(pos),
        n_clump=int((inside & (o.clump == 1)).sum()),
    )
    zones, _ = checkShadow.zone_masks((pos - centre) / h0, geo)
    geo["zone_counts"] = {k: int(v.sum()) for k, v in zones.items()}
    in_zone = pos[np.any(list(zones.values()), axis=0)] / h0
    geo["min_face_distance_zone_h"] = float(np.minimum(in_zone, side_h - in_zone).min())
    geo["min_face_distance_star_h"] = float(
        np.minimum(stars / h0, side_h - stars / h0).min()
    )
    json.dump(geo, open("geometry.json", "w"), indent=1)
    print(
        json.dumps(
            {
                k: geo[k]
                for k in (
                    "n_gas",
                    "n_clump",
                    "zone_counts",
                    "min_face_distance_zone_h",
                    "min_face_distance_star_h",
                )
            }
        )
    )
    if min(geo["zone_counts"].values()) < checkShadow.MIN_ZONE:
        raise SystemExit(
            f"A zone holds fewer than {checkShadow.MIN_ZONE} particles: the check cannot score it."
        )
    if (
        min(geo["min_face_distance_zone_h"], geo["min_face_distance_star_h"])
        < MARGIN_H - 1e-6
    ):
        raise SystemExit(
            "A zone particle or a star is closer to a face than the margin."
        )


if __name__ == "__main__":
    main()
