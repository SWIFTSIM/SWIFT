################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
"""
Reproduction of Smith et al. 2021 (MNRAS 506, 3882), Figure 2: four ionizing
sources at the centre of a uniform 100 cm^-3 background, plus a 10 pc-radius
clump of 1e4 cm^-3 gas centered 20 pc away along +x. See the README for the
source values and how this differs from ../StromgrenSphereClump (which uses
a geometry calibrated against this code's own measured ionization front
instead of the paper's absolute scale).
"""

import h5py
import numpy as np
import argparse
from astropy import units
from astropy import constants


def parse_options():
    usage = "usage: %prog [options] file"
    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument(
        "--rho",
        type=float,
        default=100,
        help="Diffuse background gas density in atom/cm3",
    )
    parser.add_argument(
        "--mass",
        type=float,
        default=4.0,
        help="Gas particle mass in solar mass (diffuse and clump). "
        "4.0 is Hu et al. 2017's coarsest convergence-test "
        "resolution; see README for the finer levels (0.5, "
        "0.0625, 0.0078125) and why they need a lot more particles.",
    )
    parser.add_argument("--boxsize", type=float, default=0.1, help="Boxsize in kpc")
    parser.add_argument(
        "--star_mass",
        type=float,
        default=19.2,
        help="Mass of each of the 4 stars in M_sun "
        "(Q_H ~ 2.5e48/s per star, matching Smith "
        "et al.'s 4 sources)",
    )
    parser.add_argument(
        "--n_stars",
        type=int,
        default=4,
        help="Number of ionizing sources clustered at the box center",
    )
    parser.add_argument(
        "--density_factor",
        type=float,
        default=100.0,
        help="Clump density, as a multiple of --rho",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=1000.0,
        help="Diffuse background temperature in K. The clump "
        "is set to temperature/density_factor so both "
        "populations start in pressure equilibrium "
        "(n*T equal on both sides) -- must match "
        "SPH:initial_temperature: 0 in params.yml, or "
        "that parameter overwrites these values.",
    )
    parser.add_argument(
        "--clump_distance_pc",
        type=float,
        default=20.0,
        help="Distance from the star to the clump center, in parsec",
    )
    parser.add_argument(
        "--clump_radius_pc", type=float, default=10.0, help="Clump radius in parsec"
    )
    parser.add_argument(
        "--n_cells",
        type=int,
        default=3,
        help="Scheduler:max_top_level_cells this IC is meant " "to be run with",
    )
    parser.add_argument(
        "--star_type",
        type=str,
        default="single_star",
        choices=["single_star", "continuous_IMF", "SSP"],
    )
    parser.add_argument(
        "-o", dest="outputfilename", type=str, default="ICs_stromgren_clump.hdf5"
    )

    return parser.parse_args()


opt = parse_options()

UnitMass_in_cgs = 1.988409870698051e43
UnitLength_in_cgs = 3.0856775814913673e21
UnitVelocity_in_cgs = 1e5
UnitCurrent_in_cgs = 1
UnitTemp_in_cgs = 1
UnitTime_in_cgs = UnitLength_in_cgs / UnitVelocity_in_cgs

UnitMass = UnitMass_in_cgs * units.g
UnitLength = UnitLength_in_cgs * units.cm

np.random.seed(1)

rho = opt.rho * constants.m_p / units.cm**3
m = opt.mass * units.Msun

L = opt.boxsize * units.kpc
M_box = rho * L**3
N_diffuse_target = int(np.ceil((M_box / m).to(units.dimensionless_unscaled).value))

m_code = m.to(UnitMass).value
L_code = L.to(UnitLength).value
rho_code = rho.to(UnitMass / UnitLength**3).value

star_pos = np.array([L_code / 2.0, L_code / 2.0, L_code / 2.0])
# +x is HEALPix nside=1 pixel 4's center direction (see README).
r_clump_code = (opt.clump_distance_pc * units.pc).to(UnitLength).value
R_clump_code = (opt.clump_radius_pc * units.pc).to(UnitLength).value
clump_center = star_pos + np.array([r_clump_code, 0.0, 0.0])

clump_mass = (
    opt.density_factor
    * rho
    * (4.0 / 3.0)
    * np.pi
    * (opt.clump_radius_pc * units.pc) ** 3
).to(units.Msun)
print(
    f"Clump mass: {clump_mass:.4g}, clump particles: "
    f"~{int(clump_mass / (opt.mass * units.Msun))}"
)

#####################
# Diffuse background, with the clump's volume carved out.
#####################
pos_diffuse = np.random.random([N_diffuse_target, 3]) * np.array(
    [L_code, L_code, L_code]
)
d_to_clump = np.sqrt(np.sum((pos_diffuse - clump_center) ** 2, axis=1))
keep = d_to_clump > R_clump_code
pos_diffuse = pos_diffuse[keep]
N_diffuse = pos_diffuse.shape[0]

print(
    f"Diffuse background particles (after clump carve-out): {N_diffuse} "
    f"(target before carve-out: {N_diffuse_target})"
)

#####################
# Dense clump, uniform within a sphere of radius R_clump around clump_center.
#####################
mass_clump_code = (
    (
        opt.density_factor
        * rho
        * (4.0 / 3.0)
        * np.pi
        * (opt.clump_radius_pc * units.pc) ** 3
    )
    .to(UnitMass)
    .value
)
N_clump = int(np.ceil(mass_clump_code / m_code))

pos_clump = np.empty([N_clump, 3])
n_filled = 0
while n_filled < N_clump:
    n_try = (N_clump - n_filled) * 2 + 16
    cand = (np.random.random([n_try, 3]) * 2.0 - 1.0) * R_clump_code
    r2 = np.sum(cand**2, axis=1)
    cand = cand[r2 <= R_clump_code**2]
    n_take = min(cand.shape[0], N_clump - n_filled)
    pos_clump[n_filled : n_filled + n_take] = cand[:n_take]
    n_filled += n_take
pos_clump += clump_center

print(f"Clump particles                                       : {N_clump}")

#####################
# Combine the two gas populations. Internal energy is set per-population so
# the clump starts in pressure equilibrium (n*T equal) with the diffuse
# background, not at the same temperature -- see u_from_temperature().
#####################
UnitVelocity = UnitVelocity_in_cgs * units.cm / units.s


def u_from_temperature(T_kelvin, X_H=0.752, gamma=5.0 / 3.0):
    """u = k_B*T / (m_p*(gamma-1)*mu) in code units; mu for neutral gas
    (T stays well below the 1e4 K ionization threshold for both
    populations here)."""
    mu = 4.0 / (1.0 + 3.0 * X_H)
    u_cgs = (constants.k_B * (T_kelvin * units.K)) / (
        constants.m_p * (gamma - 1.0) * mu
    )
    return (u_cgs / UnitVelocity**2).decompose().value


u_diffuse_code = u_from_temperature(opt.temperature)
u_clump_code = u_from_temperature(opt.temperature / opt.density_factor)
print(
    f"Diffuse temperature: {opt.temperature:.4g} K, clump temperature: "
    f"{opt.temperature / opt.density_factor:.4g} K (pressure-equilibrium pair)"
)

N = N_diffuse + N_clump
pos = np.vstack([pos_diffuse, pos_clump])
vel = np.zeros([N, 3])
mass = np.ones(N) * m_code
u_therm = np.concatenate(
    [
        np.ones(N_diffuse) * u_diffuse_code,
        np.ones(N_clump) * u_clump_code,
    ]
)
ids = np.arange(N)

h_diffuse = 3.0 * L_code / N_diffuse_target ** (1.0 / 3.0)
h_clump = 3.0 * (2.0 * R_clump_code) / max(N_clump, 1) ** (1.0 / 3.0)
h = np.concatenate(
    [
        np.ones(N_diffuse) * h_diffuse,
        np.ones(N_clump) * h_clump,
    ]
)
rho_arr = np.concatenate(
    [
        np.ones(N_diffuse) * rho_code,
        np.ones(N_clump) * opt.density_factor * rho_code,
    ]
)

print(f"Total gas particles                                   : {N}")


#####################
# The N_stars sources, clustered at the box center. SWIFT aborts if two gas
# particles share a position and warns for gravity particles (engine.c), so
# the stars are spread on a sphere far smaller than the gas resolution
# (a fraction of h_diffuse) instead of being placed at the exact same point.
#####################
def sphere_points(n, radius):
    """n points ~evenly spread on a sphere of the given radius (n=1: origin).
    Uses open-interval y sampling (never +/-1) so no point falls exactly on
    an axis -- avoids landing back on a cell boundary aligned with it."""
    if n == 1:
        return np.zeros((1, 3))
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    i = np.arange(n)
    y = 1.0 - 2.0 * (i + 0.5) / n
    r = np.sqrt(np.maximum(0.0, 1.0 - y * y))
    theta = golden_angle * i
    return radius * np.stack([r * np.cos(theta), y, r * np.sin(theta)], axis=1)


N_star = opt.n_stars
M_star_val = (opt.star_mass * units.M_sun).to(UnitMass).value
M_star = np.ones(N_star) * M_star_val
pos_star = star_pos + sphere_points(N_star, 1e-3 * h_diffuse)
vel_star = np.zeros([N_star, 3])
h_star = np.ones(N_star) * h_diffuse

star_type_map = {"single_star": 0, "continuous_IMF": 1, "SSP": 2}
star_type = star_type_map[opt.star_type]
star_particle_type = np.ones(N_star) * star_type
star_id = np.arange(N + 1, N + 1 + N_star)
star_birth_time = np.zeros(N_star)

print(f"Stars ({N_star}) clustered around (code unit)               : {star_pos}")
print(f"Clump center (code unit)                              : {clump_center}")

fileOutput = h5py.File(opt.outputfilename, "w")
print("{} saved.".format(opt.outputfilename))

grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L_code, L_code, L_code]
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, N_star, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, N_star, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = UnitLength_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = UnitMass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = UnitTime_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = UnitCurrent_in_cgs
grp.attrs["Unit temperature in cgs (U_T)"] = UnitTemp_in_cgs

grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=vel, dtype="f")
grp.create_dataset("Masses", data=mass, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u_therm, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("Densities", data=rho_arr, dtype="f")

grp = fileOutput.create_group("/PartType4")
grp.create_dataset("Coordinates", data=pos_star, dtype="d")
grp.create_dataset("Velocities", data=vel_star, dtype="f")
grp.create_dataset("Masses", data=M_star, dtype="f")
grp.create_dataset("ParticleIDs", data=star_id, dtype="L")
grp.create_dataset("SmoothingLength", data=h_star, dtype="f")
grp.create_dataset("BirthMass", data=M_star, dtype="f")
grp.create_dataset("BirthTime", data=star_birth_time, dtype="f")
grp.create_dataset("StellarParticleType", data=star_particle_type, dtype="i")

fileOutput.close()
