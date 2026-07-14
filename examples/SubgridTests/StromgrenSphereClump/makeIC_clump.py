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
Single-star variant of ../StromgrenSphere/makeIC.py: a uniform diffuse
background plus a small, dense gas clump placed near the star along +x
(HEALPix nside=1 pixel 4's center direction). See the README for why this
geometry, and for the clump-sizing rationale.

Clump placement/size are set in parsec directly (--clump_distance_pc,
--clump_radius_pc or --clump_halfangle_deg), not derived from the classical
Stromgren radius -- see the README.
"""

import h5py
import numpy as np
import argparse
from astropy import units
from astropy import constants


def parse_options():
    usage = "usage: %prog [options] file"
    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument("--rho", type=float, default=100,
                        help="Diffuse background gas density in atom/cm3")
    parser.add_argument("--mass", type=float, default=0.001,
                        help="Gas particle mass in solar mass (diffuse and clump)")
    parser.add_argument("--boxsize", type=float, default=0.004,
                        help="Boxsize in kpc")
    parser.add_argument("--star_mass", type=float, default=29.7,
                        help="Mass of the star in M_sun")
    parser.add_argument("--density_factor", type=float, default=100.0,
                        help="Clump density, as a multiple of --rho. Set to "
                             "1 for a clump-free control run.")
    parser.add_argument("--clump_distance_pc", type=float, default=0.35,
                        help="Distance from the star to the clump center, in parsec")
    parser.add_argument("--clump_halfangle_deg", type=float, default=20.0,
                        help="Angular half-size of the clump as seen from the "
                             "star (degrees); sets the clump radius unless "
                             "--clump_radius_pc is given")
    parser.add_argument("--clump_radius_pc", type=float, default=None,
                        help="Clump radius in parsec, overriding "
                             "--clump_halfangle_deg -- use to vary the "
                             "clump's size at fixed distance")
    parser.add_argument("--n_cells", type=int, default=3,
                        help="Scheduler:max_top_level_cells this IC is meant "
                             "to be run with")
    parser.add_argument("--star_type", type=str, default="single_star",
                        choices=["single_star", "continuous_IMF", "SSP"])
    parser.add_argument("-o", dest="outputfilename", type=str,
                        default="ICs_stromgren_clump.hdf5")

    return parser.parse_args()


# Measured with --density_factor 1 (no clump) at the default rho/star_mass/
# mass; not physical constants -- re-measure if those change (see README).
MEASURED_FRONT_REACH_PC = 0.87
MEASURED_FRONT_IONIZED_MASS_MSUN = 7.2


def pick_clump_geometry(opt):
    """Place and size the clump relative to the measured ionization front
    (see README), not the classical Stromgren radius."""
    r_clump = opt.clump_distance_pc * units.pc
    if opt.clump_radius_pc is not None:
        R_clump = opt.clump_radius_pc * units.pc
        halfangle = np.arctan((R_clump / r_clump).decompose())
    else:
        halfangle = np.radians(opt.clump_halfangle_deg)
        R_clump = r_clump * np.tan(halfangle)
    V_clump = (4.0 / 3.0) * np.pi * R_clump**3

    clump_mass = (opt.density_factor * opt.rho / units.cm**3
                  * constants.m_p * V_clump).to(units.Msun)
    demand_fraction = (clump_mass / (MEASURED_FRONT_IONIZED_MASS_MSUN * units.Msun)).decompose()

    print("--- Empirically-calibrated clump sizing ---")
    print(f"Measured front reach (reference)      : {MEASURED_FRONT_REACH_PC} pc")
    print(f"Measured front ionized mass (reference): {MEASURED_FRONT_IONIZED_MASS_MSUN} Msun")
    print(f"Clump center distance from star        : {r_clump.to(units.pc):.4g}")
    print(f"Clump radius                           : {R_clump.to(units.pc):.4g}")
    print(f"Clump implied half-angle from star      : {np.degrees(halfangle):.4g} deg")
    print(f"Clump density                          : {opt.density_factor * opt.rho:.4g} atom/cm3")
    print(f"Clump mass                             : {clump_mass:.4g}")
    print(f"Clump mass / front ionized mass         : {demand_fraction:.2%}")
    if np.degrees(halfangle) > 29.0:
        print("  WARNING: clump half-angle exceeds nside=1 pixel 4's ~29 "
              "degree inscribed-cap margin -- it straddles neighbouring "
              "pixels, weakening the clean single-pixel signal. Decrease "
              "--clump_radius_pc / --clump_halfangle_deg or increase "
              "--clump_distance_pc.")
    if opt.clump_distance_pc + R_clump.to(units.pc).value > MEASURED_FRONT_REACH_PC:
        print("  WARNING: clump extends past the measured front reach -- it "
              "may never be fully reached by the ionization search. Decrease "
              "--clump_distance_pc, --clump_radius_pc, or --clump_halfangle_deg.")
    if demand_fraction < 0.1:
        print("  WARNING: clump mass is small relative to the front's total "
              "budget -- nside=0 vs nside=1 may not show a clearly "
              "distinguishable signal. Increase --density_factor.")
    if demand_fraction > 0.8:
        print("  WARNING: clump mass is very large relative to the front's "
              "total budget -- the clump may starve the diffuse gas almost "
              "completely even at nside=1. Decrease --density_factor.")

    half_width = (opt.boxsize / 2.0) * units.kpc
    if (r_clump + R_clump) > 0.5 * half_width:
        print("  WARNING: clump extends past half the box half-width -- "
              "risks periodic self-interaction or exceeding "
              "Stars:max_HII_search_radius.")

    return r_clump.to(units.kpc).value, R_clump.to(units.kpc).value


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

r_clump_kpc, R_clump_kpc = pick_clump_geometry(opt)

rho = opt.rho * constants.m_p / units.cm ** 3
m = opt.mass * units.Msun

L = opt.boxsize * units.kpc
M_box = rho * L ** 3
N_diffuse_target = int(np.ceil((M_box / m).to(units.dimensionless_unscaled).value))

m_code = m.to(UnitMass).value
L_code = L.to(UnitLength).value
rho_code = rho.to(UnitMass / UnitLength ** 3).value

star_pos = np.array([L_code / 2.0, L_code / 2.0, L_code / 2.0])
# +x is HEALPix nside=1 pixel 4's center direction (see README).
r_clump_code = (r_clump_kpc * units.kpc).to(UnitLength).value
R_clump_code = (R_clump_kpc * units.kpc).to(UnitLength).value
clump_center = star_pos + np.array([r_clump_code, 0.0, 0.0])

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

print(f"Diffuse background particles (after clump carve-out): {N_diffuse} "
      f"(target before carve-out: {N_diffuse_target})")

#####################
# Dense clump, uniform within a sphere of radius R_clump around clump_center.
#####################
mass_clump_code = (
    opt.density_factor * rho * (4.0 / 3.0) * np.pi * (R_clump_kpc * units.kpc) ** 3
).to(UnitMass).value
N_clump = int(np.ceil(mass_clump_code / m_code))

# Uniform sampling within a sphere via rejection (fast enough at this N).
pos_clump = np.empty([N_clump, 3])
n_filled = 0
while n_filled < N_clump:
    n_try = (N_clump - n_filled) * 2 + 16
    cand = (np.random.random([n_try, 3]) * 2.0 - 1.0) * R_clump_code
    r2 = np.sum(cand**2, axis=1)
    cand = cand[r2 <= R_clump_code**2]
    n_take = min(cand.shape[0], N_clump - n_filled)
    pos_clump[n_filled:n_filled + n_take] = cand[:n_take]
    n_filled += n_take
pos_clump += clump_center

print(f"Clump particles                                       : {N_clump}")

#####################
# Combine the two gas populations.
#####################
N = N_diffuse + N_clump
pos = np.vstack([pos_diffuse, pos_clump])
vel = np.zeros([N, 3])
mass = np.ones(N) * m_code
u_therm = np.ones(N)
ids = np.arange(N)

h_diffuse = 3.0 * L_code / N_diffuse_target ** (1.0 / 3.0)
h_clump = 3.0 * (2.0 * R_clump_code) / max(N_clump, 1) ** (1.0 / 3.0)
h = np.concatenate([
    np.ones(N_diffuse) * h_diffuse,
    np.ones(N_clump) * h_clump,
])
rho_arr = np.concatenate([
    np.ones(N_diffuse) * rho_code,
    np.ones(N_clump) * opt.density_factor * rho_code,
])

print(f"Total gas particles                                   : {N}")
print(f"Diffuse smoothing length guess (code unit)            : {h_diffuse:e}")
print(f"Clump smoothing length guess (code unit)               : {h_clump:e}")

#####################
# The star, at the box center.
#####################
N_star = 1
M_star_val = (opt.star_mass * units.M_sun).to(UnitMass).value
M_star = np.array([M_star_val])
pos_star = np.array([star_pos])
vel_star = np.zeros([N_star, 3])
h_star = np.ones(N_star) * h_diffuse

star_type_map = {"single_star": 0, "continuous_IMF": 1, "SSP": 2}
star_type = star_type_map[opt.star_type]
star_particle_type = np.ones(N_star) * star_type
star_id = np.array([N + 1])
star_birth_time = np.zeros(N_star)

print(f"Star position (code unit)                             : {pos_star[0]}")
print(f"Clump center (code unit)                              : {clump_center}")
print(f"Clump-star direction                                  : "
      f"{(clump_center - pos_star[0]) / np.linalg.norm(clump_center - pos_star[0])}")

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
