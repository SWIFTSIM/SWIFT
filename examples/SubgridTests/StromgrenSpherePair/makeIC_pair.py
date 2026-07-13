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
Two-star variant of ../StromgrenSphere/makeIC.py: places two hand-positioned
stars (instead of one, centered) in an otherwise uniform-density gas box, so
their individually-growing HII regions eventually meet and overlap.

Why this exists: the single-star StromgrenSphere example never exercises two
concurrently-active `radiation_level` regions with geometrically-overlapping
search radii -- exactly the precondition for the cross-region gas race that
the task_lock/task_unlock locking in src/task.c
(task_type_stars_hii_ionization_feedback) protects against. This is the
smallest setup that can reproduce it under real (not synthetic) task-graph
execution.

Star separation defaults to 2x the equilibrium Stromgren radius R_st for the
fiducial 29.7 Msun star at rho=5 atom/cc (~0.0266 kpc, see
stromgren_analytic_check.py's ionizing_photon_rate/stromgren_radius) -- close
enough that sustained D-type expansion brings the two ionization fronts into
contact well before the box's periodic images would (box half-width is kept
>> R_st to avoid the periodic self-interaction pitfall documented in the
project's HomogeneousBox/StromgrenSphere test history).
"""

import h5py
import numpy as np
import argparse
from astropy import units
from astropy import constants


def parse_options():
    usage = "usage: %prog [options] file"
    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument("--rho", type=float, default=5,
                        help="Mean gas density in atom/cm3")
    parser.add_argument("--mass", type=float, default=0.1,
                        help="Gas particle mass in solar mass")
    parser.add_argument("--boxsize", type=float, default=0.16,
                        help="Boxsize in kpc (default gives a half-width of 0.08 "
                             "kpc, ~3x the fiducial 29.7 Msun star's equilibrium "
                             "Stromgren radius at rho=5 atom/cc, ~0.0266 kpc -- "
                             "enough margin to avoid periodic self-interaction "
                             "while both HII regions grow and meet)")
    parser.add_argument("--star_mass", type=float, default=29.7,
                        help="Mass of each star in M_sun")
    parser.add_argument("--n_cells", type=int, default=4,
                        help="Scheduler:max_top_level_cells this IC is meant to be "
                             "run with. Stars are placed at the centers of two "
                             "CONSECUTIVE top-level cells along x (not symmetric "
                             "about the box center, which can land exactly on a "
                             "cell boundary): the radiation task graph only links "
                             "directly-adjacent top-level cells (26-neighbour "
                             "stencil), so this is the only placement that "
                             "guarantees the two stars' regions are graph "
                             "neighbours regardless of search-radius value.")
    parser.add_argument("--star_type", type=str, default="single_star",
                        choices=["single_star", "continuous_IMF", "SSP"])
    parser.add_argument("-o", dest="outputfilename", type=str,
                        default="ICs_stromgren_pair.hdf5")

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
UnitTime = UnitTime_in_cgs * units.s

np.random.seed(1)

rho = opt.rho * constants.m_p / units.cm ** 3
m = opt.mass * units.Msun

L = opt.boxsize * units.kpc
M = rho * L ** 3
N = int(np.ceil((M / m).to(units.dimensionless_unscaled).value))

print("Boxsize                               : {}".format(L.to(units.kpc)))
print("Total gas mass                        : {}".format(M.to(units.Msun)))
print("Number of gas particles                : {}".format(N))

m = m.to(UnitMass).value
L = L.to(UnitLength).value
rho_code = rho.to(UnitMass / UnitLength ** 3).value

pos = np.random.random([N, 3]) * np.array([L, L, L])
vel = np.zeros([N, 3])
mass = np.ones(N) * m
u_therm = np.ones(N)
ids = np.arange(N)
h = np.ones(N) * 3 * L / N ** (1 / 3.0)
rho_arr = np.ones(N) * rho_code

print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))

#####################
# Two hand-placed stars, at the centers of two CONSECUTIVE top-level cells
# along x (see --n_cells help above for why not symmetric-about-center).
#####################
N_star = 2
M_star_val = (opt.star_mass * units.M_sun).to(UnitMass).value
M_star = np.array([M_star_val, M_star_val])

if opt.n_cells < 2:
    raise ValueError("--n_cells must be >= 2 to have two distinct cells.")

cell_width = L / opt.n_cells
mid_cell_lo = opt.n_cells // 2 - 1  # index of the first of the two middle cells
x_a = (mid_cell_lo + 0.5) * cell_width
x_b = (mid_cell_lo + 1.5) * cell_width
yz = L / 2.0
pos_star = np.array([
    [x_a, yz, yz],
    [x_b, yz, yz],
])
separation_code = x_b - x_a

print("Cell width (code unit)                : {:e}".format(cell_width))
print("Star separation (code unit)           : {:e}".format(separation_code))
print("Star positions                        : {}".format(pos_star))
print("Mass of each star (internal units)    : {:e}".format(M_star_val))

vel_star = np.zeros([N_star, 3])
h_star = np.ones(N_star) * 3 * L / N ** (1 / 3.0)

star_type_map = {"single_star": 0, "continuous_IMF": 1, "SSP": 2}
star_type = star_type_map[opt.star_type]
star_particle_type = np.ones(N_star) * star_type
star_id = np.array([N + 1, N + 2])
star_birth_time = np.zeros(N_star)

fileOutput = h5py.File(opt.outputfilename, "w")
print("{} saved.".format(opt.outputfilename))

grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]
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
