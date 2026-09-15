################################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
#               2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

import os
import h5py
import numpy as np
import argparse
from astropy import units
from astropy import constants


def parse_options():
    parser = argparse.ArgumentParser(
        description="Uniform periodic glass box with one star at its centre; "
        "optional uniform bulk velocity (gas and star) and extra star velocity."
    )
    parser.add_argument(
        "--rho", type=float, default=1e3, help="Mean gas density in atom/cm3"
    )
    parser.add_argument(
        "--mass", type=float, default=0.1, help="Gas particle mass in solar mass"
    )
    parser.add_argument(
        "--level", type=int, default=5, help="Resolution level: N = (2**l)**3"
    )
    parser.add_argument(
        "--star_mass", type=float, default=29.7, help="Mass of the star in M_sun"
    )
    parser.add_argument(
        "--star_type",
        type=str,
        default="single_star",
        choices=["single_star", "continuous_IMF", "SSP"],
        help="Type of the star for the GEAR model.",
    )
    parser.add_argument(
        "--bulk-velocity",
        type=float,
        nargs=3,
        default=[0.0, 0.0, 0.0],
        metavar=("VX", "VY", "VZ"),
        help="Uniform velocity (km/s) given to every gas particle and to the star.",
    )
    parser.add_argument(
        "--star-velocity",
        type=float,
        nargs=3,
        default=[0.0, 0.0, 0.0],
        metavar=("VX", "VY", "VZ"),
        help="Velocity (km/s) of the star relative to the gas, added to --bulk-velocity.",
    )
    parser.add_argument(
        "-o",
        dest="outputfilename",
        type=str,
        default="box.hdf5",
        help="output filename",
    )
    return parser.parse_args()


opt = parse_options()

UnitMass_in_cgs = 1.988409870698051e43  # 10^10 M_sun in grams
UnitLength_in_cgs = 3.0856775814913673e21  # kpc in centimeters
UnitVelocity_in_cgs = 1e5  # km/s in centimeters per second
UnitCurrent_in_cgs = 1  # Amperes
UnitTemp_in_cgs = 1  # Kelvin
UnitTime_in_cgs = UnitLength_in_cgs / UnitVelocity_in_cgs

UnitMass = UnitMass_in_cgs * units.g
UnitLength = UnitLength_in_cgs * units.cm

N = (2**opt.level) ** 3
rho = opt.rho * constants.m_p / units.cm**3
m = opt.mass * units.Msun
L = ((N * m / rho) ** (1 / 3.0)).to(UnitLength).value
m = m.to(UnitMass).value
rho = rho.to(UnitMass / UnitLength**3).value

glass_n = 2**opt.level
glass_filename = f"glassCube_{glass_n}.hdf5"
if not os.path.exists(glass_filename):
    raise RuntimeError(f"Missing {glass_filename}. Run ./getGlass.sh {glass_n} first.")
with h5py.File(glass_filename, "r") as glass:
    pos = glass["/PartType0/Coordinates"][:, :]
    h_glass = glass["/PartType0/SmoothingLength"][:]
if pos.shape[0] != N:
    raise RuntimeError(f"{glass_filename} has {pos.shape[0]} particles, expected {N}.")
eps = 1e-6
pos = (pos - pos.min()) / (pos.max() - pos.min() + eps) * L
h = h_glass * L

# Velocities are in km/s, which is the internal velocity unit.
v_bulk = np.array(opt.bulk_velocity, dtype=float)
v_star = v_bulk + np.array(opt.star_velocity, dtype=float)
vel = np.tile(v_bulk, (N, 1))

print(f"Boxsize (internal units)              : {L:.6e}")
print(f"Number of gas particles               : {N}")
print(f"Inter-particle distance (internal)    : {L / N ** (1 / 3.0):.6e}")
print(f"Gas bulk velocity (km/s)              : {v_bulk}")
print(f"Star velocity (km/s)                  : {v_star}")

pos_star = np.ones([1, 3]) * L / 2
vel_star = v_star.reshape(1, 3)
M_star = [(opt.star_mass * units.M_sun).to(UnitMass).value]
h_star = np.ones(1) * 3 * L / N ** (1 / 3.0)
star_type = {"single_star": 0, "continuous_IMF": 1, "SSP": 2}[opt.star_type]

with h5py.File(opt.outputfilename, "w") as f:
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [L, L, L]
    grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 1, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 1, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFileOutputsPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["Dimension"] = 3

    grp = f.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = UnitLength_in_cgs
    grp.attrs["Unit mass in cgs (U_M)"] = UnitMass_in_cgs
    grp.attrs["Unit time in cgs (U_t)"] = UnitTime_in_cgs
    grp.attrs["Unit current in cgs (U_I)"] = UnitCurrent_in_cgs
    grp.attrs["Unit temperature in cgs (U_T)"] = UnitTemp_in_cgs

    grp = f.create_group("/PartType0")
    grp.create_dataset("Coordinates", data=pos, dtype="d")
    grp.create_dataset("Velocities", data=vel, dtype="f")
    grp.create_dataset("Masses", data=np.ones(N) * m, dtype="f")
    grp.create_dataset("SmoothingLength", data=h, dtype="f")
    grp.create_dataset("InternalEnergy", data=np.ones(N), dtype="f")
    grp.create_dataset("ParticleIDs", data=np.arange(N), dtype="L")
    grp.create_dataset("Densities", data=np.ones(N) * rho, dtype="f")

    grp = f.create_group("/PartType4")
    grp.create_dataset("Coordinates", data=pos_star, dtype="d")
    grp.create_dataset("Velocities", data=vel_star, dtype="f")
    grp.create_dataset("Masses", data=M_star, dtype="f")
    grp.create_dataset("ParticleIDs", data=[N + 1], dtype="L")
    grp.create_dataset("SmoothingLength", data=h_star, dtype="f")
    grp.create_dataset("BirthMass", data=M_star, dtype="f")
    grp.create_dataset("BirthTime", data=np.zeros(1), dtype="f")
    grp.create_dataset("StellarParticleType", data=np.ones(1) * star_type, dtype="i")

print(f"{opt.outputfilename} saved.")
