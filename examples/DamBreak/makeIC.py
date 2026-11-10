###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 SWIFT developers.
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
##############################################################################

import argparse

import h5py
import numpy as np

# Generates a swift IC file for a quasi-2D dam break in a closed tank: a
# hydrostatic column of Tillotson water resting against the x_min wall
# collapses under constant gravity (pointing along -z). All in SI units.

parser = argparse.ArgumentParser(description="Dam break initial conditions")
parser.add_argument(
    "-n",
    "--n_part",
    type=int,
    default=1000000,
    help="Target number of particles (rounded to fill a cubic lattice)",
)
args = parser.parse_args()

# Parameters
column_width = 0.5  # Column size along x
column_depth = 0.25  # Approximate column (and tank) size along y
column_height = 1.0  # Approximate column size along z
tank_length = 2.0  # Box size along x
tank_height = 2.0  # Box size along z
g = 5.0e3  # Gravitational acceleration, must match ConstantPotential:g_cgs
u_init = 1.0  # Specific internal energy (cold: negligible thermal pressure)
matID = 102  # Tillotson water
fileOutputName = "dam_break.hdf5"

# Tillotson water parameters (see src/equation_of_state/planetary/tillotson.h)
rho_0 = 998.0
A = 2.18e9
B = 1.325e10
# ---------------------------------------------------

# Cubic lattice, half a spacing away from the walls. The width is matched
# exactly; the depth and height are rounded to a whole number of spacings.
dx_target = np.cbrt(column_width * column_depth * column_height / args.n_part)
N_x = max(1, round(column_width / dx_target))
dx = column_width / N_x
N_y = max(1, round(column_depth / dx))
N_z = max(1, round(column_height / dx))
tank_depth = N_y * dx
column_height = N_z * dx
numPart = N_x * N_y * N_z
x = (np.arange(N_x) + 0.5) * dx
y = (np.arange(N_y) + 0.5) * dx
z = (np.arange(N_z) + 0.5) * dx
A2_coords = np.stack(np.meshgrid(x, y, z, indexing="ij"), axis=-1).reshape(-1, 3)

# Hydrostatic pressure, inverted for the density with the cold compressed
# Tillotson pressure P = A mu + B mu^2, where mu = rho / rho_0 - 1
A1_P = rho_0 * g * (column_height - A2_coords[:, 2])
A1_mu = (-A + np.sqrt(A * A + 4.0 * B * A1_P)) / (2.0 * B)
A1_rho = rho_0 * (1.0 + A1_mu)

A1_m = A1_rho * dx**3
A1_h = np.full(numPart, dx)
A1_u = np.full(numPart, u_init)
A1_mat = np.full(numPart, matID)
A1_ids = np.arange(1, numPart + 1)
A2_vel = np.zeros((numPart, 3))

c_s = np.sqrt(A / rho_0)
print(f"{numPart} particles ({N_x} x {N_y} x {N_z}), spacing {dx:.4g} m")
print(f"Free-fall time sqrt(H / g) = {np.sqrt(column_height / g):.4g} s")
print(f"Mach number 2 sqrt(g H) / c_s = {2 * np.sqrt(g * column_height) / c_s:.3g}")

# Write ICs to file
with h5py.File(fileOutputName, "w") as f:
    # Header
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [tank_length, tank_depth, tank_height]
    grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFileOutputsPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["Dimension"] = 3

    # Units
    grp = f.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 100.0
    grp.attrs["Unit mass in cgs (U_M)"] = 1000.0
    grp.attrs["Unit time in cgs (U_t)"] = 1.0
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    # Particle group
    grp = f.create_group("/PartType0")
    grp.create_dataset("Coordinates", data=A2_coords, dtype="d")
    grp.create_dataset("Velocities", data=A2_vel, dtype="f")
    grp.create_dataset("Masses", data=A1_m.reshape((numPart, 1)), dtype="f")
    grp.create_dataset("Density", data=A1_rho.reshape((numPart, 1)), dtype="f")
    grp.create_dataset("SmoothingLength", data=A1_h.reshape((numPart, 1)), dtype="f")
    grp.create_dataset("InternalEnergy", data=A1_u.reshape((numPart, 1)), dtype="f")
    grp.create_dataset("ParticleIDs", data=A1_ids.reshape((numPart, 1)), dtype="L")
    grp.create_dataset("MaterialIDs", data=A1_mat.reshape((numPart, 1)), dtype="i")
