###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2021 Federico Stasyszyn (fstasyszyn@unc.edu.ar)
#               2023 Orestis Karapiperis (karapiperis@strw.leidenuniv.nl)
#
# This file generates an IC file to run a 3D Brio & Wu shock tube
# (Brio & Wu, 1998 - https://ui.adsabs.harvard.edu/abs/1988JCoPh..75..400B/abstract).
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

import argparse as ap
import h5py
import numpy as np

# Problem Parameters
gamma = 2.0  # Gas adiabatic index
x_min = -1.0
x_max = -x_min
rho_L = 1.0  # Density left state
rho_R = 0.125  # Density right state
P_L = 1.0  # Pressure left state
P_R = 0.1  # Pressure right state
Bx = 0.75  # Magnetic field, x direction
By_L = 1.0  # Magnetic field, y direction, left state
By_R = -1.0  # Magnetic field, y direction, right state

# File to save ICs to
fileName = "BrioWu.hdf5"

# Read in number of BCC cubes to stack on either side of the discontinuity
parser = ap.ArgumentParser(
    description="Stack given number of BCC unit cells to generate Brio & Wu shock tube"
)

parser.add_argument(
    "-n",
    "--numcubes",
    help="Number of BCC unit cells to stack on either side. Default: 11",
    default=11,
    type=int,
)

args = parser.parse_args()

# Simulation box attributes
boxSize = x_max - x_min
scale = boxSize / (
    2 * args.numcubes
)  # Length of y/z side to keep length of tube along x direction to 1
vol = boxSize * scale * scale

# Stack BCC cubes
glass_L = h5py.File("BCCglassCube_48.hdf5", "r")
glass_R = h5py.File("BCCglassCube_24.hdf5", "r")

pos_unit_L = glass_L["/PartType0/Coordinates"][:, :] * scale
h_unit_L = glass_L["/PartType0/SmoothingLength"][:] * scale

pos_unit_R = glass_R["/PartType0/Coordinates"][:, :] * scale
h_unit_R = glass_R["/PartType0/SmoothingLength"][:] * scale

pos_L, pos_R = pos_unit_L, pos_unit_R
h_L, h_R = h_unit_L, h_unit_R

for ii in range(1, args.numcubes):
    pos_L, pos_R = (
        np.vstack((pos_L, pos_unit_L + np.array([ii * scale, 0.0, 0.0]))),
        np.vstack((pos_R, pos_unit_R + np.array([ii * scale, 0.0, 0.0]))),
    )
    h_L, h_R = np.vstack((h_L, h_unit_L)), np.vstack((h_R, h_unit_R))

pos = np.append(pos_L + np.array([x_min, 0.0, 0.0]), pos_R, axis=0)
h = np.append(h_L, h_R)

numPart_L = np.size(h_L)
numPart_R = np.size(h_R)
numPart = np.size(h)

vol_L = boxSize / 2 * scale * scale
vol_R = boxSize / 2 * scale * scale

# Initialise extra arrays
v = np.zeros((numPart, 3))
B = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)

# Instatiate the extra arrays
B[:, 0] = Bx

for i in range(numPart):
    x = pos[i, 0]

    if x < 0:  # Left state
        m[i] = rho_L * vol_L / numPart_L
        u[i] = P_L / (rho_L * (gamma - 1.0))
        B[i, 1] = By_L
    else:  # Right state
        m[i] = rho_R * vol_R / numPart_R
        u[i] = P_R / (rho_R * (gamma - 1.0))
        B[i, 1] = By_R

# Shift particles
pos[:, 0] -= x_min

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize, scale, scale]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")

file.close()
