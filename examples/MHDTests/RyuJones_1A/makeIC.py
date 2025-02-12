###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
x_min = -4.0
x_max = -x_min
rho = 1.0
vx_L = 10.0  # Velocity left state
vx_R = -10.0  # Velocity right state
P_L = 20.0  # Pressure left state
P_R = 1.0  # Pressure right state

Bx = 0.0  # 5.0/np.sqrt(4*np.pi)
By = 0.0  # 5.0/np.sqrt(4*np.pi)

# File to save ICs to
fileName = "RyuJones_1A.hdf5"

# Read in number of BCC cubes to stack on either side of the discontinuity
parser = ap.ArgumentParser(
    description="Stack given number of BCC unit cells to generate Brio & Wu shock tube"
)

parser.add_argument(
    "-n",
    "--numcubes",
    help="Number of BCC unit cells to stack on either side. Default: 80",
    default=80,
    type=int,
)

args = parser.parse_args()

# Simulation box attributes
boxSize = x_max - x_min
scale = boxSize / (2 * args.numcubes)
vol = boxSize * scale * scale

glass = h5py.File("BCCglassCube_32.hdf5", "r")

pos_unit = glass["/PartType0/Coordinates"][:, :] * scale
h_unit = glass["/PartType0/SmoothingLength"][:] * scale

pos_L, pos_R = pos_unit, pos_unit
h_L, h_R = h_unit, h_unit

for ii in range(1, args.numcubes):
    pos_L, pos_R = (
        np.vstack((pos_L, pos_unit + np.array([ii * scale, 0.0, 0.0]))),
        np.vstack((pos_R, pos_unit + np.array([ii * scale, 0.0, 0.0]))),
    )
    h_L, h_R = np.vstack((h_L, h_unit)), np.vstack((h_R, h_unit))

pos = np.append(pos_L + np.array([x_min, 0.0, 0.0]), pos_R, axis=0)
h = np.append(h_L, h_R)

numPart_L = np.size(h_L)
numPart_R = np.size(h_R)
numPart = np.size(h)

# Generate extra arrays
v = np.zeros((numPart, 3))
B = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)

m[:] = rho * vol / numPart

# Instantiate extra arrays
B[:, 0] = Bx
B[:, 1] = By

for i in range(numPart):
    x = pos[i, 0]

    if x < 0:  # Left State
        u[i] = P_L / (rho * (gamma - 1.0))
        v[i, 0] = vx_L
    else:  # Right State
        u[i] = P_R / (rho * (gamma - 1.0))
        v[i, 0] = vx_R

    # Generate additional Riemann Problem around periodic boundary to make test runnable
    if x < -3.0:
        v[i, 0] = vx_L * (x - x_min)
    elif x > 3.0:
        v[i, 0] = -vx_R * (x - x_max)

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
