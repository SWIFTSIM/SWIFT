################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#               2021 Federico Stasyszyn (fstasyszyn@unc.edu.ar)
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

import h5py
import argparse
import numpy as np

# Generates a swift IC file for the Orszag-Tang vortex in a periodic box

# Read in number of BCC cells to stack along x and y
argparser = argparse.ArgumentParser()
argparser.add_argument(
    "-r",
    "--replications",
    help="Number of BCC unit cells to stack along the x and y directions. Default: 4",
    default=4,
    type=int,
)
args = argparser.parse_args()

# Parameters
gamma = 5.0 / 3.0 # Gas adiabatic index
v0 = 1.0          # Velocity
B0 = 1.0 / (np.sqrt(4.0 * np.pi)) # Magnetic field
P0 = gamma * B0 * B0    # Pressure
rho0 = gamma * P0 / v0  # Density

fileOutputName = "OrszagTangVortex.hdf5"

# Stack unit cells
glass = h5py.File("BCCglassCube_32.hdf5", "r")
pos_unit_cell = glass["/PartType0/Coordinates"][:, :]
h_unit_cell = glass["/PartType0/SmoothingLength"][:]

N_unit_cell = len(h_unit_cell)

times = args.replications

cx, cy, cz = times, times, 1

lx, ly, lz = 1.0, 1.0, 1.0 / times

pos = np.zeros((int(N_unit_cell * cx * cy * cz), 3))
h = np.zeros(int(N_unit_cell * cx * cy * cz))
N = N_unit_cell * cx * cy * cz

c0 = 0
c1 = N_unit_cell
for i in range(cx):
    for j in range(cy):
        for k in range(cz):
            pos[c0:c1, 0] = pos_unit_cell[:, 0] + i
            pos[c0:c1, 1] = pos_unit_cell[:, 1] + j
            pos[c0:c1, 2] = pos_unit_cell[:, 2] + k
            h[c0:c1] = h_unit_cell[:]
            c0 += N_unit_cell
            c1 += N_unit_cell

# Rescale to desired length
pos *= lx / cx
h *= lx / cx
vol = lx * ly * lz

# Generate extra arrays
v = np.zeros((N, 3))
B = np.zeros((N, 3))
A = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho0 * vol / N
u = np.ones(N) * P0 / (rho0 * (gamma - 1.0))

# Instantiate extra arrays
v[:, 0] = - v0 * np.sin(2.0 * np.pi * pos[:, 1])
v[:, 1] = v0 * np.sin(2.0 * np.pi * pos[:, 0])

B[:, 0] = -B0 * np.sin(2.0 * np.pi * pos[:, 1])
B[:, 1] = B0 * np.sin(4.0 * np.pi * pos[:, 0])

A[:, 2] = B0 * (
    np.cos(2.0 * np.pi * pos[:, 1]) / (2.0 * np.pi)
    + np.cos(4.0 * np.pi * pos[:, 0]) / (4.0 * np.pi)
)

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [lx, ly, lz]
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")

fileOutput.close()
