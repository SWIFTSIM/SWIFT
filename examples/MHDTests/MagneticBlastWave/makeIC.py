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

# Generates a swift IC file for the Magnetic Blast Wave test in a periodic box

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
gamma = 5.0 / 3.0
R_0 = 0.1
rho_0 = 1.0
P_in_0 = 10.0
P_out_0 = 0.1
B_0 = 1 / np.sqrt(2)

fileOutputName = "MagneticBlastWave.hdf5"

# Stack unit cells
glass = h5py.File("glassCube_32.hdf5", "r")
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
ids = np.linspace(1, N, N)
m = np.ones(N) * rho_0 * vol / N
u = np.ones(N)

# Instantiate extra arrays
B[:, 0] = B_0
B[:, 1] = B_0
rot = np.sqrt((pos[:, 0] - 0.5 * lx) ** 2 + (pos[:, 1] - 0.5 * ly) ** 2)
u[rot < R_0] *= P_in_0 / (rho_0 * (gamma - 1))
u[rot >= R_0] = P_out_0 / (rho_0 * (gamma - 1))

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

fileOutput.close()
