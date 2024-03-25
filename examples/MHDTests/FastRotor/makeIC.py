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

# Generates a swift IC file for the Fast Rotor MHD test in a periodic box

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
gamma = 7.0 / 5.0
R_0 = 0.1
rho_in_0 = 10.0
rho_out_0 = 1.0
v_0 = 2.0
P_0 = 1.0
B_0 = 2.5 / np.sqrt(np.pi)

fileOutputName = "FastRotor.hdf5"

# Retrieve unit cells
glass = h5py.File("glassCube_32.hdf5", "r")
pos_unit_cell = glass["/PartType0/Coordinates"][:, :]
h_unit_cell = glass["/PartType0/SmoothingLength"][:]

N_unit_cell = len(h_unit_cell)

# Stack unit cells to generate ambient medium
times = args.replications

cx_out, cy_out, cz_out = times, times, 1

lx, ly, lz = 1.0, 1.0, 1.0 / times

N_out = N_unit_cell * cx_out * cy_out * cz_out

pos_out = np.zeros((int(N_out), 3))
h_out = np.zeros(N_out)

c0 = 0
c1 = N_unit_cell

for i in range(0, cx_out):
    for j in range(0, cy_out):
        for k in range(0, cz_out):
            pos_out[c0:c1, 0] = pos_unit_cell[:, 0] + i
            pos_out[c0:c1, 1] = pos_unit_cell[:, 1] + j
            pos_out[c0:c1, 2] = pos_unit_cell[:, 2] + k
            h_out[c0:c1] = h_unit_cell[:]
            c0 = c0 + N_unit_cell
            c1 = c1 + N_unit_cell

# Rescale to desired length
pos_out[:, 0] *= lx / cx_out
pos_out[:, 1] *= ly / cy_out
pos_out[:, 2] *= lz / cz_out
h_out *= lx / cx_out

# Stack unit cells to generate dense disk
ratio = np.cbrt(rho_in_0 / rho_out_0)

cx_in, cy_in, cz_in = (
    int(np.ceil(ratio * cx_out)),
    int(np.ceil(ratio * cy_out)),
    int(np.ceil(ratio * cz_out)),
)

lx_in = lx + (np.ceil(ratio * cx_out) / (ratio * cx_out) - 1.0) * lx
ly_in = ly + (np.ceil(ratio * cy_out) / (ratio * cy_out) - 1.0) * ly
lz_in = lz + (np.ceil(ratio * cz_out) / (ratio * cz_out) - 1.0) * lz

N_in = N_unit_cell * cx_in * cy_in * cz_in

pos_in = np.zeros((int(N_in), 3))
h_in = np.zeros(N_in)

c0 = 0
c1 = N_unit_cell

for i in range(0, cx_in):
    for j in range(0, cy_in):
        for k in range(0, cz_in):
            pos_in[c0:c1, 0] = pos_unit_cell[:, 0] + i
            pos_in[c0:c1, 1] = pos_unit_cell[:, 1] + j
            pos_in[c0:c1, 2] = pos_unit_cell[:, 2] + k
            h_in[c0:c1] = h_unit_cell[:]
            c0 = c0 + N_unit_cell
            c1 = c1 + N_unit_cell

# Rescale to desired length
pos_in[:, 0] = pos_in[:, 0] * lx_in / cx_in
pos_in[:, 1] = pos_in[:, 1] * ly_in / cy_in
pos_in[:, 2] = pos_in[:, 2] * lz_in / cz_in
h_in = h_in / cx_in

# Combine dense disk with medium
vol = lx * ly * lz

rot_out = np.sqrt((pos_out[:, 0] - 0.5 * lx) ** 2 + (pos_out[:, 1] - 0.5 * ly) ** 2)
rot_in = np.sqrt((pos_in[:, 0] - 0.5 * lx) ** 2 + (pos_in[:, 1] - 0.5 * ly) ** 2)

pos_out_f = pos_out[rot_out ** 2 >= R_0 ** 2]
pos_in_f = pos_in[rot_in ** 2 < R_0 ** 2]

h_out_f = h_out[rot_out ** 2 > R_0 ** 2]
h_in_f = h_in[rot_in ** 2 < R_0 ** 2]
h_in_f = h_in_f[pos_in_f[:, 2] < lz]

pos_in_f = pos_in_f[pos_in_f[:, 2] < lz]

pos = np.append(pos_out_f, pos_in_f, axis=0)
h = np.append(h_out_f, h_in_f, axis=0)
N = len(pos)

# Generate extra arrays
v = np.zeros((N, 3))
B = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho_out_0 * vol / N_out
u = np.ones(N)

rot = np.sqrt((pos[:, 0] - 0.5 * lx) ** 2 + (pos[:, 1] - 0.5 * ly) ** 2)

v[:, 0][rot < R_0] = -v_0 * (pos[:, 1][rot < R_0] - 0.5 * ly) / rot[rot < R_0]
v[:, 1][rot < R_0] = v_0 * (pos[:, 0][rot < R_0] - 0.5 * lx) / rot[rot < R_0]

B[:, 0] = B_0

u[rot < R_0] *= P_0 / (rho_in_0 * (gamma - 1))
u[rot >= R_0] *= P_0 / (rho_out_0 * (gamma - 1))

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
