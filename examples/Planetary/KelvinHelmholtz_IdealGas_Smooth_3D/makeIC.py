###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
#               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

import h5py
import numpy as np

# Generates a swift IC file for the Kelvin-Helmholtz test in a periodic box

# Parameters
N_l = 128  # Particles along one edge in the low-density region
N_depth = 18  # Particles in z direction in low-density region
gamma = 5.0 / 3.0  # Gas adiabatic index
P1 = 2.5  # Central region pressure
P2 = 2.5  # Outskirts pressure
rho1 = 2  # Central region density
rho2 = 1  # Outskirts density
v1 = 0.5  # Central region velocity
v2 = -0.5  # Outskirts velocity
boxsize_l = 1  # size of simulation box in x and y dimension
boxsize_depth = boxsize_l * N_depth / N_l  # size of simulation box in z dimension
fileOutputName = "kelvin_helmholtz.hdf5"

# Parameters for smoothing of interfaces
vm = (v2 - v1) / 2
rhom = (rho2 - rho1) / 2
delta = 0.025
# ---------------------------------------------------

numPart = N_l * N_l * N_depth

# Now construct two lattices of particles in the two regions
A2_coords = np.empty((numPart, 3))
A2_vel = np.zeros((numPart, 3))

A1_mat = np.zeros(numPart)
A1_m = np.empty(numPart)
A1_rho = np.empty(numPart)
A1_u = np.empty(numPart)
A1_h = np.ones(numPart) * boxsize_l / N_l
A1_ids = np.linspace(1, numPart, numPart)

# Set up particles
for i in range(N_depth):
    for j in range(N_l):
        for k in range(N_l):
            index = i * N_l ** 2 + j * N_l + k

            x = (j / float(N_l) + 1.0 / (2.0 * N_l)) * boxsize_l
            y = (k / float(N_l) + 1.0 / (2.0 * N_l)) * boxsize_l
            z = (i / float(N_depth) + 1.0 / (2.0 * N_depth)) * boxsize_depth

            A2_coords[index, 0] = x
            A2_coords[index, 1] = y
            A2_coords[index, 2] = z

            if 0.0 <= y <= 0.25:
                A1_rho[index] = rho2 - rhom * np.exp((y - 0.25) / delta)
                A2_vel[index, 0] = v2 - vm * np.exp((y - 0.25) / delta)
                A1_m[index] = A1_rho[index] / N_l ** 3
                A1_u[index] = P2 / (A1_rho[index] * (gamma - 1.0))

            elif 0.25 <= y <= 0.5:
                A1_rho[index] = rho1 + rhom * np.exp((0.25 - y) / delta)
                A2_vel[index, 0] = v1 + vm * np.exp((0.25 - y) / delta)
                A1_m[index] = A1_rho[index] / N_l ** 3
                A1_u[index] = P1 / (A1_rho[index] * (gamma - 1.0))

            elif 0.5 <= y <= 0.75:
                A1_rho[index] = rho1 + rhom * np.exp((y - 0.75) / delta)
                A2_vel[index, 0] = v1 + vm * np.exp((y - 0.75) / delta)
                A1_m[index] = A1_rho[index] / N_l ** 3
                A1_u[index] = P1 / (A1_rho[index] * (gamma - 1.0))

            elif 0.75 <= y <= 1:
                A1_rho[index] = rho2 - rhom * np.exp((0.75 - y) / delta)
                A2_vel[index, 0] = v2 - vm * np.exp((0.75 - y) / delta)
                A1_m[index] = A1_rho[index] / N_l ** 3
                A1_u[index] = P2 / (A1_rho[index] * (gamma - 1.0))

# Finally add the velocity perturbation
vel_perturb_factor = 0.01 * (v1 - v2)
A2_vel[:, 1] = vel_perturb_factor * np.sin(
    2 * np.pi * A2_coords[:, 0] / (0.5 * boxsize_l)
)

# Write ICs to file
with h5py.File(fileOutputName, "w") as f:
    # Header
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [boxsize_l, boxsize_l, boxsize_depth]
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
    grp.attrs["Unit length in cgs (U_L)"] = 1.0
    grp.attrs["Unit mass in cgs (U_M)"] = 1.0
    grp.attrs["Unit time in cgs (U_t)"] = 1.0
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    # Particle group
    grp = f.create_group("/PartType0")
    ds = grp.create_dataset("Coordinates", (numPart, 3), "d")
    ds[()] = A2_coords
    ds = grp.create_dataset("Velocities", (numPart, 3), "f")
    ds[()] = A2_vel
    ds = grp.create_dataset("Masses", (numPart, 1), "f")
    ds[()] = A1_m.reshape((numPart, 1))
    ds = grp.create_dataset("Density", (numPart, 1), "f")
    ds[()] = A1_rho.reshape((numPart, 1))
    ds = grp.create_dataset("SmoothingLength", (numPart, 1), "f")
    ds[()] = A1_h.reshape((numPart, 1))
    ds = grp.create_dataset("InternalEnergy", (numPart, 1), "f")
    ds[()] = A1_u.reshape((numPart, 1))
    ds = grp.create_dataset("ParticleIDs", (numPart, 1), "L")
    ds[()] = A1_ids.reshape((numPart, 1))
    ds = grp.create_dataset("MaterialIDs", (numPart, 1), "i")
    ds[()] = A1_mat.reshape((numPart, 1))
