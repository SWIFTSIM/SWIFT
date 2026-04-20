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

# Generates a swift IC file for the Fluid--Solid interaction test in a periodic box

# Parameters
N_l = 128  # Particles along one edge in the low-density region
N_depth = 18  # Particles in z direction in low-density region
rho1 = 1  # Central region density
rho2 = 1  # Outskirts density
v1 = 1000  # Central region velocity
v2 = 0  # Outskirts velocity
boxsize_l = 1  # size of simulation box in x and y dimension
boxsize_depth = boxsize_l * N_depth / N_l  # size of simulation box in z dimension
fileOutputName = "fluid_solid.hdf5"

# ---------------------------------------------------

num_part = N_l * N_l * N_depth

# Set up grid of particles
i = np.arange(N_l)
j = np.arange(N_l)
k = np.arange(N_depth)
ii, jj, kk = np.meshgrid(i, j, k, indexing="ij")
coords = np.empty((num_part, 3))
coords[:, 0] = (jj.ravel() / N_l + 1.0 / (2.0 * N_l)) * boxsize_l
coords[:, 1] = (kk.ravel() / N_l + 1.0 / (2.0 * N_l)) * boxsize_l
coords[:, 2] = (ii.ravel() / N_depth + 1.0 / (2.0 * N_depth)) * boxsize_depth

# Set up other arrays
vel = np.zeros((num_part, 3))
rho = np.empty(num_part)
m = np.empty(num_part)
mat = np.empty(num_part)
h = np.full(num_part, boxsize_l / N_l)
ids = np.arange(1, num_part + 1)
u = np.zeros(num_part)

# region masks
mask_out = (coords[:, 1] >= 0.25) & (coords[:, 1] <= 0.75)
mask_in = ~mask_out

# density and mass
rho[mask_out] = rho1
rho[mask_in] = rho2
m = rho / (N_l ** 3)

# shear velocity
vel[mask_out, 0] = v1
vel[mask_in, 0] = v2

# material IDs
mat[mask_out] = 500
mat[mask_in] = 501

with h5py.File(fileOutputName, "w") as f:

    hdr = f.create_group("/Header")
    hdr.attrs["BoxSize"] = np.array([boxsize_l, boxsize_l, boxsize_depth])
    hdr.attrs["NumPart_Total"] = np.array([num_part, 0, 0, 0, 0, 0])
    hdr.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.int32)
    hdr.attrs["NumPart_ThisFile"] = np.array([num_part, 0, 0, 0, 0, 0])
    hdr.attrs["Time"] = 0.0
    hdr.attrs["NumFilesPerSnapshot"] = 1
    hdr.attrs["MassTable"] = np.zeros(6)
    hdr.attrs["Flag_Entropy_ICs"] = np.zeros(6, dtype=np.int32)
    hdr.attrs["Dimension"] = 3

    units = f.create_group("/Units")
    units.attrs["Unit length in cgs (U_L)"] = 1.0
    units.attrs["Unit mass in cgs (U_M)"] = 1.0
    units.attrs["Unit time in cgs (U_t)"] = 1.0
    units.attrs["Unit current in cgs (U_I)"] = 1.0
    units.attrs["Unit temperature in cgs (U_T)"] = 1.0

    part = f.create_group("/PartType0")

    part.create_dataset("Coordinates", data=coords)
    part.create_dataset("Velocities", data=vel.astype(np.float32))
    part.create_dataset("Masses", data=m.reshape(-1, 1).astype(np.float32))
    part.create_dataset("Density", data=rho.reshape(-1, 1).astype(np.float32))
    part.create_dataset("SmoothingLength", data=h.reshape(-1, 1).astype(np.float32))
    part.create_dataset("InternalEnergy", data=u.reshape(-1, 1).astype(np.float32))
    part.create_dataset("ParticleIDs", data=ids.reshape(-1, 1))
    part.create_dataset("MaterialIDs", data=mat.reshape(-1, 1))