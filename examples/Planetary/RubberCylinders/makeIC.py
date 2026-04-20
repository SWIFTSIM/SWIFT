###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Generates a swift IC file for

# Parameters
L = 200             # Number of particles on the side
density = 1.0       # Gas central density
c_s = 85200         # Sound speed
vx = 0.059 * c_s
vy = 0.0
fileOutputName = "cylinders.hdf5"
L_depth = 15
boxsize = 20

# Cylinder parameters
outer_r = 4
inner_r = 3
initial_offset = 50 * boxsize / L 
# ---------------------------------------------------

# Set up grid of particles
numPart_grid = L * L * L_depth
i = np.arange(L)
j = np.arange(L)
k = np.arange(L_depth)
ii, jj, kk = np.meshgrid(i, j, k, indexing="ij")

pos_grid = np.empty((numPart_grid, 3), dtype=np.float64)
pos_grid[:, 0] = (ii.ravel() / L) * boxsize
pos_grid[:, 1] = (jj.ravel() / L) * boxsize
pos_grid[:, 2] = (kk.ravel() / L) * boxsize

h_grid = np.full(numPart_grid, (boxsize / L) * 1.2348)
m_grid = np.full(numPart_grid, density * (boxsize / L)**3)

# Cylinder masks 
select_cylinderA = np.logical_and(
    (pos_grid[:, 0] - initial_offset - boxsize/2)**2 + (pos_grid[:, 1] - boxsize/2)**2 <= outer_r**2,
    (pos_grid[:, 0] - initial_offset - boxsize/2)**2 + (pos_grid[:, 1] - boxsize/2)**2 >= inner_r**2,
)
select_cylinderB = np.logical_and(
    (pos_grid[:, 0] + initial_offset - boxsize/2)**2 + (pos_grid[:, 1] - boxsize/2)**2 <= outer_r**2,
    (pos_grid[:, 0] + initial_offset - boxsize/2)**2 + (pos_grid[:, 1] - boxsize/2)**2 >= inner_r**2,
)

# Set up velocities
vel = np.zeros((numPart_grid, 3))
vel[select_cylinderA, 0] = -vx
vel[select_cylinderB, 0] = vx

# Select particles in the cylinders
select_both = np.logical_or(select_cylinderA, select_cylinderB)
pos = pos_grid[select_both, :]
h = h_grid[select_both]
m = m_grid[select_both]
vel = vel[select_both, :]

# Set up other arrays
numPart = np.size(h)
ids = np.linspace(1, numPart, numPart)
mat = 500 * np.ones(numPart)
rho = np.ones(numPart)
u = np.zeros(numPart)

# File
with h5py.File(fileOutputName, "w") as f:
    hdr = f.create_group("/Header")
    hdr.attrs["BoxSize"] = np.array([boxsize, boxsize, boxsize * L_depth / L])
    hdr.attrs["NumPart_Total"] = np.array([numPart, 0, 0, 0, 0, 0])
    hdr.attrs["NumPart_Total_HighWord"] = np.zeros(6, dtype=np.int32)
    hdr.attrs["NumPart_ThisFile"] = np.array([numPart, 0, 0, 0, 0, 0])
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

    part.create_dataset("Coordinates", data=pos)
    part.create_dataset("Velocities", data=vel.astype(np.float32))
    part.create_dataset("Masses", data=m.reshape(-1, 1).astype(np.float32))
    part.create_dataset("Density", data=rho.reshape(-1, 1).astype(np.float32))
    part.create_dataset("SmoothingLength", data=h.reshape(-1, 1).astype(np.float32))
    part.create_dataset("InternalEnergy", data=u.reshape(-1, 1).astype(np.float32))
    part.create_dataset("ParticleIDs", data=ids.reshape(-1, 1))
    part.create_dataset("MaterialIDs", data=mat.reshape(-1, 1))
