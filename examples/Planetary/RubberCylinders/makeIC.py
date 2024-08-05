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
from numpy import *

# Generates a swift IC file for

# Parameters
L = 200             # Number of particles on the side
rho = 1.0           # Gas central density
c_s = 85200         # Sound speed
vx = 0.059 * c_s
vy = 0.0
fileOutputName = "cylinders.hdf5"
L_depth = 15
# ---------------------------------------------------

boxsize = 20

numPart_grid = L * L * L_depth

pos_grid = zeros((numPart_grid, 3))
for i in range(L):
    for j in range(L):
    	for k in range(L_depth):
    		index = i * L * L_depth + j * L_depth + k
	    	pos_grid[index, 0] = (i / (float(L))) * boxsize
	    	pos_grid[index, 1] = (j / (float(L))) * boxsize
	    	pos_grid[index, 2] = (k / (float(L))) * boxsize
h_grid = ones(numPart_grid) * (boxsize / L) * 1.2348
m_grid = ones(numPart_grid) * rho * (boxsize / L)**3


# Remove the particles not in cylinder
# See Monaghan 2000
initial_offset = 50 * boxsize / L 
outer_r = 4
inner_r = 3

select_cylinderA = logical_and(
    (pos_grid[:, 0] - initial_offset - boxsize/2)**2 + (pos_grid[:, 1] - boxsize/2)**2 <= outer_r**2,
    (pos_grid[:, 0] - initial_offset - boxsize/2)**2 + (pos_grid[:, 1] - boxsize/2)**2 >= inner_r**2,
)


select_cylinderB = logical_and(
    (pos_grid[:, 0] + initial_offset - boxsize/2)**2 + (pos_grid[:, 1] - boxsize/2)**2 <= outer_r**2,
    (pos_grid[:, 0] + initial_offset - boxsize/2)**2 + (pos_grid[:, 1] - boxsize/2)**2 >= inner_r**2,
)

vel = zeros((numPart_grid, 3))
vel[select_cylinderA, 0] = -vx
vel[select_cylinderB, 0] = vx
vel[:, 1] = 0

select_both = logical_or(select_cylinderA, select_cylinderB)

pos = pos_grid[select_both, :]
h = h_grid[select_both]
m = m_grid[select_both]
vel = vel[select_both, :]


numPart = size(h)
ids = linspace(1, numPart, numPart)
mat = 500 * ones(numPart)
rho = ones(numPart)
u = zeros(numPart)

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [boxsize, boxsize, boxsize * L_depth / L]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
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
ds = grp.create_dataset("Coordinates", (numPart, 3), "d")
ds[()] = pos
ds = grp.create_dataset("Velocities", (numPart, 3), "f")
ds[()] = vel
ds = grp.create_dataset("Masses", (numPart, 1), "f")
ds[()] = m.reshape((numPart, 1))
ds = grp.create_dataset("Density", (numPart, 1), "f")
ds[()] = rho.reshape((numPart, 1))
ds = grp.create_dataset("SmoothingLength", (numPart, 1), "f")
ds[()] = h.reshape((numPart, 1))
ds = grp.create_dataset("InternalEnergy", (numPart, 1), "f")
ds[()] = u.reshape((numPart, 1))
ds = grp.create_dataset("ParticleIDs", (numPart, 1), "L")
ds[()] = ids.reshape((numPart, 1))
ds = grp.create_dataset("MaterialIDs", (numPart, 1), "i")
ds[()] = mat.reshape((numPart, 1))

fileOutput.close()
