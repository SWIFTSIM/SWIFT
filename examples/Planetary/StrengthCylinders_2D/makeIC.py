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
L = 200  # Number of particles on the side
gamma = 5.0 / 3.0  # Gas adiabatic index
rho = 1.0  # Gas central density
P = 2.5  # Gas central pressure
vx = 0.059
vy = 0.0
fileOutputName = "cylinders.hdf5"
# ---------------------------------------------------

boxsize = 20
vol = boxsize**2

numPart_grid = L * L

L_grid = int(sqrt(numPart_grid))

pos_grid = zeros((numPart_grid, 3))
for i in range(L_grid):
    for j in range(L_grid):
        index = i * L_grid + j
        pos_grid[index, 0] = (i / (float(L_grid))) * boxsize
        pos_grid[index, 1] = (j / (float(L_grid))) * boxsize
h_grid = ones(numPart_grid) * (1.0 / L_grid) * 1.2348
m_grid = ones(numPart_grid) * vol * rho / numPart_grid
u_grid = ones(numPart_grid) * P / (rho * (gamma - 1.0))


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
u = u_grid[select_both]
m = m_grid[select_both]
vel = vel[select_both, :]


numPart = size(h)
ids = linspace(1, numPart, numPart)
mat = zeros(numPart)
rho = ones(numPart)

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [boxsize, boxsize, 0.2]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 2

# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1/85200.0
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
ds = grp.create_dataset("Densities", (numPart, 1), "f")
ds[()] = rho.reshape((numPart, 1))
ds = grp.create_dataset("SmoothingLengths", (numPart, 1), "f")
ds[()] = h.reshape((numPart, 1))
ds = grp.create_dataset("InternalEnergies", (numPart, 1), "f")
ds[()] = u.reshape((numPart, 1))
ds = grp.create_dataset("ParticleIDs", (numPart, 1), "L")
ds[()] = ids.reshape((numPart, 1))
ds = grp.create_dataset("MaterialIDs", (numPart, 1), "i")
ds[()] = mat.reshape((numPart, 1))

fileOutput.close()
