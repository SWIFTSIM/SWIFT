###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
#               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

# Generates a swift IC file for the Square test in a periodic box

# Parameters
L = 2 * 64  # Number of particles on the side
gamma = 5.0 / 3.0  # Gas adiabatic index
rho0 = 4.0  # Gas central density
rho1 = 1.0  # Gas outskirt density
P0 = 2.5  # Gas central pressure
P1 = 2.5  # Gas central pressure
vx = 0.0  # Random velocity for all particles
vy = 0.0
fileOutputName = "square.hdf5"
# ---------------------------------------------------

vol = 1.0

numPart = L * L

pos_x = arange(0, 1, 1.0 / L)
xv, yv = meshgrid(pos_x, pos_x)
pos = zeros((numPart, 3), dtype=float)
pos[:, 0] = xv.flatten()
pos[:, 1] = yv.flatten()

# Now we can get 2d masks!
inside = logical_and.reduce([xv < 0.75, xv > 0.25, yv < 0.75, yv > 0.25])

mass_in = rho0 / numPart
mass_out = rho1 / numPart

m = ones_like(xv) * mass_out
m[inside] = mass_in
m = m.flatten()

h = ones_like(m) / L

u_in = P0 / ((gamma - 1) * rho0)
u_out = P1 / ((gamma - 1) * rho1)
u = ones_like(xv) * u_out
u[inside] = u_in
u = u.flatten()
vel = zeros_like(pos)
ids = arange(numPart)
mat = zeros(numPart)


# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [vol, vol, 0.2]
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
ds = grp.create_dataset("SmoothingLength", (numPart, 1), "f")
ds[()] = h.reshape((numPart, 1))
ds = grp.create_dataset("InternalEnergy", (numPart, 1), "f")
ds[()] = u.reshape((numPart, 1))
ds = grp.create_dataset("ParticleIDs", (numPart, 1), "L")
ds[()] = ids.reshape((numPart, 1))
ds = grp.create_dataset("MaterialIDs", (numPart, 1), "i")
ds[()] = mat.reshape((numPart, 1))


fileOutput.close()
