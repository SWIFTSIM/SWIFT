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

# Generates a swift IC file for the Sedov blast test in a periodic cubic box

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
n0 = 0.1  # cm^-3
mp = 1.67e-24  # g
kB = 1.381e-16  # erg K^-1
pc = 3.086e18  # cm
rho0 = n0 * mp  # Background density
T0 = 1.0e4  # Background pressure
E0 = 1.0e51  # Energy of the explosion
N_inject = 32  # Number of particles in which to inject energy
L = 256.0 * pc
fileName = "sedov.hdf5"

print(rho0, "cm s^-3")

uL = 1.0e3 * pc
uM = 1.99e30
uv = 1.0e10
ut = uL / uv
uU = uv ** 2
print("ut:", ut / 3.154e7, "yr")

# ---------------------------------------------------
glass = h5py.File("glassCube_64.hdf5", "r")

# Read particle positions and h from the glass
pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

pos *= L
h *= L

for i in range(3):
    pos[pos[:, i] < 0.0, i] += L
    pos[pos[:, i] >= L, i] -= L

numPart = size(h)

# newpos = zeros((numPart*8,3))
# newh = zeros(numPart*8)
# for ix in range(2):
#  for iy in range(2):
#    for iz in range(2):
#      offset = ix*4+iy*2+iz
#      newpos[offset*numPart:(offset+1)*numPart,:] = 0.5*pos[:,:]
#      newpos[offset*numPart:(offset+1)*numPart,0] += 0.5*ix*L
#      newpos[offset*numPart:(offset+1)*numPart,1] += 0.5*iy*L
#      newpos[offset*numPart:(offset+1)*numPart,2] += 0.5*iz*L
#      newh[offset*numPart:(offset+1)*numPart] = 0.5*h[:]

# pos = newpos
# h = newh
print(h, h.min(), h.max())

numPart = size(h)
vol = L ** 3

# Generate extra arrays
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
u = zeros(numPart)
r = zeros(numPart)

r = sqrt(
    (pos[:, 0] - 0.5 * L) ** 2 + (pos[:, 1] - 0.5 * L) ** 2 + (pos[:, 2] - 0.5 * L) ** 2
)
m[:] = rho0 * vol / numPart
# u[:] = P0 / (rho0 * (gamma - 1))
u[:] = kB * T0 / ((gamma - 1.0) * mp)

print(u.mean(), E0 / (N_inject * m[0]))
print(E0 / (N_inject * m[0]) / u.mean())

# Make the central particle detonate
index = argsort(r)
u[index[0:N_inject]] = E0 / (N_inject * m[0])

pos /= uL
h /= uL
L /= uL
m /= uM
u /= uU

print("L:", L)
print("m:", m.mean())
print("h:", h.mean())
print("u:", u.mean(), u.min(), u.max())
print("pos:", pos.min(), pos.max())

# --------------------------------------------------

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]
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
grp.attrs["Unit length in cgs (U_L)"] = uL
grp.attrs["Unit mass in cgs (U_M)"] = uM
grp.attrs["Unit time in cgs (U_t)"] = ut
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")

file.close()
