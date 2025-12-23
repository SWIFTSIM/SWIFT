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

# Generates a swift IC file for the Balsaro-Kim test in a periodic cubic box

# some constants
gamma = 5.0 / 3.0  # Gas adiabatic index
mp = 1.67e-24  # g
pc = 3.086e18  # cm

# unit system
uL = 1e3 * pc   # cm, 1 kpc
uM = 6.757e41   # g, 3.398×10^8 M_sun
ut = 8.071e14   # s, 2.5577×10^7  years
uv = 38.23e5    # cm / s, 38.23 km/s
uU = uv**2      # (cm / s)^2
uA = 1e10

# paramters
rho0 = 1        # Background density 
P0   = 0.3      # Baground pressure 
L    = 0.2      # boxsize

B0  = 3.64731873              # microgauss
B0 /= 1e7 * uM / (ut**2 * uA) # code units


fileName = "BalsaraKim.hdf5"

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

print('smoothing lengths (all,min,max):', h, h.min(), h.max())

numPart = size(h)
vol = L ** 3

# Generate extra arrays
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
u = zeros(numPart)
r = zeros(numPart)
B = zeros((numPart, 3))
A = zeros((numPart, 3))

m[:] = rho0 * vol / numPart
u[:] = P0 / ((gamma - 1.0) * rho0)
B[:,0] = B0


print('particle mass : %.3f' % ( m[0]))
print(u.mean())

A[:,1] = -B0 * (pos[:,2] - 0.5 * L) / 2
A[:,2] =  B0 * (pos[:,1] - 0.5 * L) / 2

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
grp.attrs["Unit current in cgs (U_I)"] = uA
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")

file.close()
