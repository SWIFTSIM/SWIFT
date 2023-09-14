###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
##############################################################################

import h5py
from numpy import *
# Parameters
gamma = 2.0  # Gas adiabatic index
x_min = 0.0
x_max = 1.0
P = 1.0
rho = 1.0
fileName = "CurlB.hdf5"


# ---------------------------------------------------
boxSide = x_max - x_min

glass = h5py.File("glassCube_16.hdf5", "r")


pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]
numPart = size(h)

#vol_L = 1.0 * 1.0 * boxSide / 2.0
#vol_R = 1.0 * 1.0 * boxSide / 2.0
vol = boxSide * boxSide * boxSide
# Generate extra arrays
v = zeros((numPart, 3))
b = zeros((numPart, 3))
vp = zeros((numPart, 3))
epa = zeros(numPart)
epb = zeros(numPart)

ids = linspace(1, numPart, numPart)
m = zeros(numPart)
u = zeros(numPart)

for i in range(numPart):
    x = pos[i, 0]
    y = pos[i, 1]
    z = pos[i,2]
    u[i]= P / (rho * (gamma -1.0))
    m[i]= rho * vol / numPart
    v[i,0] = 0.0
    v[i,1] = 0.0
    v[i,2] = 0.0
    b[i,0] = 0.5*(-y-0.5)
    b[i,1] = 0.5*(x+0.5)
    b[i,2] = 0.0
    vp[i, 0] = 0.0
    vp[i, 1] = 0.0
    vp[i, 2] = 0.0
    epa[i]=0.0
    epb[i]=0.0

# Shift particles
pos[:, 0] -= x_min
# b[:,:]  *= sqrt(4.0*3.14159265)
# vp[:,:] *= sqrt(4.0*3.14159265)

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSide, 1.0, 1.0]
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
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
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
grp.create_dataset("MagneticFluxDensity", data=b, dtype="f")
grp.create_dataset("VecPot", data=vp, dtype="f")
grp.create_dataset("EPalpha", data=epa, dtype="f")
grp.create_dataset("EPbeta", data=epb, dtype="f")


file.close()
