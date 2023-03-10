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

# Generates a swift IC file for the BrioWu in a periodic box
times = 12  # Number pf Cubes smashed in each side

# Parameters
gamma = 2.0  # Gas adiabatic index
x_min = -1.0
x_max =  1.0
rho_L = 1.0  # Density left state
rho_R = 0.125  # Density right state
v_L = 0.0  # Velocity left state
v_R = 0.0  # Velocity right state
P_L = 1.0  # Pressure left state
P_R = 0.1  # Pressure right state
fileName = "BrioWu_HCP.hdf5"
#fileName = "BrioWu.hdf5"


# ---------------------------------------------------
boxSide = x_max - x_min

#glass_L = h5py.File("glassCube_64.hdf5", "r")
#glass_R = h5py.File("glassCube_32.hdf5", "r")
glass_L = h5py.File("glassCube_64.hdf5", "r")
glass_R = h5py.File("glassCube_32.hdf5", "r")

glass_L = h5py.File("BCC_64.hdf5", "r")
glass_R = h5py.File("BCC_32.hdf5", "r")

pos_L = glass_L["/PartType0/Coordinates"][:, :]
pos_R = glass_R["/PartType0/Coordinates"][:, :]
h_L = glass_L["/PartType0/SmoothingLength"][:]
h_R = glass_R["/PartType0/SmoothingLength"][:]
print(size(h_L), size(h_R))
pos_LL = append(pos_L, pos_L + array([1.0, 0.0, 0.0]), axis=0)
pos_RR = append(pos_R, pos_R + array([1.0, 0.0, 0.0]), axis=0)
h_LL = append(h_L, h_L)
h_RR = append(h_R, h_R)

for i in range(2, times):
    pos_LL = append(pos_LL, pos_L + array([float(i), 0.0, 0.0]), axis=0)
    pos_RR = append(pos_RR, pos_R + array([float(i), 0.0, 0.0]), axis=0)
    h_LL = append(h_LL, h_L)
    h_RR = append(h_RR, h_R)


pos = append(pos_LL - array([times, 0.0, 0.0]), pos_RR, axis=0) / times
h = append(h_LL, h_RR) / times

numPart_L = size(h_LL)
numPart_R = size(h_RR)
numPart = size(h)

#vol_L = 1.0 * 1.0 * boxSide / 2.0
#vol_R = 1.0 * 1.0 * boxSide / 2.0
vol_L = 1.0 * 1.0 * boxSide / 2.0 / times / times
vol_R = 1.0 * 1.0 * boxSide / 2.0 / times / times

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

    if x < 0:  # left
        u[i] = P_L / (rho_L * (gamma - 1.0))
        m[i] = rho_L * vol_L / numPart_L
        v[i, 0] = v_L
        b[i, 0] = 0.75
        b[i, 1] = 1.0
        b[i, 2] = 0.0
        epa[i] = pos[i, 2]
        epb[i] = -0.75 * pos[i, 1] + pos[i, 0]
        vp[i, 0] = 0.0
        vp[i, 1] = -0.75 * pos[i, 2]
        vp[i, 2] = -pos[i, 0]
    else:  # right
        u[i] = P_R / (rho_R * (gamma - 1.0))
        m[i] = rho_R * vol_R / numPart_R
        v[i, 0] = v_R
        b[i, 0] = 0.75
        b[i, 1] = -1.0
        b[i, 2] = 0.0
        epa[i] = pos[i, 2]
        epb[i] = -0.75 * pos[i, 1] - pos[i, 0]
        vp[i, 0] = 0.0
        vp[i, 1] = -0.75 * pos[i, 2]
        vp[i, 2] = pos[i, 0]

# Shift particles
pos[:, 0] -= x_min
# b[:,:]  *= sqrt(4.0*3.14159265)
# vp[:,:] *= sqrt(4.0*3.14159265)

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSide, 1.0/times, 1.0/times]
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
