################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
################################################################################

import h5py
import numpy as np

# Generates a swift IC file for the OrzagTang vortex in a periodic box

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
B0 = 1.0 / (np.sqrt(4.0 * np.pi))  # Bfield
P0 = gamma * pow(B0, 2)  # Pressure
rho0 = gamma * P0  # Density

fileOutputName = "OrszagTangVortex.hdf5"

# ---------------------------------------------------

# glass = h5py.File("glassCube_64.hdf5", 'r')
# glass = h5py.File("glassCube_32.hdf5", 'r')
glass = h5py.File("glassCube_16.hdf5", "r")
# glass = h5py.File("glassCube_8.hdf5", 'r')
pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

Nold = len(h)
times = 32
############ NEW
cx = times
cy = times
cz = 1

lx = 1.0
ly = 1.0
lz = 1.0 / float(times)  # becasue 2D

pnew = np.zeros((int(Nold * cx * cy * cz), 3))
hnew = np.zeros(int(Nold * cx * cy * cz))
N = Nold * cx * cy * cz

k = 0
c0 = 0
c1 = Nold
for i in range(0, cx):
    for j in range(0, cy):
        for k in range(0, cz):
            print(i, j, k)
            pnew[c0:c1, 0] = pos[:, 0] + i
            pnew[c0:c1, 1] = pos[:, 1] + j
            pnew[c0:c1, 2] = pos[:, 2] + k
            hnew[c0:c1] = h[:]
            c0 = c0 + Nold
            c1 = c1 + Nold

print(len(pnew[:, 1]), " / ", N)
pos = pnew
pnew = 0
pos[:, 0] = pos[:, 0] * lx / cx  # -0.5
pos[:, 1] = pos[:, 1] * ly / cy  # -0.5
pos[:, 2] = pos[:, 2] * lz / cz
h = hnew / cx
hnew = 0
vol = 1.0
vol = lx * ly * lz
# Generate extra arrays
v = np.zeros((N, 3))
b = np.zeros((N, 3))
vp = np.zeros((N, 3))
epa = np.zeros(N)
epb = np.zeros(N)
ids = np.linspace(1, N, N)
m = np.ones(N) * rho0 * vol / N
u = np.ones(N) * P0 / (rho0 * (gamma - 1.0))


v[:, 0] = -np.sin(2.0 * np.pi * pos[:, 1])
v[:, 1] = np.sin(2.0 * np.pi * pos[:, 0])
v[:, 2] = 0.0
b[:, 0] = -B0 * np.sin(2.0 * np.pi * pos[:, 1])
b[:, 1] = B0 * np.sin(4.0 * np.pi * pos[:, 0])
b[:, 2] = 0.0
# epa[:]  = B0*(np.cos(2.*np.pi*pos[:,1])/(2.*np.pi) + np.cos(4.*np.pi*pos[:,0])/(4.*np.pi))
# epb[:]  = pos[:,2]

vp[:, 0] = 0.0
vp[:, 1] = 0.0
vp[:, 2] = B0 * (
    np.cos(2.0 * np.pi * pos[:, 1]) / (2.0 * np.pi)
    + np.cos(4.0 * np.pi * pos[:, 0]) / (4.0 * np.pi)
)

pos[:, 0] = pos[:, 0]  # +0.5
pos[:, 1] = pos[:, 1]  # +0.5

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [lx, ly, lz]  #####
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
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
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("MagneticFluxDensity", data=b, dtype="f")
grp.create_dataset("MagneticVectorPotential", data=vp, dtype="f")

fileOutput.close()
