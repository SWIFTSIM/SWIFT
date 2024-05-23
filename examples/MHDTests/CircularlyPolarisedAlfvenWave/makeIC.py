#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt

# Parameters

rho = 1.0
P = 0.1
gamma = 5.0 / 3.0

B0 = 0.1
wavelen = 1.0
wavenum = 2.0 * np.pi / wavelen

fileOutputName = "CircularlyPolarisedAlfvenWave.hdf5"

###---------------------------###

glass = h5py.File("FCCglassCube_32.hdf5", "r")

unit_cell = glass["/PartType0/Coordinates"][:, :]
h_unit_cell = glass["/PartType0/SmoothingLength"][:]

N_unit_cell = len(h_unit_cell)
times = 2

###---------------------------###

cx = times
cy = 1
cz = 1

Lx = 3.0
Ly = 0.5 * Lx
Lz = 0.5 * Lx

N = int(N_unit_cell * cx * cy * cz)

pos = np.zeros((N, 3))
h = np.zeros(N)

c0 = 0
c1 = N_unit_cell

for i in range(0, cx):
    for j in range(0, cy):
        for k in range(0, cz):
            pos[c0:c1, 0] = unit_cell[:, 0] + i
            pos[c0:c1, 1] = unit_cell[:, 1] + j
            pos[c0:c1, 2] = unit_cell[:, 2] + k
            h[c0:c1] = h_unit_cell[:]
            c0 = c0 + N_unit_cell
            c1 = c1 + N_unit_cell

pos[:, 0] = pos[:, 0] * Lx / cx
pos[:, 1] = pos[:, 1] * Ly / cy
pos[:, 2] = pos[:, 2] * Lz / cz
h = h * Lx / cx

vol = Lx * Ly * Lz

###---------------------------###

sina = 2.0 / 3.0
cosa = np.sqrt(1 - sina ** 2)

sinb = 2.0 / np.sqrt(5)
cosb = np.sqrt(1 - sinb ** 2)

Rotation = np.array(
    [
        [cosa * cosb, -sinb, -sina * cosb],
        [cosa * sinb, cosb, -sina * sinb],
        [sina, 0, cosa],
    ]
)

# Rotationinv = np.linalg.inv(Rotation)

# x = np.matmul(Rotationinv, pos.T)
x1 = (pos[:, 0] + 2 * pos[:, 1] + 2 * pos[:, 2]) / 3
v = np.zeros((3, N))
B = np.zeros((3, N))
A = np.zeros((3, N))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho * vol / N
u = np.ones(N) * P / (rho * (gamma - 1))

v[0, :] = 0.0
v[1, :] = 0.1 * np.sin(wavenum * x1)
v[2, :] = 0.1 * np.cos(wavenum * x1)

B[0, :] = 1.0
B[1, :] = B0 * np.sin(wavenum * x1)
B[2, :] = B0 * np.cos(wavenum * x1)

A[1, :] = B0 / wavenum * np.sin(wavenum * x1)
A[2, :] = B0 / wavenum * np.cos(wavenum * x1)

v = np.matmul(Rotation, v)
B = np.matmul(Rotation, B)
A = np.matmul(Rotation, A)

v = v.T
B = B.T
A = A.T

###---------------------------###

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [Lx, Ly, Lz]
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
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")

fileOutput.close()
