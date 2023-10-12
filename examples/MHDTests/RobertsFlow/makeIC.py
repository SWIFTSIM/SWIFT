#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt

# Parameters

rho = 1.0
cs2 = 3025.0
L = 1.0
k = 2 * np.pi / L
V0 = 16 * np.pi

gamma = 5.0 / 3.0
u0 = cs2 / (gamma * (gamma - 1))

fileOutputName = "RobertsFlow.hdf5"

###---------------------------###

glass = h5py.File("../FCCglassCube_64.hdf5", "r")

pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

N = len(h)
vol = L ** 3

###---------------------------###

v = np.zeros((N, 3))
B = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho * vol / N
u = np.ones(N) * u0

v[:, 0] = V0 * np.sin(k * pos[:, 1]) * np.cos(k * pos[:, 0])
v[:, 1] = -V0 * np.sin(k * pos[:, 0]) * np.cos(k * pos[:, 1])
v[:, 2] = V0 * np.cos(k * pos[:, 1]) * np.cos(k * pos[:, 0]) * np.sqrt(2)
# Why amplitude is k**2?
B[:, 0] = k * k * np.sin(k * pos[:, 1]) * np.sin(k * pos[:, 2])


###---------------------------###

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]  #####
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

fileOutput.close()
