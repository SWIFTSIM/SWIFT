#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt

# Parameters

r_0 = 1.0 / np.sqrt(8)
rho_0 = 1.0
P_0 = 6.0
Bx_0 = 1.0 / np.sqrt(4 * np.pi)
Bz_0 = 1.0 / np.sqrt(4 * np.pi)
gamma = 5.0 / 3.0

fileOutputName = "Monopole.hdf5"

###---------------------------###

glass = h5py.File("FCCglassCube_48.hdf5", "r")

pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

N = len(h)

###---------------------------###

lx = 2.0
ly = 2.0
lz = 2.0

lx_c = r_0  # lx/2
ly_c = r_0  # ly/2
lz_c = lz / 2

pos[:, 0] = pos[:, 0] * lx
pos[:, 1] = pos[:, 1] * ly
pos[:, 2] = pos[:, 2] * lz
h = h * lx

vol = lx * ly * lz

###---------------------------###

rot = np.sqrt(
    (pos[:, 0] - lx_c) ** 2 + (pos[:, 1] - ly_c) ** 2
)  # + (pos[:,2]-lz_c)**2)
v = np.zeros((N, 3))
B = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho_0 * vol / N
u = np.ones(N) * P_0 / (rho_0 * (gamma - 1))

v[:, 0] = 1.0
v[:, 1] = 1.0

B[:, 0][rot < r_0] = Bx_0 * (
    (rot[rot < r_0] / r_0) ** 8 - 2.0 * (rot[rot < r_0] / r_0) ** 4 + 1.0
)
B[:, 2] = Bz_0

###---------------------------###

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
grp.create_dataset("MagneticFluxDensity", data=B, dtype="f")
# grp.create_dataset("VecPot", data = vp, dtype = 'f')
# grp.create_dataset("EPalpha", data = epa, dtype = 'f')
# grp.create_dataset("EPbeta" , data = epb, dtype = 'f')

fileOutput.close()
