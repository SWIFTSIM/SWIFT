#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt

# Parameters

r_in = 0.1
rho_in_0 = 10.0
rho_out_0 = 1.0
P_0 = 1.0
B_0 = 2.5 / np.sqrt(np.pi)
omega_0 = 20.0
gamma = 1.4

fileInputName = "FastRotor_temp_0000.hdf5"
fileOutputName = "FastRotor.hdf5"

###---------------------------###

particles = h5py.File(fileInputName, "r")

[lx, ly, lz] = particles["/Header"].attrs.get("BoxSize")
N = particles["/Header"].attrs.get("NumPart_Total")[0]
print(N)
N = int(N)

rho = particles["/PartType0/Densities"][:]
pos = particles["/PartType0/Coordinates"][:, :]
v = particles["/PartType0/Velocities"][:, :]
B = particles["/PartType0/MagneticFluxDensities"][:, :]
ids = particles["/PartType0/ParticleIDs"][:]
m = particles["/PartType0/Masses"][:]
u = P_0 / (rho * (gamma - 1))
h = particles["/PartType0/SmoothingLengths"][:]


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
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
# grp.create_dataset("VecPot", data = vp, dtype = 'f')
# grp.create_dataset("EPalpha", data = epa, dtype = 'f')
# grp.create_dataset("EPbeta" , data = epb, dtype = 'f')

fileOutput.close()
