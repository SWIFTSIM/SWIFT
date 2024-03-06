#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import argparse as ap
import numpy as np
import matplotlib.pyplot as plt

parser = ap.ArgumentParser()

parser.add_argument(
    "-p",
    "--positionsFile",
    help="File with partcile positions to generate Circularly Polarised Alfven Wave Test ICs",
    default="CPLglassCube_32x_16y_16z.hdf5",
    type=str,
)

args = parser.parse_args()

# Parameters
rho = 1.0
P = 0.1
gamma = 5.0 / 3.0

v_par = 0
B_par = 1
v0 = 0.1
B0 = 0.1

theta = 30 # degrees

fileOutputName = "CircularlyPolarisedAlfvenWave.hdf5"

# Retrieve positions file
glass = h5py.File("glassPlane_32.hdf5", "r")
glassBoxSize = glass["/Header"].attrs["BoxSize"]
pos = glass["/PartType0/Coordinates"][:, :]
h   = glass["/PartType0/SmoothingLength"][:]

# Duplicate along y
shifted = pos + np.array([0, 1, 0])
pos = np.append(pos, shifted, axis=0)
h = np.append(h, h)
glassBoxSize[1] *= 2.                

# Get the final number of particles
N = len(h)

# Expected dimensions
Lx = 1. / np.cos(30 * np.pi / 180.)
Ly = 1. / np.sin(30 * np.pi / 180.)

# Scaling factors
scaling = [Lx / glassBoxSize[0], Ly / glassBoxSize[1]]

print("The scaling is : ", scaling)

pos[:,0] *= scaling[0]
pos[:,1] *= scaling[1]
h *= np.sqrt(scaling[0]**2 + scaling[1]**2)

# Info about volume
vol = Lx * Ly 
print("The box size is : ", Lx, Ly)
print("The max (x,y,z) is : ", max(pos[:,0]), max(pos[:,1]), max(pos[:,2]))

# Set the masses, ids and energies
m = np.ones(N) * rho * vol / N
ids = np.linspace(1, N, N)
u = np.ones(N) * P / (rho * (gamma - 1))

# Set the v and B fields
v = np.zeros((3, N))
B = np.zeros((3, N))

r_par = pos[:,0] * np.cos(theta * np.pi / 180.) + pos[:,1] * np.sin(theta * np.pi / 180.)

v_per = v0 * np.sin(2. * np.pi * r_par)
B_per = B0 * np.sin(2. * np.pi * r_par)

v[0, :] = v_par * np.cos(theta * np.pi / 180.) - v_per * np.sin(theta * np.pi / 180.)
v[1, :] = v_par * np.sin(theta * np.pi / 180.) + v_per * np.cos(theta * np.pi / 180.)
v[2, :] = v0 * np.cos(2. * np.pi * r_par)

B[0, :] = B_par * np.cos(theta * np.pi / 180.) - B_per * np.sin(theta * np.pi / 180.)
B[1, :] = B_par * np.sin(theta * np.pi / 180.) + B_per * np.cos(theta * np.pi / 180.)
B[2, :] = B0 * np.cos(2. * np.pi * r_par)

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [Lx, Ly, 0.]
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
#grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")

fileOutput.close()
