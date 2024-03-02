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
    ddefault="CPLglassCube_32x_16y_16z.hdf5",
    type=str,
)

args = parser.parse_args()

# Parameters
rho = 1.0
P = 0.1
gamma = 5.0 / 3.0

v0 = 0.1
B0 = 0.1
wavelen = 1.0 
wavenum = 2.0 * np.pi / wavelen

fileOutputName = "CircularlyPolarisedAlfvenWave.hdf5"

# Retrieve positions file
glass = h5py.File(args.positionsFile, "r")

pos = glass["/PartType0/Coordinates"][:, :]
h   = glass["/PartType0/SmoothingLength"][:]
N = len(h)

posFileBoxSize = glass["/Header"].attrs["BoxSize"]

print("The positions' file box size is : ", posFileBoxSize)
print("The shape of the positions array is : ", pos.shape) 

# Rescale positions and smoothing lenths to desired scale
Lx = 2.0                # 2.0 / np.sqrt(3.0)        # 3.0
Ly = 0.5 * np.sqrt(3.0) # 2.0                       # scaling * posFileBoxSize[1]
Lz = np.sqrt(2.0/3.0)   # 2.0 * np.sqrt(2.0) / 3.0  # scaling * posFileBoxSize[2]

scaling = [Lx / posFileBoxSize[0], Ly / posFileBoxSize[1], Lz / posFileBoxSize[2]]

print("The scaling is : ", scaling)

pos[:,0] *= scaling[0]
pos[:,1] *= scaling[1]
pos[:,2] *= scaling[2]

h *= scaling[0]

vol = Lx * Ly * Lz

print("The box size is : ", Lx, Ly, Lz)

# Initialise and instantiate additional arrays

sina = 0.0 
cosa = np.sqrt(1 - sina ** 2)

sinb = 0.0 # 0.5
cosb = np.sqrt(1 - sinb ** 2)

print("sina and sinb are : ", sina, sinb)

Rot_inv = np.array(
    [
        [cosa * cosb, -sinb, -sina * cosb],
        [cosa * sinb, cosb, -sina * sinb],
        [sina, 0, cosa],
    ]
)

Rot = Rot_inv.T

#x1 = (pos[:, 0] + 2 * pos[:, 1] + 2 * pos[:, 2]) / 3
x1 = np.dot(Rot, pos.T)[0,:]

print("The max (x,y,z) is : ", max(pos[:,0]), max(pos[:,1]), max(pos[:,2]))
print("The max x1 is : ", max(x1))

v = np.zeros((3, N))
B = np.zeros((3, N))
A = np.zeros((3, N))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho * vol / N
u = np.ones(N) * P / (rho * (gamma - 1))

v[0, :] = 0.0
v[1, :] = v0 * np.sin(wavenum * x1)
v[2, :] = v0 * np.cos(wavenum * x1)

B[0, :] = 1.0
B[1, :] = B0 * np.sin(wavenum * x1)
B[2, :] = B0 * np.cos(wavenum * x1)

A[1, :] = B0 / wavenum * np.sin(wavenum * x1)
A[2, :] = B0 / wavenum * np.cos(wavenum * x1)

v = np.dot(Rot_inv, v)
B = np.dot(Rot_inv, B)
A = np.dot(Rot_inv, A)

v = v.T
B = B.T
A = A.T

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
