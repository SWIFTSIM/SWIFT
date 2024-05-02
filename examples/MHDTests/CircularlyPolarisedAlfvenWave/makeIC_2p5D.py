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

v_par = 0
B_par = 1
v0 = 0.1
B0 = 0.1
theta = 30 # degrees
Nside = 32

fileOutputName = "CircularlyPolarisedAlfvenWave.hdf5"

# Generate square grid
#N = Nside * Nside
#grid_x = np.linspace(0, 1 - 1./Nside, Nside)
#grid_y = np.linspace(0, 1 - 1./Nside, Nside)
#xv, yv = np.meshgrid(grid_x, grid_y)
#pos = np.array((xv.flatten(), yv.flatten(), np.zeros(N))).T
#h = np.ones(N) *  1.2 * (1. / Nside)
#glassBoxSize = np.array([1, 1, 0])

N = Nside * Nside
grid_x = np.linspace(1./(2 * Nside), 1 - 1./(2 * Nside), Nside)
grid_y = np.linspace(1./(2 * Nside), 1 - 1./(2 * Nside), Nside)
xv, yv = np.meshgrid(grid_x, grid_y)
pos = np.array((xv.flatten(), yv.flatten(), np.zeros(N))).T

# Shift every second row
mask = pos[:,1] % (2. / Nside) == (3. / (2 * Nside))
pos[mask, 0] += 1. / (2 * Nside)

# Shift everything to avoid particles on the very edge
pos[:, 0] -= 1. / (4 * Nside)

# Compute sensible values for h
h = np.ones(N) * (1. / (2 * Nside))

# Simulation volume
glassBoxSize = np.array([1, 1, 0])

# Duplicate along y
shifted = pos + np.array([0, 1, 0])
pos = np.append(pos, shifted, axis=0)
h = np.append(h, h)
glassBoxSize[1] *= 2.
N *= 2

print(pos)

# Expected dimensions
Lx = 1. / np.cos(theta * np.pi / 180.)
Ly = 1. / np.sin(theta * np.pi / 180.)

# Scaling factors
scaling = [Lx / glassBoxSize[0], Ly / glassBoxSize[1]]

print("The scaling is : ", scaling)

pos[:, 0] *= scaling[0]
pos[:, 1] *= scaling[1]
h *= np.sqrt(scaling[0]**2 + scaling[1]**2)

# Info about volume
vol = Lx * Ly 
print("The box size is : ", Lx, Ly)
print("The max (x,y,z) is : ", max(pos[:,0]), max(pos[:,1]), max(pos[:,2]))

# Set the masses, ids and energies
m = np.ones(N) * rho * vol / N / 0.995881
ids = np.linspace(1, N, N)
u = np.ones(N) * P / (rho * (gamma - 1))

# P = u * rho * (g - 1)
# P / (rho) = u * (g - 1)
# P / (rho * (g - 1)) = u

# Set the v and B fields
v = np.zeros((N, 3))
B = np.zeros((N, 3))

r_par = pos[:,0] * np.cos(theta * np.pi / 180.) + pos[:,1] * np.sin(theta * np.pi / 180.)

v_per = v0 * np.sin(2. * np.pi * r_par)
B_per = B0 * np.sin(2. * np.pi * r_par)

v[:, 0] = v_par * np.cos(theta * np.pi / 180.) - v_per * np.sin(theta * np.pi / 180.)
v[:, 1] = v_par * np.sin(theta * np.pi / 180.) + v_per * np.cos(theta * np.pi / 180.)
v[:, 2] = v0 * np.cos(2. * np.pi * r_par)

B[:, 0] = B_par * np.cos(theta * np.pi / 180.) - B_per * np.sin(theta * np.pi / 180.)
B[:, 1] = B_par * np.sin(theta * np.pi / 180.) + B_per * np.cos(theta * np.pi / 180.)
B[:, 2] = B0 * np.cos(2. * np.pi * r_par)

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
grp.attrs["Dimension"] = 2

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
