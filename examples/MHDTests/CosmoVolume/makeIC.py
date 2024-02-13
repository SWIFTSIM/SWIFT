#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Parameters
shift_ICs = False

B0 = (
    1.56908e-3
)  # Physical magnetic field at z=63 corresponding to a comoving seed of 1e-12 Gauss (g/A*s^2), in units of 1e10 M_sun/(1e10 A * (Mpc/ 1e5 km/s))
UI = 1e10  # Unit current

# Files to read from and write to
fileInputName = sys.argv[1]
fileOutputName = sys.argv[2]

# Set up magnetic field and vector potential

infile = h5py.File(fileInputName, "r")

head = infile["/Header"]
pos_in = infile["/PartType0/Coordinates"][:, :]

BoxSize = head.attrs["BoxSize"]
afact = head.attrs["Time"]
N_in = head.attrs["NumPart_Total"][0]
infile.close()

wavelen = BoxSize / 10.0
wavenum = 2.0 * np.pi / wavelen

B = np.zeros((N_in, 3))
A = np.zeros((N_in, 3))

A0 = B0 / wavenum

A[:, 0] = A0 * (np.sin(pos_in[:, 2] * wavenum) + np.cos(pos_in[:, 1] * wavenum))
A[:, 1] = A0 * (np.sin(pos_in[:, 0] * wavenum) + np.cos(pos_in[:, 2] * wavenum))
A[:, 2] = A0 * (np.sin(pos_in[:, 1] * wavenum) + np.cos(pos_in[:, 0] * wavenum))
B[:, :] = wavenum * A[:, :]

os.system("cp " + fileInputName + " " + fileOutputName)

# File
fileOutput = h5py.File(fileOutputName, "a")

# Particle group
grp = fileOutput.require_group("/PartType0")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")

# shifting ICs
if shift_ICs:
    shift_vector = np.array([0.5, 0.5, 0]) * BoxSize
    pos_in += shift_vector
    for k in range(3):
        mask_right = pos_in[:, k] > BoxSize
        mask_left = pos_in[:, k] < 0
        pos_in[mask_right][:, k] -= BoxSize
        pos_in[mask_left][:, k] += BoxSize
    fileOutput["/PartType0/Coordinates"][:, :] = pos_in


# Change current unit to something more resonable
unitSystem = fileOutput["/Units"]
unitSystem.attrs.modify("Unit current in cgs (U_I)", UI)

fileOutput.close()
