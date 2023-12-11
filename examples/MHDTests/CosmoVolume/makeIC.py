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


fileOutputName = sys.argv[2]
fileInputName = sys.argv[1]

###---------------------------###

infile = h5py.File(fileInputName, "r")

head = infile["/Header"]
pos_in = infile["/PartType0/Coordinates"][:, :]

BoxSize = head.attrs["BoxSize"]
afact = head.attrs["Time"]
N_in = head.attrs["NumPart_Total"][0]
print(BoxSize)
print(N_in)
infile.close()

Bini = 1e-12
wavelen = BoxSize / 10.0
wavenum = 2.0 * np.pi / wavelen

B = np.zeros((N_in, 3))
A = np.zeros((N_in, 3))

Aini = Bini / wavenum

A[:, 0] = Aini * (np.sin(pos_in[:, 2] * wavenum) + np.cos(pos_in[:, 1] * wavenum))
A[:, 1] = Aini * (np.sin(pos_in[:, 0] * wavenum) + np.cos(pos_in[:, 2] * wavenum))
A[:, 2] = Aini * (np.sin(pos_in[:, 1] * wavenum) + np.cos(pos_in[:, 0] * wavenum))
B[:, :] = wavenum * A[:, :]

os.system("cp " + fileInputName + " " + fileOutputName)
# File
fileOutput = h5py.File(fileOutputName, "a")


# Particle group
grp = fileOutput.require_group("/PartType0")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")

fileOutput.close()
