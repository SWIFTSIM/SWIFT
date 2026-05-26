#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import sys
import os

# Unit current
UI = 2.088357e14

# Physical seed magnetic field in units of nG (Marinacci 2015, Pillepich 2018)
B0 = 0.01

# Files to read from and write to
fileInputName = sys.argv[1]
fileOutputName = sys.argv[2]

# Read in relevant information                                                                                                                                               
infile = h5py.File(fileInputName, "r")
head = infile["/Header"]
N_in = head.attrs["NumPart_Total"][0]
infile.close()

# Set magnetic field
B = np.zeros((N_in, 3))
B[:,2] = B0

os.system("cp " + fileInputName + " " + fileOutputName)

# File
fileOutput = h5py.File(fileOutputName, "r+")

# Set seed B field
parts = fileOutput["/PartType0"]
B_exists = "MagneticFluxDensities" in parts.keys()

if B_exists:
    magneticFluxDensities = parts["MagneticFluxDensities"]
    magneticFluxDensities[...] = B
else:
    parts.create_dataset("MagneticFluxDensities", data=B, dtype="f")

# Change current unit to something more resonable
unitSystem = fileOutput["/Units"]
unitSystem.attrs.modify("Unit current in cgs (U_I)", UI)

fileOutput.close()
