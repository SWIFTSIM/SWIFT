#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

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

# Other variables
wavelen = BoxSize / 10.0
wavenum = 2.0 * np.pi / wavelen

B = np.zeros((N_in, 3))
A = np.zeros((N_in, 3))

# Set the magnitude of the uniform seed magnetic field
B0_Gaussian_Units = 1e-3 # 1e-3 micro Gauss
CONST_MU0_CGS = 4 * np.pi
B0_cgs = np.sqrt(CONST_MU0_CGS / (4.0 * np.pi)) * B0_Gaussian_Units
B[:,0] = B0_cgs # gauss 

# Now set the magnetic field unit
const_unit_magnetic_field_in_cgs = 1e-7  # 1muG
const_unit_mass_in_cgs = 1.9891e43
const_unit_length_in_cgs = 3.08567758e21
const_unit_velocity_in_cgs = 1e5

# Derived quantities
const_unit_time_in_cgs = const_unit_length_in_cgs / const_unit_velocity_in_cgs
const_unit_current_in_cgs = const_unit_mass_in_cgs / (
    const_unit_magnetic_field_in_cgs * const_unit_time_in_cgs ** 2
)
UI = const_unit_current_in_cgs

#A0 = B0 / wavenum * afact

#A[:, 0] = A0 * (np.sin(pos_in[:, 2] * wavenum) + np.cos(pos_in[:, 1] * wavenum))
#A[:, 1] = A0 * (np.sin(pos_in[:, 0] * wavenum) + np.cos(pos_in[:, 2] * wavenum))
#A[:, 2] = A0 * (np.sin(pos_in[:, 1] * wavenum) + np.cos(pos_in[:, 0] * wavenum))
#B[:, 0] = B0 * (np.sin(pos_in[:, 2] * wavenum) + np.cos(pos_in[:, 1] * wavenum))
#B[:, 1] = B0 * (np.sin(pos_in[:, 0] * wavenum) + np.cos(pos_in[:, 2] * wavenum))
#B[:, 2] = B0 * (np.sin(pos_in[:, 1] * wavenum) + np.cos(pos_in[:, 0] * wavenum))

os.system("cp " + fileInputName + " " + fileOutputName)

# File
fileOutput = h5py.File(fileOutputName, "a")

# Particle group
grp = fileOutput.require_group("/PartType0")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")


# Change current unit to something more resonable
unitSystem = fileOutput["/Units"]
unitSystem.attrs.modify("Unit current in cgs (U_I)", UI)
fileOutput.close()
