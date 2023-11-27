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
afact   = head.attrs["Time"]
N_in = head.attrs["NumPart_Total"][0]
print(BoxSize)
print(N_in)
infile.close()

Bini = 1E-12

B = np.zeros((N_in,3))
A = np.zeros((N_in,3))


wavelen = BoxSize / 10.0
wavenum = 2.0 * np.pi / wavelen
Aini = Bini / wavenum ;

print(wavelen,wavenum)

A[:,0] = Aini * (np.sin(pos_in[:,2]*wavenum) + np.cos(pos_in[:,1]*wavenum))
A[:,1] = Aini * (np.sin(pos_in[:,0]*wavenum) + np.cos(pos_in[:,2]*wavenum))
A[:,2] = Aini * (np.sin(pos_in[:,1]*wavenum) + np.cos(pos_in[:,0]*wavenum))
B[:,0] = Bini * (np.sin(pos_in[:,2]*wavenum) + np.cos(pos_in[:,1]*wavenum))
B[:,1] = Bini * (np.sin(pos_in[:,0]*wavenum) + np.cos(pos_in[:,2]*wavenum))
B[:,2] = Bini * (np.sin(pos_in[:,1]*wavenum) + np.cos(pos_in[:,0]*wavenum))

print(min(A[:,0]),max(A[:,0]))
print(min(B[:,0]),max(B[:,0]))

#B[:,0] = np.where(pos_in[:,1] < BoxSize * 0.5, Bini, -Bini)
#B[:,1] = 0.0 
#B[:,2] = 0.0
#A[:,0] = 0.0
#A[:,1] = 0.0 
#A[:,2] = np.where(pos_in[:,1] < BoxSize *0.5, Bini * pos_in[:,1], Bini * (BoxSize-pos_in[:,1]))

#print(min(A[:,2]),max(A[:,2]))
#print(min(B[:,0]),max(B[:,0]))

os.system("cp "+fileInputName+" "+fileOutputName)
# File
fileOutput = h5py.File(fileOutputName, "a")

# Units
#print(list(fileOutput.keys()))
#grpu = fileOutput.require_group(["/Units"])
grpu = fileOutput["/Units"]
print(grpu.attrs["Unit length in cgs (U_L)"])
print(grpu.attrs["Unit mass in cgs (U_M)"])
print(grpu.attrs["Unit time in cgs (U_t)"]) 
print(grpu.attrs["Unit current in cgs (U_I)"])
print(grpu.attrs["Unit temperature in cgs (U_T)"])

grpu.attrs["Unit current in cgs (U_I)"] = 4.688e6  # amperes to be Gauss

# Particle group
grp = fileOutput.require_group("/PartType0")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")

fileOutput.close()
