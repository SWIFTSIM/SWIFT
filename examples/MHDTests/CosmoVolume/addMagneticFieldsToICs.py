#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt

# Parameters
B0 = 1.56908e-3 # Physical magnetic field at z=63 corresponding to a comoving seed of 1e-12 Gauss (g/A*s^2), in units of 1e10 M_sun/(1e10 A * (Mpc/ 1e5 km/s))
UI = 1e10       # Unit current

# File to read from
file = h5py.File(sys.argv[1], "r+")   

parts = file["/PartType0"]
Npart = parts["ParticleIDs"].len()

# Instantiate magnetic field    
B = np.zeros((Npart,3))
B[:,2] = B0

B_exists = "MagneticFluxDensities" in parts.keys()

if B_exists:
    magneticFluxDensities = parts["MagneticFluxDensities"]
    magneticFluxDensities[...] = B
else:  
    parts.create_dataset("MagneticFluxDensities", data=B, dtype="f")

# Change current unit to something more resonable    
unitSystem = file["/Units"]
unitSystem.attrs.modify('Unit current in cgs (U_I)', UI)

file.close()
