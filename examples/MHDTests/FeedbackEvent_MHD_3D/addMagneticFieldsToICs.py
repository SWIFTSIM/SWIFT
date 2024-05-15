#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt

import unyt

'''
cosmo_units = unyt.UnitSystem(
        "cosmological",
        unyt.Mpc,
        unyt.unyt_quantity(1e10, units=unyt.Solar_Mass),
        unyt.unyt_quantity(1.0, units=unyt.s * unyt.Mpc / unyt.km).to(unyt.Gyr),
        unyt.K,
        unyt.rad,
        unyt.unyt_quantity(1e10, units=unyt.A),
)

# Parameters
B0 = 1e-13 * unyt.g / (unyt.statA * unyt.s * unyt.s) # Specified in terms of Gauss
B0.convert_to_base(cosmo_units)
'''

B0 = 0.004788438032
UI = 1e10  # Unit current

# File to read from
file = h5py.File(sys.argv[1], "r+")

parts = file["/PartType0"]
Npart = parts["ParticleIDs"].len()

# Instantiate magnetic field
B = np.zeros((Npart, 3))
B[:, 2] = B0

B_exists = "MagneticFluxDensities" in parts.keys()

if B_exists:
    magneticFluxDensities = parts["MagneticFluxDensities"]
    magneticFluxDensities[...] = B
else:
    parts.create_dataset("MagneticFluxDensities", data=B, dtype="f")

# Change current unit to something more resonable
unitSystem = file["/Units"]
unitSystem.attrs.modify("Unit current in cgs (U_I)", UI)

file.close()
