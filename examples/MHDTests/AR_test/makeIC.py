###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2021 Federico Stasyszyn (fstasyszyn@unc.edu.ar)
#               2023 Orestis Karapiperis (karapiperis@strw.leidenuniv.nl)
#
# This file generates an IC file to run a 3D Brio & Wu shock tube
# (Brio & Wu, 1998 - https://ui.adsabs.harvard.edu/abs/1988JCoPh..75..400B/abstract).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

import argparse as ap
import h5py
import numpy as np

# Problem Parameters
gamma = 5.0 / 3.0 # Gas adiabatic index
rho   = 1.0  # Density
P = 1.0  # Pressure 
v0 = 0.1
B0 = 0.01  # Magnetic field

# File to save ICs to
fileName = "AR_test.hdf5"

# Simulation box attributes
boxSize = 1.0
vol = boxSize ** 3

# Stack BCC cubes
glass = h5py.File("BCCglassCube_32.hdf5", "r")

pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

numPart = np.size(h)

# Initialise extra arrays
v = np.zeros((numPart, 3))
B = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.ones(numPart)
u = np.ones(numPart)

# Instatiate the extra arrays
x = pos[:, 0]

v[:, 0] = v0 * np.sin(2.0 * np.pi * x)
B[:, 1] = B0 * np.sin(np.pi * x)

m *= rho * vol / numPart
u *= P / (rho * (gamma - 1.0))

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize, boxSize, boxSize]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")

file.close()
