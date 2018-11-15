###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
################################################################################
from galpy.potential import NFWPotential
from galpy.orbit import Orbit
from galpy.util import bovy_conversion
import numpy as np
import matplotlib.pyplot as plt
from astropy import units
import h5py as h5

C = 8.0
M_200 = 2.0
N_PARTICLES = 3
print("Initial conditions written to 'test_nfw.hdf5'")

pos = np.zeros((3, 3))
pos[0, 2] = 50.0
pos[1, 2] = 10.0
pos[2, 2] = 2.0
pos = pos + 500.0
vel = np.zeros((3, 3))
vel[0, 1] = 348.0
vel[1, 1] = 466.9
vel[2, 1] = 348.0
ids = np.array([1.0, 2.0, 3.0])
mass = np.array([1.0, 1.0, 1.0])

# File
file = h5.File("circularorbitshernquist.hdf5", "w")

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 3.086e21
grp.attrs["Unit mass in cgs (U_M)"] = 1.988e33
grp.attrs["Unit time in cgs (U_t)"] = 3.086e16
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = 1000.0
grp.attrs["NumPart_Total"] = [0, N_PARTICLES, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, N_PARTICLES, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = 1

# Particle group
grp1 = file.create_group("/PartType1")
ds = grp1.create_dataset("Velocities", (N_PARTICLES, 3), "f", data=vel)

ds = grp1.create_dataset("Masses", (N_PARTICLES,), "f", data=mass)

ds = grp1.create_dataset("ParticleIDs", (N_PARTICLES,), "L", data=ids)

ds = grp1.create_dataset("Coordinates", (N_PARTICLES, 3), "d", data=pos)

file.close()
