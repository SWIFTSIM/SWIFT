###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
import numpy as np
import h5py as h5

box_size = 1000.0
N_PARTICLES = 1
print("Initial conditions written to 'IC.hdf5'")

pos = np.zeros((1, 3))
pos[0, 0] = 0.1
pos[0, 1] = 0.1
pos[0, 2] = 100


# pos += np.array(
#    [box_size / 2, box_size / 2, box_size / 2]
# )  # shifts the particles to the center of the box

vel = np.zeros((1, 3))
vel[0, 0] = 80
vel[0, 1] = 0
vel[0, 2] = 0


ids = np.array([1.0])
mass = np.array([1.0]) * 1e-10

# File
file = h5.File("IC.hdf5", "w")

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 3.086e21
grp.attrs["Unit mass in cgs (U_M)"] = 1.98848e43
grp.attrs["Unit time in cgs (U_t)"] = 3.086e16
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = box_size
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
