###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

import h5py
import numpy as np
from shutil import copyfile

# Add sink particles to the isolated galaxy

fileName = "3e11-star-only-DM-halo-galaxy.hdf5"
output = "isolated_galaxy.hdf5"

L = 13756  # Size of the box (internal units)
L_sink = 1000 # Size of the sink particle area (L_sink < L)
N = 1000  # Number of sink particles
max_vel = 100  # Maximal velocity of the sink particles (in km / s)
mass = 0.000142  # Mass of a sink particle (internal units)
min_id = 360001  # minimal ids of the other particles

# ---------------------------------------------------

copyfile(fileName, output)

# File
file = h5py.File(output, 'a')

pos = np.random.rand(N, 3) * L_sink + 0.5 * (L - L_sink)
vel = 2 * (np .random.rand(N, 3) - 0.5) * max_vel
m = mass * np.ones(N)
# Generate extra arrays
ids = np.linspace(min_id, min_id + N, N)

# --------------------------------------------------

# Header
grp = file["Header"]
numpart = grp.attrs["NumPart_Total"][:]
numpart[3] = N
grp.attrs["NumPart_Total"] = numpart
grp.attrs["NumPart_ThisFile"] = numpart

# Units
grp = file["Units"]
uv = grp.attrs["Unit length in cgs (U_L)"]
uv /= grp.attrs["Unit time in cgs (U_t)"]

vel *= 1e5 / uv

# Particle group
grp = file.create_group("/PartType3")
grp.create_dataset('Coordinates', data=pos, dtype='d')
grp.create_dataset('Velocities', data=vel, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')

file.close()
