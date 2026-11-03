################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

import h5py
import numpy as np
import unyt
import swiftsimio as sw

# Generates a swift IC file for the Kelvin-Helmholtz vortex in a periodic box

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
P0 = 2.5  # Pressure
rho0 = 1.0  # Density
d = 0.0317  # Thickness of the transition layer
B = 0.0005  # Amplitude of the seed velocity

fileOutputName = "kelvinHelmholtzGrowthRate.hdf5"

# ---------------------------------------------------

glass = h5py.File("glassPlane_128.hdf5", "r")
pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

N = len(h)
vol = 1.0

# Generate extra arrays
v = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho0 * vol / N
u = np.ones(N) * P0 / (rho0 * (gamma - 1.0))

v[pos[:, 1] <= 0.25 - d, 0] = -0.5
v[(pos[:, 1] < 0.25 + d) & (pos[:, 1] > 0.25 - d), 0] = (
    -0.5
    + 0.5 * (pos[(pos[:, 1] < 0.25 + d) & (pos[:, 1] > 0.25 - d), 1] + d - 0.25) / d
)
v[(pos[:, 1] <= 0.75 - d) & (pos[:, 1] >= 0.25 + d), 0] = 0.5
v[(pos[:, 1] < 0.75 + d) & (pos[:, 1] > 0.75 - d), 0] = (
    0.5 - 0.5 * (pos[(pos[:, 1] < 0.75 + d) & (pos[:, 1] > 0.75 - d), 1] + d - 0.75) / d
)
v[pos[:, 1] >= 0.75 + d, 0] = -0.5

v[:, 1] = (
    B
    * np.sin(4.0 * np.pi * pos[:, 0])
    * (
        np.exp(-((pos[:, 1] - 0.25) ** 2) / 32.0 / d**2)
        + np.exp(-((pos[:, 1] - 0.75) ** 2) / 32.0 / d**2)
    )
)

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [1.0, 1.0], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=2,
)

w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileOutputName)
