###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Generates a distorted 2D grid with a density profile that balances out the
# external sine wave potential if run with an isothermal equation of state.

import numpy as np
import h5py
import unyt
import swiftsimio as sw

# constant thermal energy
# the definition below assumes the same thermal energy is defined in const.h,
# and also that the code was configured with an adiabatic index of 5./3.
uconst = 20.2615290634
cs2 = 2.0 * uconst / 3.0
A = 10.0

fileName = "sineWavePotential.hdf5"
numPart_1D = 50
boxSize = [1.0, 1.0]
numPart = numPart_1D**2

coords = np.zeros((numPart, 3))
v = np.zeros((numPart, 3))
m = np.zeros(numPart) + 1.0
h = np.zeros(numPart) + 2.0 / numPart
u = np.zeros(numPart) + uconst
ids = np.arange(1, numPart + 1, dtype="L")
rho = np.zeros(numPart)

# first set the positions, as we try to do a reasonable volume estimate to
# set the masses
for i in range(numPart_1D):
    for j in range(numPart_1D):
        coords[numPart_1D * i + j, 0] = (i + 0.5) / numPart_1D
        coords[numPart_1D * i + j, 1] = (j + 0.5) / numPart_1D

V = 1.0 / numPart
for i in range(numPart):
    rho[i] = 1000.0 * np.exp(
        -0.5 * A / np.pi / cs2 * np.cos(2.0 * np.pi * coords[i, 0])
    )
    m[i] = rho[i] * V

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        boxSize, unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=2,
)

w.gas.coordinates = coords * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids

w.write(fileName)

# Extra fields not handled by the Writer
with h5py.File(fileName, "r+") as f:
    f["PartType0"].create_dataset("Density", data=rho, dtype="f")
