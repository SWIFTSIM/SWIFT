###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

"""Initial conditions for the driven subsonic turbulence test of
Bauer & Springel (2012, MNRAS 423, 2558; arXiv:1109.4413), section 3:
a unit periodic box of gas with rho = 1 and isothermal sound speed c_s = 1,
N^3 particles on a Cartesian lattice at rest. The gas is nearly isothermal
(gamma = 1.001, see the README), so the internal energy is u = c_s^2 / (gamma - 1).

Usage: python3 makeIC.py [N] [gamma] [eta]   (defaults: 64, 1.001, 1.2348)
Writes drivenTurbulence.hdf5.
"""

import sys
import numpy as np
import unyt
import swiftsimio as sw

N = int(sys.argv[1]) if len(sys.argv) > 1 else 64
gamma = float(sys.argv[2]) if len(sys.argv) > 2 else 1.001
eta = float(sys.argv[3]) if len(sys.argv) > 3 else 1.2348
fileName = "drivenTurbulence.hdf5"

L = 1.0  # Box size
rho = 1.0  # Density
cs = 1.0  # Isothermal sound speed
u = cs * cs / (gamma - 1.0)  # P = (gamma - 1) rho u = rho c_s^2

numPart = N**3
g = (np.arange(N) + 0.5) * L / N
x, y, z = np.meshgrid(g, g, g, indexing="ij")
pos = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)

w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [L, L, L], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = np.zeros((numPart, 3)) * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = np.full(numPart, rho * L**3 / numPart) * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = np.full(numPart, eta * L / N) * unyt.cm * w.a.comoving
w.gas.internal_energy = np.full(numPart, u) * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = np.arange(1, numPart + 1, dtype=np.uint64)

w.write(fileName)

print(f"Wrote {fileName}: {N}^3 particles, u = {u:g}, h = {eta * L / N:g}")
