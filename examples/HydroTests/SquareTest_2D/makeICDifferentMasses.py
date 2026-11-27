###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
#               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

import unyt
import swiftsimio as sw
import numpy as np

# Generates a swift IC file for the Square test in a periodic box

# Parameters
L = 2 * 64  # Number of particles on the side
gamma = 5.0 / 3.0  # Gas adiabatic index
rho0 = 4.0  # Gas central density
rho1 = 1.0  # Gas outskirt density
P0 = 2.5  # Gas central pressure
P1 = 2.5  # Gas central pressure
vx = 0.0  # Random velocity for all particles
vy = 0.0
fileOutputName = "square.hdf5"
# ---------------------------------------------------

vol = 1.0

numPart = L * L

pos_x = np.arange(0, 1, 1.0 / L)
xv, yv = np.meshgrid(pos_x, pos_x)
pos = np.zeros((numPart, 3), dtype=float)
pos[:, 0] = xv.flatten()
pos[:, 1] = yv.flatten()

# Now we can get 2d masks!
inside = np.logical_and.reduce([xv < 0.75, xv > 0.25, yv < 0.75, yv > 0.25])

mass_in = rho0 / numPart
mass_out = rho1 / numPart

m = np.ones_like(xv) * mass_out
m[inside] = mass_in
m = m.flatten()

h = np.ones_like(m) / L

u_in = P0 / ((gamma - 1) * rho0)
u_out = P1 / ((gamma - 1) * rho1)
u = np.ones_like(xv) * u_out
u[inside] = u_in
u = u.flatten()

vel = np.zeros_like(pos)

ids = np.arange(1, numPart + 1)


# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [vol, vol], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=2,
)

w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = vel * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileOutputName)
