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
"""Shearing shock (Cullen & Dehnen 2010, section 5.3).

Two streams of uniform gas collide head-on at x = L/2 (v_x = +dv for x < L/2,
-dv beyond, periodic box, so the streams also recede from each other at the
box edge) while a strong sinusoidal shear v_y = s sin(2 pi x / L) with
s = 100 dv is superposed. The shear carries no information about the shock:
in the exact solution every particle keeps its v_y (no pressure gradient in
y) and the x-dynamics is the 1D Riemann problem of the collision. The test
asks whether a shear-aware viscosity switch or limiter still captures the
shock when the shear rate exceeds the compression rate by two orders of
magnitude, and how much shear momentum is diffused in the process.

The box is N unit glass cubes along x, with dv = sound speed (Mach 1 streams).
"""

import argparse

import unyt
import swiftsimio as sw
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "-n", "--ncubes", type=int, default=4, help="number of glass cubes along x"
)
parser.add_argument(
    "--glass", default="glassCube_32.hdf5", help="3D glass file (unit cube)"
)
parser.add_argument(
    "--dv", type=float, default=1.0, help="stream speed (the sound speed is 1)"
)
parser.add_argument(
    "-s",
    "--shear",
    type=float,
    default=10.0,
    help="amplitude of v_y (paper: 100 dv, which heats the gas by 10-100x through the diffused shear in every SPH scheme)",
)
parser.add_argument("-o", "--output", default="shearingShock.hdf5")
args = parser.parse_args()

gamma = 5.0 / 3.0
rho0 = 1.0
cs = 1.0
P0 = rho0 * cs**2 / gamma

glass = sw.load(args.glass)
pos_g = glass.gas.coordinates.value
h_g = glass.gas.smoothing_lengths.value
n_g = len(h_g)
L = float(args.ncubes)
pos = np.concatenate([pos_g + np.array([i, 0.0, 0.0]) for i in range(args.ncubes)])
h = np.tile(h_g, args.ncubes)
N = len(h)
x = pos[:, 0]

v = np.zeros((N, 3))
v[:, 0] = np.where(x < 0.5 * L, args.dv, -args.dv)
v[:, 1] = args.shear * np.sin(2.0 * np.pi * x / L)
m = np.ones(N) * rho0 / n_g
u = np.ones(N) * P0 / (rho0 * (gamma - 1.0))
ids = np.arange(1, N + 1)

w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [L, 1.0, 1.0], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)
w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = np.asarray(ids).astype(np.uint64)
w.write(args.output)

print(
    f"Shearing shock: {N} particles in a {L:g} x 1 x 1 box, dv = {args.dv}, shear amplitude s = {args.shear}"
)
