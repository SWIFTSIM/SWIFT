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
"""Linear sound wave in a 3D glass (periodic unit box).

A plane wave along x with n wavelengths across the box (k = 2 pi n / L),
relative amplitude A: rho = rho0 (1 + A cos kx), P = P0 (1 + gamma A cos kx)
(adiabatic), with rho0 = 1 and P0 = 1 / gamma so that the sound speed is 1 and
the period of the n = 1 wave is 1. The density perturbation is imposed by the
Lagrangian displacement xi = -(A / k) sin kx of the glass particles.

  travelling (default): v_x = A cs cos kx, the wave moves in +x at the sound
      speed; the phase of the mode measures the numerical dispersion, its
      amplitude the damping.
  --standing: v = 0, the wave oscillates in place with the frequency cs k;
      the decay of the mode energy over many periods measures the damping.

A bulk velocity --boost adds a uniform v_x (Galilean-invariance check).
"""

import argparse

import h5py
import unyt
import swiftsimio as sw
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--glass", default="glassCube_32.hdf5", help="3D glass file")
parser.add_argument(
    "--amplitude", type=float, default=1e-3, help="relative density amplitude A"
)
parser.add_argument("--nwave", type=int, default=1, help="wavelengths across the box")
parser.add_argument(
    "--standing",
    action="store_true",
    help="standing wave (v = 0) instead of a travelling one",
)
parser.add_argument(
    "--boost",
    type=float,
    default=0.0,
    help="uniform bulk velocity along x in units of cs",
)
parser.add_argument("-o", "--output", default="soundWave.hdf5")
args = parser.parse_args()

gamma = 5.0 / 3.0
rho0 = 1.0
cs = 1.0
P0 = rho0 * cs**2 / gamma
u0 = P0 / ((gamma - 1.0) * rho0)
A = args.amplitude

glass = sw.load(args.glass)
pos = glass.gas.coordinates.value
h = glass.gas.smoothing_lengths.value
L = float(glass.metadata.boxsize[0].value)
N = len(h)
k = 2.0 * np.pi * args.nwave / L

# Lagrangian displacement: drho / rho = -d xi / dx = A cos kx
x0 = pos[:, 0]
pos[:, 0] = np.mod(x0 - (A / k) * np.sin(k * x0), L)
x = pos[:, 0]

m = np.ones(N) * rho0 * L**3 / N
u = u0 * (1.0 + (gamma - 1.0) * A * np.cos(k * x))
v = np.zeros((N, 3))
if not args.standing:
    v[:, 0] = A * cs * np.cos(k * x)
v[:, 0] += args.boost * cs
ids = np.arange(1, N + 1)

w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [L, L, L], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
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

# The parameters of the wave, for the analysis scripts (not handled by the
# Writer)
with h5py.File(args.output, "r+") as f:
    grp = f["Header"]
    grp.attrs["SoundWaveAmplitude"] = A
    grp.attrs["SoundWaveNWave"] = args.nwave
    grp.attrs["SoundWaveStanding"] = int(args.standing)
    grp.attrs["SoundWaveBoost"] = args.boost * cs
    grp.attrs["SoundWaveSoundSpeed"] = cs

kind = "standing" if args.standing else "travelling"
print(
    f"{kind} sound wave: {N} particles, A = {A}, {args.nwave} wavelength(s) per box "
    f"({round(N ** (1 / 3)) / args.nwave:.0f} particles per wavelength), period {L / (cs * args.nwave)}"
)
