################################################################################
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
################################################################################
"""Keplerian ring of Cullen & Dehnen (2010, MNRAS 408, 669), Section 4.3.

A gaseous ring with Gaussian surface density Sigma(r) ~ exp(-(r - r0)^2 /
(2 w^2)), r0 = 10, w = 2.5, orbits a central point mass with GM = 1000 (orbital
period 2 pi and velocity 10 at the ring centre). The gas is cold, c_s = 0.01 =
1e-3 v_phi, so that pressure forces are negligible and the ring should stay in
its initial configuration; numerical shear viscosity spreads it and breaks it
up through the viscous instability. N = 9745 equal-mass particles are placed
with the method of Cartwright, Stamatellos & Whitworth (2009): radii from the
inverse of the cumulative mass profile, azimuths on a golden-angle spiral, so
that the particle distribution is quasi-uniform at the local density without
lattice directions. The paper follows the ring for more than six orbits
(t = 40).

SWIFT's point-mass potential is softened; the default softening of 0.1 (1% of
the ring radius) is negligible for the ring and the circular velocities are
those of the softened potential. The ring sits in a 50 x 50 periodic box. The
Gaussian is truncated at 2 widths by default (95% of the mass, r in [5, 15]):
the paper's untruncated tails place isolated particles many spacings apart,
which the kernel cannot resolve and which fall onto the point mass and break
the ring (see --truncate).

Usage: makeIC.py [-n 9745] [-c 0.01] [-e 0.1] [-m 1000] [--truncate 2] [-o initial_conditions.hdf5]
"""
import argparse

import h5py
import numpy as np

parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
parser.add_argument("-n", "--nparts", type=int, default=9745, help="number of particles (paper: 9745)")
parser.add_argument("-c", "--soundspeed", type=float, default=0.01, help="sound speed (paper: 0.01)")
parser.add_argument("-e", "--softening", type=float, default=0.1,
                    help="softening of the point mass (the paper has none); PointMassPotential:softening")
parser.add_argument("-m", "--gravitymass", type=float, default=1000.0, help="GM of the point mass (paper: 1000)")
parser.add_argument("-r", "--radius", type=float, default=10.0, help="ring radius (paper: 10)")
parser.add_argument("-w", "--width", type=float, default=2.5, help="Gaussian width of the ring (paper: 2.5)")
parser.add_argument("-L", "--boxsize", type=float, default=50.0, help="box size")
parser.add_argument("-g", "--gamma", type=float, default=5.0 / 3.0, help="adiabatic index")
parser.add_argument("--truncate", type=float, default=2.0,
                    help="truncate the Gaussian at this many widths from the ring centre (0: no truncation). "
                         "The paper's ring has no truncation, but the isolated particles of its tails "
                         "(surface density below a few per cent of the peak, many spacings apart) cannot be "
                         "resolved by the kernel: they orbit much faster than the ring, fall onto the point "
                         "mass within an orbit and act on the ring's inner edge through kernels the size of "
                         "the gap (CD10 report the resulting central features as artefacts; here they break "
                         "the ring). The default of 2 widths keeps 95%% of the mass as a resolved annulus; "
                         "3 widths (0.3%% of the mass removed) still leaves isolated tail particles.")
parser.add_argument("--hcap", type=float, default=3.0,
                    help="cap on the initial smoothing length (Gaussian tails); keep SPH:h_max consistent")
parser.add_argument("-o", "--output", default="initial_conditions.hdf5")
args = parser.parse_args()

N = args.nparts
r0, w, L = args.radius, args.width, args.boxsize
GM, eps = args.gravitymass, args.softening
centre = 0.5 * L

# --- Radii: invert the cumulative mass M(<r) of Sigma(r) r dr numerically
n_w = args.truncate if args.truncate > 0 else 6.0
r_grid = np.linspace(max(0.0, r0 - n_w * w), r0 + n_w * w, 200001)
sigma = np.exp(-0.5 * ((r_grid - r0) / w) ** 2)
cum = np.concatenate(([0.0], np.cumsum(0.5 * (sigma[1:] * r_grid[1:] + sigma[:-1] * r_grid[:-1]) * np.diff(r_grid))))
cum /= cum[-1]
mass_fraction = (np.arange(N) + 0.5) / N
r = np.interp(mass_fraction, cum, r_grid)

# --- Azimuths: golden-angle spiral (Cartwright, Stamatellos & Whitworth 2009)
golden_angle = np.pi * (3.0 - np.sqrt(5.0))
theta = np.mod(np.arange(N) * golden_angle, 2.0 * np.pi)

x = centre + r * np.cos(theta)
y = centre + r * np.sin(theta)
pos = np.zeros((N, 3))
pos[:, 0], pos[:, 1] = x, y

# --- Circular velocities of the softened point-mass potential
v_phi = np.sqrt(GM * r * r / (r * r + eps * eps) ** 1.5)
vel = np.zeros((N, 3))
vel[:, 0] = -v_phi * np.sin(theta)
vel[:, 1] = v_phi * np.cos(theta)

# --- Equal masses (total ring mass 1, irrelevant without self-gravity),
#     internal energy from the sound speed, smoothing length from the local
#     surface density
m = np.full(N, 1.0 / N)
u = np.full(N, args.soundspeed**2 / (args.gamma * (args.gamma - 1.0)))
# Sigma normalised to a ring mass of 1: Sigma(r) ~ exp(..) / (2 pi r0 w sqrt(2 pi)) for w << r0
sigma_local = np.exp(-0.5 * ((r - r0) / w) ** 2) / (2.0 * np.pi * r0 * w * np.sqrt(2.0 * np.pi))
# Capped for the few particles in the Gaussian tails (their neighbours are far
# apart; SWIFT's cell size follows the largest h and must fit the box)
h = np.minimum(1.2348 * np.sqrt(m / sigma_local), args.hcap)
ids = np.arange(1, N + 1)

v_centre = np.sqrt(GM / r0)
print(f"Cullen & Dehnen (2010) ring: N = {N}, r0 = {r0}, width = {w}, GM = {GM}")
print(f"  v_phi(r0) = {v_centre:.3f}, orbital period at r0 = {2.0 * np.pi * r0 / v_centre:.3f}, "
      f"c_s = {args.soundspeed} (Mach {v_centre / args.soundspeed:.0f}), softening = {eps}")
print(f"  Gaussian truncated at {n_w:g} widths: r in [{r.min():.2f}, {r.max():.2f}]")

with h5py.File(args.output, "w") as f:
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [L, L, L]
    grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFilesPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = 0
    grp.attrs["Dimension"] = 2

    grp = f.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 1.0
    grp.attrs["Unit mass in cgs (U_M)"] = 1.0
    grp.attrs["Unit time in cgs (U_t)"] = 1.0
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    grp = f.create_group("/PartType0")
    grp.create_dataset("Coordinates", data=pos, dtype="d")
    grp.create_dataset("Velocities", data=vel, dtype="f")
    grp.create_dataset("Masses", data=m, dtype="f")
    grp.create_dataset("SmoothingLength", data=h, dtype="f")
    grp.create_dataset("InternalEnergy", data=u, dtype="f")
    grp.create_dataset("ParticleIDs", data=ids, dtype="L")
