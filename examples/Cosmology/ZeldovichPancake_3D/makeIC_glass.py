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
"""Zel'dovich pancake on a glass instead of a lattice.

Same problem as makeIC.py (1D sinusoidal perturbation of wavelength 64 Mpc/h
in an Einstein-de Sitter universe of baryons, caustic at z_c = 1, temperature
T_i at the mean density at z = 100) but the Lagrangian positions q are those
of a glass file. This breaks the lattice symmetry which otherwise keeps the
particles in columns along the collapse direction (and SPH in its "column
regime" once the transverse spacing exceeds the kernel support).

A glass has particle-scale density noise, and the cold gas of this test is
Jeans-unstable at every resolved scale: at T_i = 100 K the comoving Jeans
length is ~1 kpc at z = 1, against a particle spacing of 2 Mpc. Making the
gas Jeans-stable until the caustic would need T_i ~ 5e10 K and reduce the
infall to Mach 3, so that is not an option. Instead the run starts later:
the noise grows linearly with the expansion factor, so starting at z_start
(default 10, zfac = 0.18) instead of z = 100 limits its growth to a factor
(1 + z_start) / (1 + z_end) ~ 6 instead of 50. The 1D Zel'dovich solution is
exact until the caustic whatever the starting redshift, and the thermal state
is the adiabatic one of the z = 100 start (T_i (1 + z_start)^2 / 101^2 at
mean density), so the analytic solution of the lattice test applies
unchanged.

Usage: makeIC_glass.py [--glass glassCube_32.hdf5] [--z_start 10]
                       [--z_end 0.9] [--T_i 100] [-o zeldovichPancake.hdf5]
The run must start at z_start: Cosmology:a_begin = 1 / (1 + z_start).
"""
import argparse

import h5py
import numpy as np

parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
parser.add_argument("--glass", default="glassCube_32.hdf5", help="glass file (unit cube)")
parser.add_argument("--z_start", type=float, default=10.0, help="starting redshift of the run")
parser.add_argument("--z_end", type=float, default=0.9, help="final redshift (diagnostics only)")
parser.add_argument("--T_i", type=float, default=100.0,
                    help="temperature at mean density at z = 100 in K (as makeIC.py)")
parser.add_argument("-o", "--output", default="zeldovichPancake.hdf5")
args = parser.parse_args()

# Parameters (as in makeIC.py)
T_i = args.T_i  # Temperature of the gas at mean density at z = 100 (in K)
z_c = 1.0  # Redshift of caustic formation (non-linear collapse)
z_ref = 100.0  # Redshift at which T_i is defined (start of the lattice test)
z_i = args.z_start  # Starting redshift of this run
gamma = 5.0 / 3.0  # Gas adiabatic index

# Some units
Mpc_in_m = 3.08567758e22
Msol_in_kg = 1.98848e30
Gyr_in_s = 3.08567758e19
mH_in_kg = 1.6737236e-27

# Some constants
kB_in_SI = 1.38064852e-23
G_in_SI = 6.67408e-11

# Some useful variables in h-full units
H_0 = 1.0 / Mpc_in_m * 10**5  # h s^-1
rho_0 = 3.0 * H_0**2 / (8 * np.pi * G_in_SI)  # h^2 kg m^-3
lambda_i = 64.0 / H_0 * 10**5  # h^-1 m (= 64 h^-1 Mpc)
x_min = -0.5 * lambda_i
x_max = 0.5 * lambda_i

# SI system of units
unit_l_in_si = Mpc_in_m
unit_m_in_si = Msol_in_kg * 1.0e10
unit_t_in_si = Gyr_in_s
unit_v_in_si = unit_l_in_si / unit_t_in_si
unit_u_in_si = unit_v_in_si**2

# ---------------------------------------------------

# Lagrangian positions from the glass
with h5py.File(args.glass, "r") as f:
    glass = f["/PartType0/Coordinates"][:]
    glass_box = f["/Header"].attrs["BoxSize"]
glass /= np.asarray(glass_box, dtype=float).ravel()[0]
glass = np.mod(glass, 1.0)
numPart = len(glass)
numPart_1D = int(round(numPart ** (1.0 / 3.0)))

boxSize = x_max - x_min
delta_x = boxSize / numPart_1D  # mean inter-particle separation
k_i = 2.0 * np.pi / lambda_i
zfac = (1.0 + z_c) / (1.0 + z_i)
m_i = boxSize**3 * rho_0 / numPart

# Temperature at mean density at the start (adiabatic from z_ref)
T_start = T_i * ((1.0 + z_i) / (1.0 + z_ref)) ** 2


# ---------------------------------------------------
# Diagnostics: Jeans scale of the noise and infall Mach number
def temperature(z, rho_over_mean):
    """Adiabatic temperature at redshift z in a region of comoving density
    rho_over_mean (relative to the mean)."""
    return T_i * rho_over_mean ** (2.0 / 3.0) * ((1.0 + z) / (1.0 + z_ref)) ** 2


def jeans_length_comoving(z, rho_over_mean):
    c_s2 = gamma * kB_in_SI * temperature(z, rho_over_mean) / mH_in_kg
    rho_phys = rho_0 * (1.0 + z) ** 3 * rho_over_mean
    return np.sqrt(c_s2 * np.pi / (G_in_SI * rho_phys)) * (1.0 + z)


zfac_end = (1.0 + z_c) / (1.0 + args.z_end)
rho_void_end = 1.0 / (1.0 + zfac_end)
v_infall_c = H_0 * (1.0 + z_c) / np.sqrt(1.0 + z_c) / k_i  # peculiar, at the caustic
c_s_c = np.sqrt(gamma * kB_in_SI * temperature(z_c, 1.0) / mH_in_kg)
print(f"Glass: {args.glass}, {numPart} particles ({numPart_1D}^3), spacing {delta_x / Mpc_in_m:.3f} Mpc/h")
print(f"Start at z = {z_i:g} (zfac = {zfac:.3f}), T = {T_start:.3g} K at mean density "
      f"(T_i = {T_i:g} K at z = {z_ref:g})")
print(f"  comoving Jeans length / spacing: {jeans_length_comoving(z_i, 1.0) / delta_x:.1e} at the start, "
      f"{jeans_length_comoving(args.z_end, rho_void_end) / delta_x:.1e} in the void at z = {args.z_end:g}")
print(f"  linear growth factor of the glass noise to z = {args.z_end:g}: {(1.0 + z_i) / (1.0 + args.z_end):.1f}")
print(f"  infall velocity at the caustic / sound speed: {v_infall_c / c_s_c:.0f}")
print(f"  run with Cosmology:a_begin = {1.0 / (1.0 + z_i):.6g}")

# ---------------------------------------------------
# Particles
q = glass[:, 0] * boxSize + x_min  # Lagrangian x in [x_min, x_max)
coords = np.empty((numPart, 3))
coords[:, 0] = np.mod(q - zfac * np.sin(k_i * q) / k_i - x_min, boxSize)
coords[:, 1] = glass[:, 1] * boxSize
coords[:, 2] = glass[:, 2] * boxSize

v = np.zeros((numPart, 3))
v[:, 0] = -H_0 * (1.0 + z_c) / np.sqrt(1.0 + z_i) * np.sin(k_i * q) / k_i

T = T_start * (1.0 / (1.0 - zfac * np.cos(k_i * q))) ** (2.0 / 3.0)
u = kB_in_SI * T / (gamma - 1.0) / mH_in_kg
h = np.full(numPart, 1.2348 * delta_x)
m = np.full(numPart, m_i)
ids = np.arange(1, numPart + 1)

# Unit conversion
coords /= unit_l_in_si
v /= unit_v_in_si
m /= unit_m_in_si
h /= unit_l_in_si
u /= unit_u_in_si
boxSize /= unit_l_in_si

# File
with h5py.File(args.output, "w") as file:
    grp = file.create_group("/Header")
    grp.attrs["BoxSize"] = [boxSize, boxSize, boxSize]
    grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFilesPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = 0
    grp.attrs["Dimension"] = 3
    grp.attrs["Zeldovich_z_start"] = z_i

    grp = file.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 100.0 * unit_l_in_si
    grp.attrs["Unit mass in cgs (U_M)"] = 1000.0 * unit_m_in_si
    grp.attrs["Unit time in cgs (U_t)"] = 1.0 * unit_t_in_si
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    grp = file.create_group("/PartType0")
    grp.create_dataset("Coordinates", data=coords, dtype="d")
    grp.create_dataset("Velocities", data=v, dtype="f")
    grp.create_dataset("Masses", data=m, dtype="f")
    grp.create_dataset("SmoothingLength", data=h, dtype="f")
    grp.create_dataset("InternalEnergy", data=u, dtype="f")
    grp.create_dataset("ParticleIDs", data=ids, dtype="L")
