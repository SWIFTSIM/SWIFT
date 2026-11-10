###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
import unyt
import swiftsimio as sw
import numpy as np

# Some constants
solar_mass_cgs = 1.988410e33
kpc_in_cm = 3.08567758149e21
mp_cgs = 1.67e-24
boltzmann_k_cgs = 1.38e-16

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
rho_cgs = mp_cgs  # Background density
u0_cgs = 1.2e12  # Desired initial internal energy (1.2e12 ~ 10^4K)
P_cgs = rho_cgs * u0_cgs * (gamma - 1.0)  # Background pressure
fileName = "bh_swallowing.hdf5"

# Units
unit_l_cgs = 3.08567758149e24  # kpc
unit_m_cgs = 1.988410e43  # 10^10 Msun
unit_v_cgs = 1e5  # km / s
unit_A_cgs = 1.0
unit_T_cgs = 1.0
unit_t_cgs = unit_l_cgs / unit_v_cgs

boxsize_cgs = 10.0 * kpc_in_cm
vol_cgs = boxsize_cgs**3

# ---------------------------------------------------
glass = sw.load("glassCube_32.hdf5")

# Read particle positions and h from the glass
pos = glass.gas.coordinates.value
eps = 1e-6
pos = (pos - pos.min()) / (pos.max() - pos.min() + eps) * boxsize_cgs / unit_l_cgs
h = glass.gas.smoothing_lengths.value * boxsize_cgs / unit_l_cgs

numPart = np.size(h)

# Generate extra arrays
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m_cgs = np.zeros(numPart)
u_cgs = np.zeros(numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)

m_cgs[:] = rho_cgs * vol_cgs / numPart
u_cgs[:] = P_cgs / (rho_cgs * (gamma - 1))

# BHs
bh_pos = np.zeros((2, 3))
bh_pos[0, :] = 0.5 * boxsize_cgs / unit_l_cgs

# diff = [4.812127e-03, 4.908179e-03, -4.878537e-03] - bh_pos[0,:]

# print diff

bh_pos[1, :] = 0.5 * boxsize_cgs / unit_l_cgs + np.diff

print(bh_pos)

# BHs don't move
bh_v = np.zeros((2, 3))
bh_v[:, :] = 0.0

# increase mass to keep it at center
bh_m_cgs = np.ones(2) * m_cgs[0]
bh_ids = np.linspace(numPart + 1, numPart + 2, 2, dtype="L")
bh_h = np.ones(2) * [h.max()]

print(bh_ids)

# --------------------------------------------------

# Check quantities are correct for debugging
print(("boxsize kpc " + str(boxsize_cgs / kpc_in_cm)))
print(("density cm^-3 " + str(rho_cgs / mp_cgs)))
print(
    (
        "initial temperature K "
        + str(u_cgs[0] / boltzmann_k_cgs * ((gamma - 1) * rho_cgs))
    )
)

# Convert to internal units
bh_m = bh_m_cgs / unit_m_cgs
m[:] = m_cgs / unit_m_cgs
u[:] = u_cgs * unit_v_cgs**-2
boxsize = boxsize_cgs / unit_l_cgs


# --------------------------------------------------

# File
UL = unyt.Unit(unit_l_cgs * unyt.cm)
UM = unyt.Unit(unit_m_cgs * unyt.g)
Ut = unyt.Unit(unit_t_cgs * unyt.s)
units = unyt.UnitSystem("bh_swallowing", UL, UM, Ut)

w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [boxsize] * 3, UL, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

# Particle group
w.gas.coordinates = pos * UL * w.a.comoving
w.gas.velocities = v * UL / Ut * w.a.comoving**0
w.gas.masses = m * UM * w.a.comoving**0
w.gas.smoothing_lengths = h * UL * w.a.comoving
w.gas.internal_energy = u * UL**2 / Ut**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

# black hole group
w.black_holes.coordinates = bh_pos * UL * w.a.comoving
w.black_holes.velocities = bh_v * UL / Ut * w.a.comoving**0
w.black_holes.masses = bh_m * UM * w.a.comoving**0
w.black_holes.particle_ids = bh_ids[::-1].astype(np.uint64)

w.write(fileName)

# Add the extra fields
with h5py.File(fileName, "r+") as f:
    f["PartType5"].create_dataset("SmoothingLength", data=bh_h, dtype="f")
