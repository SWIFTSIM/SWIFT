###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Maarten Elion (elion@lorentz.leidenuniv.nl)
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
from numpy import *

# Generates a swift IC file for the Balsaro-Kim test in a periodic cubic box

# some constants
gamma = 5.0 / 3.0  # Gas adiabatic index
mp = 1.67e-24  # g
pc = 3.086e18  # cm

# unit system
uL = 1e3 * pc  # cm, 1 kpc
uM = 6.757e40  # g
ut = 2.552e15  # s
uv = 12.09e5  # cm / s
uU = uv**2  # (cm / s)^2
uA = 1e10

# paramters
rho0 = 1  # Background density
P0 = 0.3  # Baground pressure
L = 0.2  # boxsize

B0 = 3.64731873  # microgauss
B0 /= 1e7 * uM / (ut**2 * uA)  # code units


fileName = "BalsaraKim.hdf5"

# ---------------------------------------------------
glass = sw.load("glassCube_64.hdf5")

# Read particle positions and h from the glass
pos = glass.gas.coordinates.value
h = glass.gas.smoothing_lengths.value

pos *= L
h *= L

for i in range(3):
    pos[pos[:, i] < 0.0, i] += L
    pos[pos[:, i] >= L, i] -= L

numPart = size(h)

print("smoothing lengths (all,min,max):", h, h.min(), h.max())

numPart = size(h)
vol = L**3

# Generate extra arrays
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
u = zeros(numPart)
r = zeros(numPart)
B = zeros((numPart, 3))
A = zeros((numPart, 3))

m[:] = rho0 * vol / numPart
u[:] = P0 / ((gamma - 1.0) * rho0)
B[:, 0] = B0


print("particle mass : %.3f" % (m[0]))
print(u.mean())

A[:, 1] = -B0 * (pos[:, 2] - 0.5 * L) / 2
A[:, 2] = B0 * (pos[:, 1] - 0.5 * L) / 2

print("L:", L)
print("m:", m.mean())
print("h:", h.mean())
print("u:", u.mean(), u.min(), u.max())
print("pos:", pos.min(), pos.max())

# --------------------------------------------------

# File
UL = unyt.Unit(uL * unyt.cm)
UM = unyt.Unit(uM * unyt.g)
Ut = unyt.Unit(ut * unyt.s)
units = unyt.UnitSystem("BalsaraKim", UL, UM, Ut)

w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [L, L, L], UL, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

w.gas.coordinates = pos * UL * w.a.comoving
w.gas.velocities = v * UL / Ut * w.a.comoving**0
w.gas.masses = m * UM * w.a.comoving**0
w.gas.smoothing_lengths = h * UL * w.a.comoving
w.gas.internal_energy = u * UL**2 / Ut**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(uint64)

w.write(fileName)

# Add the magnetic fields and the current unit (not handled by the Writer)
with h5py.File(fileName, "r+") as f:
    f["Units"].attrs["Unit current in cgs (U_I)"] = uA
    grp = f["PartType0"]
    grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
    grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")
