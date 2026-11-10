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
import numpy as np

# Generates a swift IC file for the Balsaro-Kim test in a periodic cubic box

# some constants
gamma = 5.0 / 3.0  # Gas adiabatic index
mp = 1.67e-24  # g
pc = 3.086e18  # cm
Msun = 1.9891e33  # g
Myr = 1e6 * 60 * 60 * 24 * 365  # s

# unit system
uL = 1e3 * pc  # cm, 1 kpc
uM = 6.757e40  # g, 3.398×10^8 M_sun
ut = 2.552e15  # s, 2.5577×10^7  years
uv = 12.1e5  # cm / s, 38.23 km/s
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

numPart = np.size(h)

print("smoothing lengths (all,min,max):", h, h.min(), h.max())

numPart = np.size(h)
vol = L**3

# Generate extra arrays
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)
r = np.zeros(numPart)
B = np.zeros((numPart, 3))
A = np.zeros((numPart, 3))

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

### Initialise star particles ###

injection_data = np.loadtxt("injection_coords.txt", delimiter=",", skiprows=1)

numStars = injection_data.shape[0]

pos_stars = injection_data[:, 1:4] * L  # positions in the file are relative to boxsize
v_stars = np.zeros((numStars, 3))
ids_stars = np.linspace(numPart + 1, numPart + numStars, numStars)
m_stars = np.zeros(numStars)
h_stars = np.zeros(numStars)
z_stars = np.zeros(numStars)
n0_stars = np.zeros(numStars)

m_stars[:] = (
    57.6 * Msun / uM
)  # corresponds to an available energy of ~E_sn per star part
birthtimes = injection_data[:, 0] * Myr / ut
z_stars[:] = 0.1 * 0.0127
n0_stars[:] = 0.67 * uL**3
h_stars[:] = 5 * pc / uL

# --------------------------------------------------

# File
UL = unyt.Unit(uL * unyt.cm)
UM = unyt.Unit(uM * unyt.g)
Ut = unyt.Unit(ut * unyt.s)
units = unyt.UnitSystem("balsara_kim", UL, UM, Ut)

w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [L, L, L], UL, comoving=True, scale_factor=1.0, scale_exponent=1
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

# Stars group
w.stars.coordinates = pos_stars * UL * w.a.comoving
w.stars.velocities = v_stars * UL / Ut * w.a.comoving**0
w.stars.masses = m_stars * UM * w.a.comoving**0
w.stars.smoothing_lengths = h_stars * UL * w.a.comoving
w.stars.particle_ids = ids_stars.astype(np.uint64)

w.write(fileName)

# Add the extra fields
with h5py.File(fileName, "r+") as f:
    # The Writer always writes a current unit of 1; the B fields are in uA
    f["Units"].attrs["Unit current in cgs (U_I)"] = [uA]

    grp = f["PartType0"]
    grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
    grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")

    grp = f["PartType4"]
    grp.create_dataset("StellarFormationTime", data=birthtimes, dtype="f")
