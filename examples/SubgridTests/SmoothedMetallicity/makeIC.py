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
import numpy as np
import unyt
import swiftsimio as sw

# Generates a swift IC file for the Smoothed Metallicity test
# in a periodic cubic box

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
rho0 = 1.0  # Background density
P0 = 1e-6  # Background pressure
Nelem = 10  # Gear: 10, EAGLE: 9
low_metal = -6  # Low iron fraction
high_metal = -5.5  # high iron fraction
max_shift = 1  # Shift between the different elements
sigma_metal = 0.2  # relative standard deviation for the metallicities
fileName = "smoothed_metallicity.hdf5"

# shift all metals in order to obtain nicer plots
low_metal = [low_metal] * Nelem + np.linspace(0, max_shift, Nelem)
low_metal = 10**low_metal

high_metal = [high_metal] * Nelem + np.linspace(0, max_shift, Nelem)
high_metal = 10**high_metal

# ---------------------------------------------------
glass = sw.load("glassCube_32.hdf5")

# Read particle positions and h from the glass
pos = glass.gas.coordinates.value
h = glass.gas.smoothing_lengths.value

numPart = np.size(h)
vol = 1.0

# Generate extra arrays
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)
r = np.zeros(numPart)
mass_frac = np.zeros((numPart, Nelem))

m[:] = rho0 * vol / numPart
u[:] = P0 / (rho0 * (gamma - 1))

# set metallicities
select = pos[:, 0] < 0.5
nber = sum(select)
mass_frac[select, :] = low_metal * (
    1 + np.random.normal(loc=0.0, scale=sigma_metal, size=(nber, Nelem))
)
nber = numPart - nber
mass_frac[~select, :] = high_metal * (
    1 + np.random.normal(loc=0.0, scale=sigma_metal, size=(nber, Nelem))
)

v[select, 2] = 1
v[~select, 2] = -1

# --------------------------------------------------

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [1.0, 1.0, 1.0], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

# Particle group
w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileName)

# Add the extra fields
with h5py.File(fileName, "r+") as f:
    grp = f["PartType0"]
    grp.create_dataset("MetalMassFraction", data=mass_frac, dtype="d")
    grp.create_dataset("ElementAbundance", data=mass_frac, dtype="d")
