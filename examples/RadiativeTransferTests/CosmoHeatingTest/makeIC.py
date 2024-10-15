#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
#               2024 Stan Verhoeve (s06verhoeve@gmail.com)
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


# -----------------------------------------------------------
# Use 10 particles in each dimension to generate a uniform box
# with a fixed amount of radiation.
# -----------------------------------------------------------

from swiftsimio import Writer
from swiftsimio.units import cosmo_units

import unyt
import numpy as np
import h5py
import yaml

# Grab scale factor at start of run
with open(r"rt_heating_test.yml") as paramfile:
    params = yaml.load(paramfile, Loader=yaml.FullLoader)
    a_begin = params["Cosmology"]["a_begin"]
    a_begin = float(a_begin)

# Reduced speed of light
f_reduce_speed_of_light = 1.0e-6

# number of particles in each dimension
n_p = 10
nparts = n_p ** 3

# filename of ICs to be generated
outputfilename = "heating_test.hdf5"
# adiabatic index
gamma = 5.0 / 3.0
# total hydrogen mass fraction
XH = 0.76
# total helium mass fraction
XHe = 0.24
# boxsize
boxsize = 1 * unyt.kpc
# initial gas temperature
initial_temperature = 1e3 * unyt.K
# particle mass
pmass = (unyt.atomic_mass_unit / unyt.cm ** 3) * (boxsize ** 3 / nparts) * a_begin ** 3
pmass = pmass.to("Msun")

# -----------------------------------------------


def internal_energy(T, mu):
    """
    Compute the internal energy of the gas for a given
    temperature and mean molecular weight
    """
    # Using u = 1 / (gamma - 1) * p / rho
    #   and p = N/V * kT = rho / (mu * m_u) * kT

    u = unyt.boltzmann_constant * T / (gamma - 1) / (mu * unyt.atomic_mass_unit)
    return u


def mean_molecular_weight(XH0, XHp, XHe0, XHep, XHepp):
    """
    Determines the mean molecular weight for given 
    mass fractions of
        hydrogen:   XH0
        H+:         XHp
        He:         XHe0
        He+:        XHep
        He++:       XHepp

    returns:
        mu: mean molecular weight [in atomic mass units]
        NOTE: to get the actual mean mass, you still need
        to multiply it by m_u, as is tradition in the formulae
    """

    # 1/mu = sum_j X_j / A_j * (1 + E_j)
    # A_H    = 1, E_H    = 0
    # A_Hp   = 1, E_Hp   = 1
    # A_He   = 4, E_He   = 0
    # A_Hep  = 4, E_Hep  = 1
    # A_Hepp = 4, E_Hepp = 2
    one_over_mu = XH0 + 2 * XHp + 0.25 * XHe0 + 0.5 * XHep + 0.75 * XHepp

    return 1.0 / one_over_mu


# assume everything is neutral initially
mu = mean_molecular_weight(XH, 0, XHe, 0.0, 0.0)
u_part = internal_energy(initial_temperature, mu)
pmass = pmass.to("Msun")


xp = unyt.unyt_array(np.zeros((nparts, 3), dtype=np.float32), boxsize.units)
dx = boxsize / n_p
ind = 0
for i in range(n_p):
    x = (i + 0.5) * dx
    for j in range(n_p):
        y = (j + 0.5) * dx
        for k in range(n_p):
            z = (k + 0.5) * dx

            xp[ind] = (x, y, z)
            ind += 1

w = Writer(cosmo_units, boxsize, dimension=3)


w.gas.coordinates = xp
w.gas.velocities = np.zeros(xp.shape, dtype=np.float32) * (unyt.km / unyt.s)
w.gas.masses = np.ones(nparts, dtype=np.float32) * pmass
w.gas.internal_energy = np.ones(nparts, dtype=np.float32) * u_part

# Generate initial guess for smoothing lengths based on MIPS
w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)

# If IDs are not present, this automatically generates
w.write(outputfilename)

# Now open file back up again and add RT data.

F = h5py.File(outputfilename, "r+")
header = F["Header"]
nparts = header.attrs["NumPart_ThisFile"][0]
parts = F["/PartType0"]

# Create initial ionization species mass fractions.
# Assume everything neutral initially
# NOTE: grackle doesn't really like exact zeroes, so
# use something very small instead.
HIdata = np.ones(nparts, dtype=np.float32) * XH
HIIdata = np.ones(nparts, dtype=np.float32) * 1e-12
HeIdata = np.ones(nparts, dtype=np.float32) * XHe
HeIIdata = np.ones(nparts, dtype=np.float32) * 1e-12
HeIIIdata = np.ones(nparts, dtype=np.float32) * 1e-12

parts.create_dataset("MassFractionHI", data=HIdata)
parts.create_dataset("MassFractionHII", data=HIIdata)
parts.create_dataset("MassFractionHeI", data=HeIdata)
parts.create_dataset("MassFractionHeII", data=HeIIdata)
parts.create_dataset("MassFractionHeIII", data=HeIIIdata)


# Add photon groups
nPhotonGroups = 3

# with this IC, the radiative cooling is negligible.
#  photon_energy = u_part * pmass * 5.0
#  photon_energy = np.arange(1, nPhotonGroups+1) * photon_energy

# Fluxes from the Iliev Test0 part3
fluxes_iliev = np.array([1.350e1, 2.779e1, 6.152e0]) * unyt.erg / unyt.s / unyt.cm ** 2
energy_density = fluxes_iliev / unyt.c
photon_energy = energy_density * boxsize ** 3 / nparts * a_begin ** 3
photon_energy = photon_energy * f_reduce_speed_of_light

photon_energy.convert_to_units(cosmo_units["energy"])
photon_fluxes = 0.333333 * unyt.c * photon_energy
photon_fluxes.convert_to_units(
    cosmo_units["energy"] * cosmo_units["length"] / cosmo_units["time"]
)


for grp in range(nPhotonGroups):
    dsetname = "PhotonEnergiesGroup{0:d}".format(grp + 1)
    energydata = np.ones(nparts, dtype=np.float32) * photon_energy[grp]
    parts.create_dataset(dsetname, data=energydata)

    dsetname = "PhotonFluxesGroup{0:d}".format(grp + 1)
    fluxdata = np.zeros((nparts, 3), dtype=np.float32)
    fluxdata[:, 0] = photon_fluxes[grp]
    parts.create_dataset(dsetname, data=fluxdata)

# close up, and we're done!
F.close()
