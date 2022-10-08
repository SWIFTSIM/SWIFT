#!/usr/bin/env python3

# ---------------------------------------------------------------------
# Add a single star in the center of a glass distribution
# The gas is set up with a mixture of hydrogen and helium
# with 75% hydrogen and 25% helium in mass.
# This initial condition is comparable to Section 5.3.2 of
# Pawlik & Schaye 2011 doi:10.1111/j.1365-2966.2010.18032.x.
# ---------------------------------------------------------------------

import h5py
import numpy as np
import unyt
from swiftsimio import Writer
from swiftsimio.units import cosmo_units

import stromgren_plotting_tools as spt

gamma = 5.0 / 3.0

# switch to replace the central gas particle with a star
# else put the star particle among gas particles
replace_gas = True


if __name__ == "__main__":

    glass = h5py.File("glassCube_64.hdf5", "r")
    parts = glass["PartType0"]
    xp = parts["Coordinates"][:]
    h = parts["SmoothingLength"][:]
    glass.close()

    r = np.sqrt(np.sum((0.5 - xp) ** 2, axis=1))

    if replace_gas:
        # replace a central gas particle with a star particle
        rmin = np.argmin(r)
        xs = xp[rmin]
        xp = np.delete(xp, rmin, axis=0)
        h = np.delete(h, rmin)
    else:
        # find particles closest to the center
        # and select a couple of them to put the star in their middle
        mininds = np.argsort(r)
        center_parts = xp[mininds[:4]]
        xs = center_parts.sum(axis=0) / center_parts.shape[0]

    # Double-check all particles for boundaries
    for i in range(3):
        mask = xp[:, i] < 0.0
        xp[mask, i] += 1.0
        mask = xp[:, i] > 1.0
        xp[mask, i] -= 1.0

    # Set up metadata
    unitL = unyt.Mpc
    edgelen = 22 * 1e-3 * unitL  # 22 so we can cut off 1kpc on each edge for image
    edgelen = edgelen.to(unitL)
    boxsize = np.array([1.0, 1.0, 1.0]) * edgelen

    xs = unyt.unyt_array(
        [np.array([xs[0] * edgelen, xs[1] * edgelen, xs[2] * edgelen])], unitL
    )
    xp *= edgelen
    h *= edgelen

    w = Writer(unit_system=cosmo_units, box_size=boxsize, dimension=3)

    # write particle positions and smoothing lengths
    w.gas.coordinates = xp
    w.stars.coordinates = xs
    w.gas.velocities = np.zeros(xp.shape) * (unitL / unyt.Myr)
    w.stars.velocities = np.zeros(xs.shape) * (unitL / unyt.Myr)
    w.gas.smoothing_length = h
    w.stars.smoothing_length = w.gas.smoothing_length[:1]

    # get gas masses
    XH = 0.75  # hydrogen mass fraction
    XHe = 0.25  # helium mass fraction
    nH = 1e-3 * unyt.cm ** (-3)
    rho_gas = nH * unyt.proton_mass / XH
    Mtot = rho_gas * edgelen ** 3
    mpart = Mtot / xp.shape[0]
    mpart = mpart.to(cosmo_units["mass"])
    w.gas.masses = np.ones(xp.shape[0], dtype=np.float64) * mpart
    w.stars.masses = np.ones(xs.shape[0], dtype=np.float64) * mpart

    # get gas internal energy for a given temperature and composition
    T = 100.0 * unyt.K
    XHI, XHII, XHeI, XHeII, XHeIII = spt.get_mass_fractions(T, XH, XHe)
    mu = spt.mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
    internal_energy = spt.internal_energy(T, mu, gamma)

    w.gas.internal_energy = np.ones(xp.shape[0], dtype=np.float64) * internal_energy

    w.write("stromgrenSphere-3D-HHe.hdf5")
