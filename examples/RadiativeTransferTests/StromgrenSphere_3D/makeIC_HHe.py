#!/usr/bin/env python3

# ---------------------------------------------------------------------
# Add a single star in the center of a glass distribution
# The gas is set up with a mixture of hydrogen and helium
# with 75% hydrogen and 25% helium in mass.
# This initial condition is comparable to Section 5.3.2 of
# Pawlik & Schaye 2011 doi:10.1111/j.1365-2966.2010.18032.x.
# ---------------------------------------------------------------------

from swiftsimio import Writer
from swiftsimio.units import cosmo_units
import unyt
import numpy as np
import h5py

gamma = 5.0 / 3.0

# switch to replace the central gas particle with a star
# else put the star particle among gas particles
replace_gas = False


def get_number_densities(Temp, XH, XHe):
    """
    Compute the number densities of all species at a given
    temperature following Katz, Hernquist, and Weinberg 1996

    Temp: temperature [unyt quantity]
    XH: total mass fraction of all hydrogen species (HI + HII)
    XHe: total mass fraction of all helium species (HeI + HeII + HeIII)
    """

    # n_H = X_H * rho_gas / m_H =
    # n_He = X_He * rho_gas / m_He = (1 - X_H) / (4 X_H) * n_H
    #      =  X_He / 4(1 - X_He) * nH = y * nH

    if XH == 0:
        nH = 0.0
        nHe = 1.0
    else:
        nH = XH
        nHe = XHe / 4

    # NOTE: This is not the ionization threshold!
    if Temp < 5000 * unyt.K:
        nH0 = nH
        nHp = 0.0
        nHe0 = nHe
        nHep = 0.0
        nHepp = 0.0

    else:

        Temp.convert_to_cgs()
        T = Temp.v
        # Recombination rate for H+ in units of cm^3 s^-1
        A_Hp = (
            8.40e-11
            / np.sqrt(T)
            * (T * 1e-3) ** (-0.2)
            * 1.0
            / (1.0 + (T * 1e-6) ** 0.7)
        )

        # Dielectronic recombination rate for He+ in units of cm^3 s^-1
        A_d = (
            1.9e-3
            / T ** 1.5
            * np.exp(-470000.0 / T)
            * (1.0 + 0.3 * np.exp(-94000.0 / T))
        )
        # Recombination rate for He+ in units of cm^3 s^-1
        A_Hep = 1.5e-10 / T ** 0.6353
        # Recombination rate for He++ in units of cm^3 s^-1
        A_Hepp = (
            3.36e-10
            / np.sqrt(T)
            * (T * 1e-3) ** (-0.2)
            * 1.0
            / (1.0 + (T * 1e-6) ** 0.7)
        )
        # collisional ionization rate for H0 in units of cm^3 s^-1
        #  G_H0 = 1.17e-10 * np.sqrt(T) * np.exp(-157809.1 / T) * 1. / (1. + np.sqrt(T*1e-5))
        G_H0 = (
            5.85e-11
            * np.sqrt(T)
            * np.exp(-157809.1 / T)
            * 1.0
            / (1.0 + np.sqrt(T * 1e-5))
        )
        # collisional ionization rate for He0 in units of cm^3 s^-1
        G_He0 = (
            2.38e-11
            * np.sqrt(T)
            * np.exp(-285335.4 / T)
            * 1.0
            / (1.0 + np.sqrt(T * 1e-5))
        )
        # collisional ionization rate for He+ in units of cm^3 s^-1
        G_Hep = (
            5.68e-12
            * np.sqrt(T)
            * np.exp(-631515.0 / T)
            * 1.0
            / (1.0 + np.sqrt(T * 1e-5))
        )

        # Katz et al. 1996 eq. 33 - 38
        # Note: We assume all photoionization rates to be zero.
        # Also, here we don't care about the actual number density, i.e.
        # about the volume, since it'll cancel out later when we compute
        # the mass fractions.

        nH0 = nH * A_Hp / (A_Hp + G_H0)
        nHp = nH - nH0
        nHep = nHe / (1.0 + (A_Hep + A_d) / G_He0 + G_Hep / A_Hepp)
        nHe0 = nHep * (A_Hep + A_d) / G_He0
        nHepp = nHep * G_Hep / A_Hepp

    # electron density
    ne = nHp + nHep + 2.0 * nHepp

    return nH0, nHp, nHe0, nHep, nHepp, ne


def get_number_densities_array(Temp, XH, XHe):
    """
    Compute the number densities of all species at a given
    temperature following Katz, Hernquist, and Weinberg 1996

    Temp: temperature [unyt array]
    XH: total mass fraction of all hydrogen species (HI + HII)
    XHe: total mass fraction of all helium species (HeI + HeII + HeIII)
    """

    # n_H = X_H * rho_gas / m_H =
    # n_He = X_He * rho_gas / m_He = (1 - X_H) / (4 X_H) * n_H
    #      =  X_He / 4(1 - X_He) * nH = y * nH

    nH = np.zeros(XH.shape, dtype=np.float64)
    nHe = np.zeros(XH.shape, dtype=np.float64)

    mask = XH == 0
    nH[mask] = 0.0
    nHe[mask] = 1.0

    inv_mask = np.logical_not(mask)
    nH[inv_mask] = XH[inv_mask]
    nHe[inv_mask] = 0.25 * XHe[inv_mask]

    nH0 = np.zeros(XH.shape, dtype=np.float64)
    nHp = np.zeros(XH.shape, dtype=np.float64)
    nHe0 = np.zeros(XH.shape, dtype=np.float64)
    nHep = np.zeros(XH.shape, dtype=np.float64)
    nHepp = np.zeros(XH.shape, dtype=np.float64)

    # NOTE: This is not the ionization threshold!
    neutral = Temp < 5000 * unyt.K

    nH0[neutral] = nH[neutral]
    nHp[neutral] = 0.0
    nHe0[neutral] = nHe[neutral]
    nHep[neutral] = 0.0
    nHepp[neutral] = 0.0

    Temp.convert_to_cgs()
    T = Temp.v
    ionized = np.logical_not(neutral)

    # Recombination rate for H+ in units of cm^3 s^-1
    A_Hp = (
        8.40e-11 / np.sqrt(T) * (T * 1e-3) ** (-0.2) * 1.0 / (1.0 + (T * 1e-6) ** 0.7)
    )
    # Dielectronic recombination rate for He+ in units of cm^3 s^-1
    A_d = 1.9e-3 / T ** 1.5 * np.exp(-470000.0 / T) * (1.0 + 0.3 * np.exp(-94000.0 / T))
    # Recombination rate for He+ in units of cm^3 s^-1
    A_Hep = 1.5e-10 / T ** 0.6353
    # Recombination rate for He++ in units of cm^3 s^-1
    A_Hepp = (
        3.36e-10 / np.sqrt(T) * (T * 1e-3) ** (-0.2) * 1.0 / (1.0 + (T * 1e-6) ** 0.7)
    )
    # collisional ionization rate for H0 in units of cm^3 s^-1
    G_H0 = (
        5.85e-11 * np.sqrt(T) * np.exp(-157809.1 / T) * 1.0 / (1.0 + np.sqrt(T * 1e-5))
    )
    # collisional ionization rate for He0 in units of cm^3 s^-1
    G_He0 = (
        2.38e-11 * np.sqrt(T) * np.exp(-285335.4 / T) * 1.0 / (1.0 + np.sqrt(T * 1e-5))
    )
    # collisional ionization rate for He+ in units of cm^3 s^-1
    G_Hep = (
        5.68e-12 * np.sqrt(T) * np.exp(-631515.0 / T) * 1.0 / (1.0 + np.sqrt(T * 1e-5))
    )

    # Katz et al. 1996 eq. 33 - 38
    # Note: We assume all photoionization rates to be zero.
    # Also, here we don't care about the actual number density, i.e.
    # about the volume, since it'll cancel out later when we compute
    # the mass fractions.

    nH0[ionized] = nH[ionized] * A_Hp[ionized] / (A_Hp[ionized] + G_H0[ionized])
    nHp[ionized] = nH[ionized] - nH0[ionized]
    nHep[ionized] = nHe[ionized] / (
        1.0
        + (A_Hep[ionized] + A_d[ionized]) / G_He0[ionized]
        + G_Hep[ionized] / A_Hepp[ionized]
    )
    nHe0[ionized] = nHep[ionized] * (A_Hep[ionized] + A_d[ionized]) / G_He0[ionized]
    nHepp[ionized] = nHep[ionized] * G_Hep[ionized] / A_Hepp[ionized]

    # electron density
    ne = nHp + nHep + 2.0 * nHepp

    return nH0, nHp, nHe0, nHep, nHepp, ne


def get_mass_fractions(T, XH, XHe):
    """
    Compute the mass fractions of all species at a
    given temperature

    T: temperature
    XH: total mass fraction of all hydrogen species (HI + HII)
    XHe: total mass fraction of all helium species (HeI + HeII + HeIII)
    """

    # first get number densities
    if isinstance(XH, np.ndarray):
        nH0, nHp, nHe0, nHep, nHepp, ne = get_number_densities_array(T, XH, XHe)
    else:
        nH0, nHp, nHe0, nHep, nHepp, ne = get_number_densities(T, XH, XHe)

    # now get mass denities in units of atomic mass units
    mH0 = nH0
    mHp = nHp
    mHe0 = 4.0 * nHe0
    mHep = 4.0 * nHep
    mHepp = 4.0 * nHepp
    # neglect electron mass contributions

    mtot = mH0 + mHp + mHe0 + mHep + mHepp  # + me

    XH0 = mH0 / mtot
    XHp = mHp / mtot
    XHe0 = mHe0 / mtot
    XHep = mHep / mtot
    XHepp = mHepp / mtot
    #  Xe = me / mtot

    return XH0, XHp, XHe0, XHep, XHepp


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


if __name__ == "__main__":

    glass = h5py.File("glassCube_64.hdf5", "r")
    parts = glass["PartType0"]
    xp = parts["Coordinates"][:]
    h = parts["SmoothingLength"][:]
    glass.close()

    r = np.sqrt(np.sum((0.5 - xp) ** 2, axis=1))

    if replace_gas == True:
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
    XHI, XHII, XHeI, XHeII, XHeIII = get_mass_fractions(T, XH, XHe)
    mu = mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
    internal_energy = internal_energy(T, mu)

    w.gas.internal_energy = np.ones(xp.shape[0], dtype=np.float64) * internal_energy

    w.write("stromgrenSphere-3D-HHe.hdf5")
