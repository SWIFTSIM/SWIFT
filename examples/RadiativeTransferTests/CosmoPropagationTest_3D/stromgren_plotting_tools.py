# Contains commonly used functions related to plotting results
# of the Strömgren Sphere Examples and other gas related functions.

import copy
import os

import numpy as np
import unyt

# species masses in atomic mass units
mamu = {"e": 0.0, "HI": 1.0, "HII": 1.0, "HeI": 4.0, "HeII": 4.0, "HeIII": 4.0}


def get_number_densities(Temp, XH, XHe):
    """
    Compute the number densities of all species at a given
    temperature following Katz, Hernquist, and Weinberg 1996

    Temp: temperature [unyt quantity]
    XH: total mass fraction of all hydrogen species (HI + HII)
    XHe: total mass fraction of all helium species (HeI + HeII + HeIII)
    """

    # n_H = X_H * rho_gas / m_H
    # then m_H = rho_gas * X_H / n_H
    # n_He = X_He * rho_gas / m_He =
    #      = (1 - X_H) * rho_gas / (4 * m_H)
    #      = (1 - X_H) * rho_gas / (4 * rho_gas * x_H / n_H)
    #      = (1 - X_H) / (4 * X_H) * n_H
    #      =  X_He / 4(1 - X_He) * n_H = y * n_H

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

    # n_H = X_H * rho_gas / m_H
    # then m_H = rho_gas * X_H / n_H
    # n_He = X_He * rho_gas / m_He =
    #      = (1 - X_H) * rho_gas / (4 * m_H)
    #      = (1 - X_H) * rho_gas / (4 * rho_gas * x_H / n_H)
    #      = (1 - X_H) / (4 * X_H) * n_H
    #      =  X_He / 4(1 - X_He) * n_H = y * n_H

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


def internal_energy(T, mu, gamma):
    """
    Compute the internal energy of the gas for a given
    temperature and mean molecular weight
    """
    # Using u = 1 / (gamma - 1) * p / rho
    #   and p = N/V * kT = rho / (mu * m_u) * kT

    u = unyt.boltzmann_constant * T / (gamma - 1) / (mu * unyt.atomic_mass_unit)
    return u


def get_soundspeed_from_internal_energy(data):
    """
    Compute the local sound speed for all particles.
    data: swiftsimio.load() object.
    """

    u = data.gas.internal_energies
    gamma = data.metadata.gas_gamma

    return np.sqrt(u * gamma / (gamma - 1))


def get_soundspeed_from_density_pressure(data):
    """
    Compute the local sound speed for all particles.
    data: swiftsimio.load() object.
    """

    gamma = data.metadata.gas_gamma
    P = data.gas.pressures
    rho = data.gas.densities

    return np.sqrt(gamma * P / rho)


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


def gas_temperature(u, mu, gamma):
    """
    Compute the gas temperature given the specific internal 
    energy u and the mean molecular weight mu
    """

    # Using u = 1 / (gamma - 1) * p / rho
    #   and p = N/V * kT = rho / (mu * m_u) * kT

    T = u * (gamma - 1) * mu * unyt.atomic_mass_unit / unyt.boltzmann_constant

    return T.to("K")


def get_snapshot_list(snapshot_basename="output", plot_all=True, snapnr=0):
    """
    Find the snapshot(s) that are to be plotted 
    and return their names as list
    """

    snaplist = []

    if plot_all:
        dirlist = os.listdir()
        for f in dirlist:
            if f.startswith(snapshot_basename) and f.endswith("hdf5"):
                snaplist.append(f)

        snaplist = sorted(snaplist)

    else:
        fname = snapshot_basename + "_" + str(snapnr).zfill(4) + ".hdf5"
        if not os.path.exists(fname):
            print("Didn't find file", fname)
            quit(1)
        snaplist.append(fname)

    return snaplist


def get_imf(scheme, data):
    """
    Get the ion mass fraction (imf) according to the scheme.
    return a class with ion mass fraction for species X, 
    including HI, HII, HeI, HeII, HeIII:
    The ion mass function can be accessed through: imf.X
    The unit is in m_X/m_tot, where m_X is the mass in species X
    and m_tot is the total gas mass.
    """
    if scheme.startswith("GEAR M1closure"):
        imf = data.gas.ion_mass_fractions
    elif scheme.startswith("SPH M1closure"):
        # atomic mass
        mass_fraction_hydrogen = data.gas.rt_element_mass_fractions.hydrogen
        imf = copy.deepcopy(data.gas.rt_species_abundances)
        named_columns = data.gas.rt_species_abundances.named_columns
        for column in named_columns:
            # abundance is in n_X/n_H unit. We convert it to mass fraction by multipling mass fraction of H
            mass_fraction = (
                getattr(data.gas.rt_species_abundances, column)
                * mass_fraction_hydrogen
                * mamu[column]
            )
            setattr(imf, column, mass_fraction)
    else:
        raise ValueError("Unknown scheme", scheme)
    return imf


def get_abundances(scheme, data):
    """
    Get the species abundance according to the scheme
    return a class with normalized number densities for abunance X, 
    including HI, HII, HeI, HeII, HeIII:
    The abundances can be accessed through: sA.X
    The unit is in n_X/n_H, where n_X is the number density of species X
    and n_H is the number density of hydrogen.
    """
    if scheme.startswith("GEAR M1closure"):
        # atomic mass
        sA = copy.deepcopy(data.gas.ion_mass_fractions)
        mass_fraction_hydrogen = (
            data.gas.ion_mass_fractions.HI + data.gas.ion_mass_fractions.HII
        )
        # abundance is in n_X/n_H unit. We convert mass fraction to abundance by dividing mass fraction of H
        sA.HI = data.gas.ion_mass_fractions.HI / mass_fraction_hydrogen / mamu["HI"]
        sA.HII = data.gas.ion_mass_fractions.HII / mass_fraction_hydrogen / mamu["HII"]
        sA.HeI = data.gas.ion_mass_fractions.HeI / mass_fraction_hydrogen / mamu["HeI"]
        sA.HeII = (
            data.gas.ion_mass_fractions.HeII / mass_fraction_hydrogen / mamu["HeII"]
        )
        sA.HeIII = (
            data.gas.ion_mass_fractions.HeIII / mass_fraction_hydrogen / mamu["HeIII"]
        )
    elif scheme.startswith("SPH M1closure"):
        sA = data.gas.rt_species_abundances
    else:
        raise ValueError("Unknown scheme", scheme)
    return sA


def trim_paramstr(paramstr):
    """
    clean up strings in the form [x,y,z,...]
    and return an array: array([x,y,z,...])
    """
    paramstr = paramstr.strip()
    if paramstr.startswith("["):
        paramstr = paramstr[1:]
    if paramstr.endswith("]"):
        paramstr = paramstr[:-1]

    params = paramstr.split(",")
    paramtrimmed = []
    for er in params:
        paramtrimmed.append(float(er))
    return paramtrimmed
