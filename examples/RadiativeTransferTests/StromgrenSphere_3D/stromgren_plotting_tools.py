import os
import numpy as np
import unyt
import copy


mamu = {"e": 0.0, "HI": 1.0, "HII": 1.0, "HeI": 4.0, "HeII": 4.0, "HeIII": 4.0}


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
    return sA

    # clean string up


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
