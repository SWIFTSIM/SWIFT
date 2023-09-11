#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

# -------------------------------------------
# Plot the gas temperature, mean molecular
# weight, and mass fractions
# -------------------------------------------

import copy
import os

import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt

# arguments for plots of results
plotkwargs = {"alpha": 0.5}
# arguments for plot of references
referenceplotkwargs = {"color": "grey", "lw": 4, "alpha": 0.6}
# arguments for legends
legendprops = {"size": 8}
# snapshot basenames
snapshot_base = "output"

# -----------------------------------------------------------------------

energy_units = unyt.Msun * unyt.kpc ** 2 / unyt.kyr ** 2
mass_units = unyt.Msun


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

    return T


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted
    and return their names as list
    """

    snaplist = []

    dirlist = os.listdir()
    for f in dirlist:
        if f.startswith(snapshot_basename) and f.endswith("hdf5"):
            snaplist.append(f)

    if len(snaplist) == 0:
        raise FileNotFoundError(
            "Didn't find any snapshots with basename '" + snapshot_basename + "'"
        )

    snaplist = sorted(snaplist)

    return snaplist


def get_ion_mass_fractions(swiftsimio_loaded_data):
    """
    Returns the ion mass fractions according to
    the used scheme.

    swiftsimio_loaded_data: the swiftsimio.load() object
    """

    data = swiftsimio_loaded_data
    meta = data.metadata
    gas = data.gas
    with_rt = True
    scheme = None
    try:
        scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))
    except KeyError:
        # allow to read in solutions with only cooling, without RT
        with_rt = False

    if with_rt:
        if scheme.startswith("GEAR M1closure"):
            imf = data.gas.ion_mass_fractions
        elif scheme.startswith("SPH M1closure"):
            # atomic mass
            mamu = {
                "e": 0.0,
                "HI": 1.0,
                "HII": 1.0,
                "HeI": 4.0,
                "HeII": 4.0,
                "HeIII": 4.0,
            }
            mass_function_hydrogen = data.gas.rt_element_mass_fractions.hydrogen
            imf = copy.deepcopy(data.gas.rt_species_abundances)
            named_columns = data.gas.rt_species_abundances.named_columns
            for column in named_columns:
                # abundance is in n_X/n_H unit. We convert it to mass fraction by multipling mass fraction of H
                mass_function = (
                    getattr(data.gas.rt_species_abundances, column)
                    * mass_function_hydrogen
                    * mamu[column]
                )
                setattr(imf, column, mass_function)
        else:
            raise ValueError("Unknown scheme", scheme)
    else:
        # try to find solutions for cooling only runs
        imf = {
            "HI": gas.hi[:],
            "HII": gas.hii[:],
            "HeI": gas.he_i[:],
            "HeII": gas.he_ii[:],
            "HeIII": gas.he_iii[:],
        }

    return imf


def get_snapshot_data(snaplist):
    """
    Extract the relevant data from the list of snapshots.

    Returns:
        numpy arrays of:
            time
            temperatures
            mean molecular weights
            mass fractions
    """

    nsnaps = len(snaplist)
    firstdata = swiftsimio.load(snaplist[0])
    with_rt = True
    try:
        ngroups = int(firstdata.metadata.subgrid_scheme["PhotonGroupNumber"])
    except KeyError:
        # allow to read in solutions with only cooling, without RT
        with_rt = False

    times = np.zeros(nsnaps) * unyt.Myr
    temperatures = np.zeros(nsnaps) * unyt.K
    mean_molecular_weights = np.zeros(nsnaps)
    mass_fractions = np.zeros((nsnaps, 5))
    internal_energies = np.zeros(nsnaps) * energy_units
    if with_rt:
        photon_energies = np.zeros((ngroups, nsnaps)) * energy_units
    else:
        photon_energies = None

    for i, snap in enumerate(snaplist):

        data = swiftsimio.load(snap)
        gamma = data.gas.metadata.gas_gamma[0]
        time = data.metadata.time
        gas = data.gas

        u = gas.internal_energies[:].to(energy_units / mass_units)
        masses = gas.masses[:].to(mass_units)
        imf = get_ion_mass_fractions(data)
        mu = mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
        T = gas_temperature(u, mu, gamma).to("K")
        um = u.to(energy_units / mass_units) * masses
        um.to(energy_units)

        times[i] = time.to("Myr")
        temperatures[i] = np.mean(T)
        mean_molecular_weights[i] = np.mean(mu)
        internal_energies[i] = np.mean(um)

        mass_fractions[i, 0] = np.mean(imf.HI)
        mass_fractions[i, 1] = np.mean(imf.HII)
        mass_fractions[i, 2] = np.mean(imf.HeI)
        mass_fractions[i, 3] = np.mean(imf.HeII)
        mass_fractions[i, 4] = np.mean(imf.HeIII)

        if with_rt:
            for g in range(ngroups):
                en = getattr(data.gas.photon_energies, "group" + str(g + 1))
                en = en[:].to(energy_units)
                photon_energies[g, i] = en.sum() / en.shape[0]

    return (
        times,
        temperatures,
        mean_molecular_weights,
        mass_fractions,
        internal_energies,
        photon_energies,
    )


if __name__ == "__main__":

    # ------------------
    # Plot figures
    # ------------------

    snaplist = get_snapshot_list(snapshot_base)
    t, T, mu, mass_fraction, u, photon_energies = get_snapshot_data(snaplist)
    with_rt = photon_energies is not None
    if with_rt:
        ngroups = photon_energies.shape[0]

    refdata = np.loadtxt("RTCoolingTestReference.txt", dtype=np.float64)
    t_ref = refdata[:, 0]
    dt_ref = refdata[:, 1]
    T_ref = refdata[:, 2]
    mu_ref = refdata[:, 3]
    rho_ref = refdata[:, 4]
    rhoHI_ref = refdata[:, 5]
    rhoHII_ref = refdata[:, 6]
    rhoHeI_ref = refdata[:, 7]
    rhoHeII_ref = refdata[:, 8]
    rhoHeIII_ref = refdata[:, 9]
    rhoe_ref = refdata[:, 10]

    t_ref *= 1e-6  # turn to Myrs
    mass_fraction_ref = np.empty((t_ref.shape[0], 5))
    mass_fraction_ref[:, 0] = rhoHI_ref / rho_ref
    mass_fraction_ref[:, 1] = rhoHII_ref / rho_ref
    mass_fraction_ref[:, 2] = rhoHeI_ref / rho_ref
    mass_fraction_ref[:, 3] = rhoHeII_ref / rho_ref
    mass_fraction_ref[:, 4] = rhoHeIII_ref / rho_ref

    fig = plt.figure(figsize=(8, 8), dpi=300)
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    ax1.semilogy(t_ref[1:], T_ref[1:], label="reference", **referenceplotkwargs)
    ax1.semilogy(t, T, label="obtained results")
    ax1.set_xlabel("time [Myr]")
    ax1.set_ylabel("gas temperature [K]")
    ax1.legend(prop=legendprops)
    ax1.grid()

    ax2.plot(t_ref, mu_ref, label="reference", **referenceplotkwargs)
    ax2.plot(t, mu, label="obtained results")
    ax2.set_xlabel("time [Myr]")
    ax2.set_ylabel("mean molecular weight")
    ax2.legend(prop=legendprops)
    ax2.grid()

    total_mass_fraction = np.sum(mass_fraction, axis=1)
    ax3.plot(t, total_mass_fraction, "k", label="total", ls="-")

    ax3.plot(
        t_ref[1:],
        mass_fraction_ref[1:, 0],
        label="reference",
        **referenceplotkwargs,
        zorder=0,
    )
    ax3.plot(t_ref[1:], mass_fraction_ref[1:, 1], **referenceplotkwargs, zorder=0)
    ax3.plot(t_ref[1:], mass_fraction_ref[1:, 2], **referenceplotkwargs, zorder=0)
    ax3.plot(t_ref[1:], mass_fraction_ref[1:, 3], **referenceplotkwargs, zorder=0)
    ax3.plot(t_ref[1:], mass_fraction_ref[1:, 4], **referenceplotkwargs, zorder=0)

    ax3.plot(t, mass_fraction[:, 0], label="HI", ls=":", **plotkwargs, zorder=1)
    ax3.plot(t, mass_fraction[:, 1], label="HII", ls="-.", **plotkwargs, zorder=1)
    ax3.plot(t, mass_fraction[:, 2], label="HeI", ls=":", **plotkwargs, zorder=1)
    ax3.plot(t, mass_fraction[:, 3], label="HeII", ls="-.", **plotkwargs, zorder=1)
    ax3.plot(t, mass_fraction[:, 4], label="HeIII", ls="--", **plotkwargs, zorder=1)
    ax3.legend(loc="upper right", prop=legendprops)
    ax3.set_xlabel("time [Myr]")
    ax3.set_ylabel("gas mass fractions [1]")
    ax3.grid()

    if with_rt:
        u.convert_to_units(energy_units)
        photon_energies.convert_to_units(energy_units)
        tot_energy = u + np.sum(photon_energies, axis=0)
        ax4.plot(
            t,
            tot_energy,
            label=f"total energy budget",
            color="k",
            ls="--",
            **plotkwargs,
        )
        for g in range(ngroups):
            ax4.plot(
                t,
                photon_energies[g, :],
                label=f"photon energies group {g + 1}",
                **plotkwargs,
            )
        ax4.plot(t, u, label="gas internal energy", **plotkwargs)
        ax4.set_xlabel("time [Myr]")
        ax4.set_ylabel(
            r"energy budget [$" + u.units.latex_representation() + "$]", usetex=True
        )
        ax4.legend(prop=legendprops)
        ax4.grid()

    plt.tight_layout()
    #  plt.show()
    plt.savefig("cooling_test.png")
