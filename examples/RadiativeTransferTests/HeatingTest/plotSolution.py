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

#-------------------------------------------
# Plot the gas temperature, mean molecular
# weight, and mass fractions
#-------------------------------------------

import numpy as np
from matplotlib import pyplot as plt
import unyt
import swiftsimio
import os

plotkwargs = {"alpha": 0.5}
referenceplotkwargs = {"color":"grey", "lw":4, "alpha":0.6}
snapshot_base = "output"

# -----------------------------------------------------------------------

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
        raise FileNotFoundError("Didn't find any snapshots with basename '"+snapshot_basename+"'")

    snaplist = sorted(snaplist)

    return snaplist


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
    ngroups = int(firstdata.metadata.subgrid_scheme["PhotonGroupNumber"])

    times = np.zeros(nsnaps) * unyt.Myr
    temperatures = np.zeros(nsnaps) * unyt.K
    mean_molecular_weights = np.zeros(nsnaps) * unyt.atomic_mass_unit
    mass_fractions = np.zeros((nsnaps, 5))
    internal_energies = np.zeros(nsnaps) * unyt.Msun * unyt.kpc ** 2 / unyt.Myr ** 2 
    photon_energies = np.zeros((ngroups, nsnaps)) * unyt.Msun * unyt.kpc ** 2 / unyt.Myr ** 2

    for i, snap in enumerate(snaplist):

        data = swiftsimio.load(snap)
        gamma = data.gas.metadata.gas_gamma[0]
        time = data.metadata.time
        gas = data.gas

        u = gas.internal_energies.to("erg/g")
        masses = gas.masses.to("Msun")
        XHI = gas.ion_mass_fractions.HI
        XHII = gas.ion_mass_fractions.HII
        XHeI = gas.ion_mass_fractions.HeI
        XHeII = gas.ion_mass_fractions.HeII
        XHeIII = gas.ion_mass_fractions.HeIII
        XH = XHI + XHII
        XHe = XHeI + XHeII + XHeIII
        mu = mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
        T = gas_temperature(u, mu, gamma).to("K")

        print(u[0])
        print(masses[0])
        um = u.to("kpc**2 / Myr**2") * masses
        print(um[0])
        #  um = um.to("1e40 * erg")
        #  print(um[0])


        times[i] = time.to("Myr")
        temperatures[i] = np.mean(T)
        mean_molecular_weights[i] = np.mean(mu)
        internal_energies[i] = np.mean(um)

        mass_fractions[i, 0] = np.mean(XHI)
        mass_fractions[i, 1] = np.mean(XHII)
        mass_fractions[i, 2] = np.mean(XHeI)
        mass_fractions[i, 3] = np.mean(XHeII)
        mass_fractions[i, 4] = np.mean(XHeIII)

        for g in range(ngroups):
            en = getattr(data.gas.photon_energies, "group" + str(g + 1))
            en = en.to("Msun * kpc**2 / Myr**2")
            photon_energies[g, i] = en.sum() / en.shape[0]


    return times, temperatures, mean_molecular_weights, mass_fractions, internal_energies, photon_energies




if __name__ == "__main__":

    # ------------------
    # Plot figures
    # ------------------

    snaplist = get_snapshot_list(snapshot_base)
    t, T, mu, mass_fraction, u, photon_energies = get_snapshot_data(snaplist)
    ngroups = photon_energies.shape[0]

    #  t_ref, dt_ref, T_ref, mu_ref, rho_ref, rhoHI_ref, rhoHII_ref, rhoHeI_ref, rhoHeII_ref, rhoHeIII_ref , rhoe_ref= np.loadtxt(
    #      "RTCoolingTestReference.txt", dtype=np.float64, unpack=True
    #  )
    #  t_ref *= 1e-6 # turn to Myrs
    #  mass_fraction_ref = np.empty((t_ref.shape[0], 5))
    #  mass_fraction_ref[:,0] = rhoHI_ref / rho_ref
    #  mass_fraction_ref[:,1] = rhoHII_ref / rho_ref
    #  mass_fraction_ref[:,2] = rhoHeI_ref / rho_ref
    #  mass_fraction_ref[:,3] = rhoHeII_ref / rho_ref
    #  mass_fraction_ref[:,4] = rhoHeIII_ref / rho_ref

    fig = plt.figure(figsize=(8, 8), dpi=300)
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    #  ax1.semilogy(t_ref[1:], T_ref[1:], label="reference", **referenceplotkwargs)
    ax1.semilogy(t, T, label="obtained results")
    ax1.set_xlabel("time [Myr]")
    ax1.set_ylabel("gas temperature [K]")
    ax1.legend()
    ax1.grid()

    #  ax2.plot(t_ref, mu_ref, label="reference", **referenceplotkwargs)
    ax2.plot(t, mu, label="obtained results")
    ax2.set_xlabel("time [Myr]")
    ax2.set_ylabel("mean molecular weight")
    ax2.legend()
    ax2.grid()

    total_mass_fraction = np.sum(mass_fraction, axis = 1)
    ax3.plot(t, total_mass_fraction, "k", label="total", ls="-")

    #  ax3.plot(t_ref[1:], mass_fraction_ref[1:,0], label="reference", **referenceplotkwargs, zorder=0)
    #  ax3.plot(t_ref[1:], mass_fraction_ref[1:,1], **referenceplotkwargs, zorder=0)
    #  ax3.plot(t_ref[1:], mass_fraction_ref[1:,2], **referenceplotkwargs, zorder=0)
    #  ax3.plot(t_ref[1:], mass_fraction_ref[1:,3], **referenceplotkwargs, zorder=0)
    #  ax3.plot(t_ref[1:], mass_fraction_ref[1:,4], **referenceplotkwargs, zorder=0)

    ax3.plot(t, mass_fraction[:,0], label="HI", ls=":", **plotkwargs, zorder=1)
    ax3.plot(t, mass_fraction[:,1], label="HII", ls="-.", **plotkwargs, zorder=1)
    ax3.plot(t, mass_fraction[:,2], label="HeI", ls=":", **plotkwargs, zorder=1)
    ax3.plot(t, mass_fraction[:,3], label="HeII", ls="-.", **plotkwargs, zorder=1)
    ax3.plot(t, mass_fraction[:,4], label="HeIII", ls="--", **plotkwargs, zorder=1)
    ax3.legend(loc="upper right")
    ax3.set_xlabel("time [Myr]")
    ax3.set_ylabel("gas mass fractions [1]")
    ax3.grid()

    tot_energy = u + np.sum(photon_energies, axis=0)
    ax4.plot(t, tot_energy, label=f"total energy budget", color="k", ls = "--", **plotkwargs)
    for g in range(ngroups):
        ax4.plot(t, photon_energies[g, :], label=f"photon energies group {g+1}", **plotkwargs)
    ax4.plot(t, u, label=f"gas internal energy", **plotkwargs)
    ax4.set_xlabel("time [Myr]")
    ax4.set_ylabel(r"energy budget [$M_\odot$ kpc$^2$ / Myr$^2$]")
    ax4.legend()
    ax4.grid()


    plt.tight_layout()
    #  plt.show()
    plt.savefig("cooling_test.png")
