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
snapshot_base = "output"
plot_errorbars = True

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
            standard deviations of temperatures
            mean molecular weights
            standard deviations of mean molecular weights
            mass fractions
            standard deviations of mass fractions
    """

    nsnaps = len(snaplist)
    times = np.zeros(nsnaps) * unyt.Myr
    temperatures = np.zeros(nsnaps) * unyt.K
    temperatures_std = np.zeros(nsnaps) * unyt.K
    mean_molecular_weights = np.zeros(nsnaps) * unyt.atomic_mass_unit
    mean_molecular_weights_std = np.zeros(nsnaps) * unyt.atomic_mass_unit
    mass_fractions = np.zeros((nsnaps, 5))
    mass_fractions_std = np.zeros((nsnaps, 5))


    for i, snap in enumerate(snaplist):

        data = swiftsimio.load(snap)
        gas = data.gas

        u = gas.internal_energies.to("erg/g")
        XHI = gas.ion_mass_fractions.HI
        XHII = gas.ion_mass_fractions.HII
        XHeI = gas.ion_mass_fractions.HeI
        XHeII = gas.ion_mass_fractions.HeII
        XHeIII = gas.ion_mass_fractions.HeIII
        XH = XHI + XHII
        XHe = XHeI + XHeII + XHeIII
        mu = mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
        T = gas_temperature(u, mu, gamma).to("K")

        gamma = data.gas.metadata.gas_gamma[0]
        time = data.meta.time


        times[i] = time

        temperatures[i] = np.mean(T)
        temperatures_std[i] = np.std(T)

        mean_molecular_weights[i] = np.mean(mean_molecular_weights)
        mean_molecular_weights_std[i] = np.std(mean_molecular_weights)

        mass_fractions[i, 0] = np.mean(XHI)
        mass_fractions[i, 1] = np.mean(XHII)
        mass_fractions[i, 2] = np.mean(XHeI)
        mass_fractions[i, 3] = np.mean(XHeII)
        mass_fractions[i, 4] = np.mean(XHeIII)

        mass_fractions_std[i, 0] = np.std(XHI)
        mass_fractions_std[i, 1] = np.std(XHII)
        mass_fractions_std[i, 2] = np.std(XHeI)
        mass_fractions_std[i, 3] = np.std(XHeII)
        mass_fractions_std[i, 4] = np.std(XHeIII)
    
    return times, temperatures, temperatures_std, mean_molecular_weights, mean_molecular_weights_std, mass_fractions, mass_fractions_std




#  u_theory, mu_theory, T_theory, XHI_theory, XHII_theory, XHeI_theory, XHeII_theory, XHeIII_theory = np.loadtxt(
#      "IonizationEquilibriumICSetupTestReference.txt", dtype=np.float64, unpack=True
#  )

if __name__ == "__main__":

    # ------------------
    # Plot figures
    # ------------------

    snaplist = get_snapshot_list(snapshot_base)
    t, T, T_std, mu, mu_std, mf, mf_std = get_snapshot_data(snaplist)

    fig = plt.figure(figsize=(12, 6), dpi=300)
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    ax1.loglog(t, T, label="obtained results")
    #  ax1.loglog(T_theory, u_theory, label="expected results", ls=":")
    ax1.set_xlabel("time")
    ax1.set_ylabel("gas temperature [K]")
    ax1.legend()
    ax1.grid()

    ax2.semilogx(t, mu, label="obtained results")
    #  ax2.semilogx(T_theory, mu_theory, label="expected results", ls=":")
    ax2.set_xlabel("time")
    ax2.set_ylabel("mean molecular weight")
    ax2.legend()
    ax2.grid()

    ax3.set_title("obtained results")
    total_mf = np.sum(mf, axis = 1)
    ax3.semilogx(time, total_mf, "k", label="total", ls="-")

    ax3.semilogx( time, mf[:,0], label="HI", ls=":", **plotkwargs, zorder=1)
    ax3.semilogx( time, mf[:,1], label="HII", ls="-.", **plotkwargs)
    ax3.semilogx( time, mf[:,2], label="HeI", ls=":", **plotkwargs)
    ax3.semilogx( time, mf[:,3], label="HeII", ls="-.", **plotkwargs)
    ax3.semilogx( time, mf[:,4], label="HeIII", ls="--", **plotkwargs)
    ax3.legend()
    ax3.set_xlabel("time")
    ax3.set_ylabel("gas mass fractions [1]")
    ax3.grid()


    #  plt.tight_layout()

    plt.show()
    #  plt.savefig("cooling_test.png")
