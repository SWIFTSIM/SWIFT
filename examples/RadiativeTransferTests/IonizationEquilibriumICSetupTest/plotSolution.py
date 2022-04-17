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


import numpy as np
from matplotlib import pyplot as plt
import unyt
import swiftsimio

plotkwargs = {"alpha": 0.5}


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


data = swiftsimio.load("output_0000.hdf5")
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
gamma = data.gas.metadata.gas_gamma[0]

# Given the particle internal energies, re-compute the
# results to cross-compare
T = gas_temperature(u, mu, gamma).to("K")
Tmin = T.v.min()
Tmax = T.v.max()


u_theory, mu_theory, T_theory, XHI_theory, XHII_theory, XHeI_theory, XHeII_theory, XHeIII_theory = np.loadtxt(
    "IonizationEquilibriumICSetupTestReference.txt", dtype=np.float64, unpack=True
)

# ------------------
# Plot figures
# ------------------

fig = plt.figure(figsize=(12, 6), dpi=300)
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

ax1.loglog(T, u, label="obtained results")
ax1.loglog(T_theory, u_theory, label="expected results", ls=":")
ax1.set_xlabel("gas temperature [K]")
ax1.set_ylabel("specific internal energy [erg/g]")
ax1.legend()
ax1.grid()

ax2.semilogx(T, mu, label="obtained results")
ax2.semilogx(T_theory, mu_theory, label="expected results", ls=":")
ax2.set_xlabel("gas temperature [K]")
ax2.set_ylabel("mean molecular weight [1]")
ax2.legend()
ax2.grid()

ax3.set_title("obtained results")
ax3.semilogx(u, XHI + XHII + XHeI + XHeII + XHeIII, "k", label="total", ls="-")

ax3.semilogx(u, XHI, label="HI", ls=":", **plotkwargs, zorder=1)
ax3.semilogx(u, XHII, label="HII", ls="-.", **plotkwargs)
ax3.semilogx(u, XHeI, label="HeI", ls=":", **plotkwargs)
ax3.semilogx(u, XHeII, label="HeII", ls="-.", **plotkwargs)
ax3.semilogx(u, XHeIII, label="HeIII", ls="--", **plotkwargs)
ax3.legend()
ax3.set_xlabel("specific internal energy [erg/g]")
ax3.set_ylabel("gas mass fractions [1]")
ax3.grid()


ax4.semilogx(
    u_theory,
    XHI_theory + XHII_theory + XHeI_theory + XHeII_theory + XHeIII_theory,
    "k",
    label="total",
    ls="-",
)

ax4.semilogx(u, XHI, label="obtained results", zorder=-1, color="k", lw=3, alpha=0.2)
ax4.semilogx(u, XHII, zorder=-1, color="k", lw=3, alpha=0.2)
ax4.semilogx(u, XHeI, zorder=-1, color="k", lw=3, alpha=0.2)
ax4.semilogx(u, XHeII, zorder=-1, color="k", lw=3, alpha=0.2)
ax4.semilogx(u, XHeIII, zorder=-1, color="k", lw=3, alpha=0.2)

ax4.semilogx(u_theory, XHI_theory, label="expected HI", ls=":", **plotkwargs)
ax4.semilogx(u_theory, XHII_theory, label="expected HII", ls="-.", **plotkwargs)
ax4.semilogx(u_theory, XHeI_theory, label="expected HeI", ls=":", **plotkwargs)
ax4.semilogx(u_theory, XHeII_theory, label="expected HeII", ls="-.", **plotkwargs)
ax4.semilogx(u_theory, XHeIII_theory, label="expected HeIII", ls="--", **plotkwargs)
ax4.legend()
ax4.grid()
#  ax4.set_xlabel("gas temperature [K]")
ax4.set_xlabel("specific internal energy [erg/g]")
ax4.set_ylabel("gas mass fractions [1]")
ax4.set_title("Expected results")

plt.tight_layout()


plt.savefig("ionization_equilibrium_test.png")
