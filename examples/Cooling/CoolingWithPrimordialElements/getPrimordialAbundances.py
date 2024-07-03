#!/usr/bin/env python3

"""
This script is an example of how to compute the
fractional abundances of different elements with
respect to Hydrogen.

For comparison, see
for example Dalgarno 2005
"""

# some parameters (Planck 2018)
Omega0 = 0.315
Omega_b = 0.049
h0 = 0.674
T0 = 2.73
fH = 0.76
z0 = 300

z = 77

# Temperature
T = T0 / (1 + 200) * (1 + z) ** 2  # Anninos & Norman 1996 eq.11

# all abundances per unit of nH
nHII = 1.2e-5 * Omega0 ** 0.5 / (h0 * Omega_b)  # Anninos & Norman 1996 eq. 9

nHeI = (1 - fH) / fH / 4
nHeII = 0  # assume neutral He
nHeIII = 0  # assume neutral He

nHM = 2e-11 * (1 + z) ** (1.76) * nHII  # Anninos & Norman 1996 eq. 12

nH2I = (
    2e-20 * fH * Omega0 ** (3 / 2.0) / (h0 * Omega_b) * (1 + z0) ** (5.1)
)  # Anninos & Norman 1996 eq. 15
nH2II = 1e-17 * (1 + z) ** (3.6) * nHII  # Anninos & Norman 1996 eq. 13

nDI = 4e-5  # Bromm 20002 from Galli & Palla 1998
nDII = 1e-8  # Bromm 20002

nHDI = 1e-3 * nH2I  # Bromm 20002


print("T                           : %g K" % T)

print("initial_nHII_to_nH_ratio    : %g" % nHII)

print("initial_nHeI_to_nH_ratio    : %g" % nHeI)
print("initial_nHeII_to_nH_ratio   : %g" % nHeII)
print("initial_nHeIII_to_nH_ratio  : %g" % nHeIII)

print("initial_nDI_to_nH_ratio     : %g" % nDI)
print("initial_nDII_to_nH_ratio    : %g" % nDII)

print("initial_nHM_to_nH_ratio     : %g" % nHM)

print("initial_nH2I_to_nH_ratio    : %g" % nH2I)
print("initial_nH2II_to_nH_ratio   : %g" % nH2II)


print("initial_nHDI_to_nH_ratio    : %g" % nHDI)
