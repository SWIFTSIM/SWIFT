################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
################################################################################

# Computes the temperature evolution of the gas in a cosmological box

# Physical constants needed for internal energy to temperature conversion
k_in_J_K = 1.38064852e-23
mH_in_kg = 1.6737236e-27

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys

# plt.style.use("../../../../tools/stylesheets/mnras.mplstyle")

snap = int(sys.argv[1])

# Read the simulation data
sim = h5py.File("snapshots/snap_%04d.hdf5" % snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
z = sim["/Cosmology"].attrs["Redshift"][0]
a = sim["/Cosmology"].attrs["Scale-factor"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"][0]
kernel = sim["/HydroScheme"].attrs["Kernel function"][0]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"][0]
eta = sim["/HydroScheme"].attrs["Kernel eta"][0]
alpha = sim["/HydroScheme"].attrs["Alpha viscosity"][0]
H_mass_fraction = sim["/HydroScheme"].attrs["Hydrogen mass fraction"][0]
H_transition_temp = sim["/HydroScheme"].attrs[
    "Hydrogen ionization transition temperature"
][0]
T_initial = sim["/HydroScheme"].attrs["Initial temperature"][0]
T_minimal = sim["/HydroScheme"].attrs["Minimal temperature"][0]
git = sim["Code"].attrs["Git Revision"]

# Cosmological parameters
H_0 = sim["/Cosmology"].attrs["H0 [internal units]"][0]
gas_gamma = sim["/HydroScheme"].attrs["Adiabatic index"][0]

unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]

unit_length_in_si = 0.01 * unit_length_in_cgs
unit_mass_in_si = 0.001 * unit_mass_in_cgs
unit_time_in_si = unit_time_in_cgs

rho = sim["/PartType0/Densities"][:]
B = sim["/PartType0/MagneticFluxDensities"][:]

# Life is better in log-space
log_rho = np.log10(rho)

normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
log_B = np.log10(normB)

# Make a 2D histogram
log_rho_min = -2
log_rho_max = 3
log_B_min = -14
log_B_max = -6

bins_x = np.linspace(log_rho_min, log_rho_max, 54)
bins_y = np.linspace(log_B_min, log_B_max, 54)
# H, _, _ = np.histogram2d(log_rho, log_B, bins=[bins_x, bins_y], normed=True)
H, _, _ = np.histogram2d(log_rho, log_B, bins=[bins_x, bins_y])


# Plot the interesting quantities
plt.figure()

# plt.pcolormesh(bins_x, bins_y, np.log10(H).T)
plt.pcolormesh(bins_x, bins_y, H.T)

plt.text(-5, 8.0, "$z=%.2f$" % z)

plt.xticks(
    [-5, -4, -3, -2, -1, 0, 1, 2, 3],
    ["", "$10^{-4}$", "", "$10^{-2}$", "", "$10^0$", "", "$10^2$", ""],
)
# plt.yticks(
#    [2, 3, 4, 5, 6, 7, 8], ["$10^{2}$", "", "$10^{4}$", "", "$10^{6}$", "", "$10^8$"]
# )
plt.plot([-2, 3], [-13, 2 / 3 * (3 + 2) - 13], lw=1.5, alpha=0.7, color="red")
plt.xlabel("${\\rm Physical~Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
plt.ylabel("${\\rm Magnetic~Field}~B[{\\rm G}]$", labelpad=0)
plt.xlim(log_rho_min, log_rho_max)
plt.ylim(log_B_min, log_B_max)

# plt.tight_layout()

plt.savefig("rhoB_%04d.png" % snap, dpi=200)
