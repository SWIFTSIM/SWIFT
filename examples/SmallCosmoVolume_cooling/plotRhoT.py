################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
from pylab import *
import h5py
import os.path

# Plot parameters
params = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 9,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": True,
    "figure.figsize": (3.15, 3.15),
    "figure.subplot.left": 0.15,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.13,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.12,
    "lines.markersize": 6,
    "lines.linewidth": 2.0,
    "text.latex.unicode": True,
}
rcParams.update(params)
rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})

snap = int(sys.argv[1])

# Read the simulation data
sim = h5py.File("snap_%04d.hdf5" % snap, "r")
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

# Primoridal ean molecular weight as a function of temperature
def mu(T, H_frac=H_mass_fraction, T_trans=H_transition_temp):
    if T > T_trans:
        return 4.0 / (8.0 - 5.0 * (1.0 - H_frac))
    else:
        return 4.0 / (1.0 + 3.0 * H_frac)


# Temperature of some primoridal gas with a given internal energy
def T(u, H_frac=H_mass_fraction, T_trans=H_transition_temp):
    T_over_mu = (gas_gamma - 1.0) * u * mH_in_kg / k_in_J_K
    ret = np.ones(np.size(u)) * T_trans

    # Enough energy to be ionized?
    mask_ionized = T_over_mu > (T_trans + 1) / mu(T_trans + 1, H_frac, T_trans)
    if np.sum(mask_ionized) > 0:
        ret[mask_ionized] = T_over_mu[mask_ionized] * mu(T_trans * 10, H_frac, T_trans)

    # Neutral gas?
    mask_neutral = T_over_mu < (T_trans - 1) / mu((T_trans - 1), H_frac, T_trans)
    if np.sum(mask_neutral) > 0:
        ret[mask_neutral] = T_over_mu[mask_neutral] * mu(0, H_frac, T_trans)

    return ret


rho = sim["/PartType0/Density"][:]
u = sim["/PartType0/InternalEnergy"][:]

# Compute the temperature
u *= unit_length_in_si ** 2 / unit_time_in_si ** 2
u /= a ** (3 * (gas_gamma - 1.0))
Temp = T(u)

# Compute the physical density
rho *= unit_mass_in_cgs / unit_length_in_cgs ** 3
rho /= a ** 3
rho /= mH_in_kg

# Life is better in log-space
log_T = np.log10(Temp)
log_rho = np.log10(rho)


# Make a 2D histogram
log_rho_min = -6
log_rho_max = 3
log_T_min = 1
log_T_max = 8

bins_x = np.linspace(log_rho_min, log_rho_max, 54)
bins_y = np.linspace(log_T_min, log_T_max, 54)
H, _, _ = histogram2d(log_rho, log_T, bins=[bins_x, bins_y], normed=True)


# Plot the interesting quantities
figure()

pcolormesh(bins_x, bins_y, np.log10(H).T)

text(-5, 8.0, "$z=%.2f$" % z)

xticks(
    [-5, -4, -3, -2, -1, 0, 1, 2, 3],
    ["", "$10^{-4}$", "", "$10^{-2}$", "", "$10^0$", "", "$10^2$", ""],
)
yticks(
    [2, 3, 4, 5, 6, 7, 8], ["$10^{2}$", "", "$10^{4}$", "", "$10^{6}$", "", "$10^8$"]
)
xlabel("${\\rm Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
ylabel("${\\rm Temperature}~T~[{\\rm K}]$", labelpad=2)
xlim(-5.2, 3.2)
ylim(1, 8.5)

savefig("rhoT_%04d.png" % snap, dpi=200)
