###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
#                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

import matplotlib

matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py
import numpy as np
import glob
import os.path


def find_indices(a, b):
    result = np.zeros(len(b))
    for i in range(len(b)):
        result[i] = ((np.where(a == b[i]))[0])[0]

    return result


# Plot parameters
params = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": True,
    "figure.figsize": (9.90, 6.45),
    "figure.subplot.left": 0.1,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.1,
    "figure.subplot.top": 0.95,
    "figure.subplot.wspace": 0.2,
    "figure.subplot.hspace": 0.2,
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
    "text.latex.unicode": True,
}
rcParams.update(params)
rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})


# Number of snapshots and elements
newest_snap_name = max(glob.glob("output_*.hdf5"), key=os.path.getctime)
n_snapshots = int(newest_snap_name.replace("output_", "").replace(".hdf5", "")) + 1
n_elements = 9

# Read the simulation data
sim = h5py.File("output_0000.hdf5", "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]
stellar_mass = sim["/PartType4/Masses"][0]

# Units
unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]
unit_temp_in_cgs = sim["/Units"].attrs["Unit temperature in cgs (U_T)"]
unit_vel_in_cgs = unit_length_in_cgs / unit_time_in_cgs
unit_energy_in_cgs = unit_mass_in_cgs * unit_vel_in_cgs * unit_vel_in_cgs
unit_length_in_si = 0.01 * unit_length_in_cgs
unit_mass_in_si = 0.001 * unit_mass_in_cgs
unit_time_in_si = unit_time_in_cgs
unit_density_in_cgs = unit_mass_in_cgs * unit_length_in_cgs ** -3
unit_pressure_in_cgs = unit_mass_in_cgs / unit_length_in_cgs * unit_time_in_cgs ** -2
unit_int_energy_in_cgs = unit_energy_in_cgs / unit_mass_in_cgs
unit_entropy_in_cgs = unit_energy_in_cgs / unit_temp_in_cgs
Gyr_in_cgs = 3.155e16
Msun_in_cgs = 1.989e33

box_energy = zeros(n_snapshots)
box_mass = zeros(n_snapshots)
box_star_mass = zeros(n_snapshots)
box_metal_mass = zeros(n_snapshots)
element_mass = zeros((n_snapshots, n_elements))
t = zeros(n_snapshots)

# Read data from snapshots
for i in range(n_snapshots):
    print("reading snapshot " + str(i))
    # Read the simulation data
    sim = h5py.File("output_%04d.hdf5" % i, "r")
    t[i] = sim["/Header"].attrs["Time"][0]
    # ids = sim["/PartType0/ParticleIDs"][:]

    masses = sim["/PartType0/Masses"][:]
    box_mass[i] = np.sum(masses)

    star_masses = sim["/PartType4/Masses"][:]
    box_star_mass[i] = np.sum(star_masses)

    metallicities = sim["/PartType0/Metallicities"][:]
    box_metal_mass[i] = np.sum(metallicities * masses)

    internal_energies = sim["/PartType0/InternalEnergies"][:]
    box_energy[i] = np.sum(masses * internal_energies)

# Plot the interesting quantities
figure()

# Box mass --------------------------------
subplot(221)
plot(
    t[1:] * unit_time_in_cgs / Gyr_in_cgs,
    (box_mass[1:] - box_mass[0]) * unit_mass_in_cgs / Msun_in_cgs,
    linewidth=0.5,
    color="k",
    marker="*",
    ms=0.5,
    label="swift",
)
xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
ylabel("Change in total gas particle mass (Msun)", labelpad=2)
ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

# Box metal mass --------------------------------
subplot(222)
plot(
    t[1:] * unit_time_in_cgs / Gyr_in_cgs,
    (box_metal_mass[1:] - box_metal_mass[0]) * unit_mass_in_cgs / Msun_in_cgs,
    linewidth=0.5,
    color="k",
    ms=0.5,
    label="swift",
)
xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
ylabel("Change in total metal mass of gas particles (Msun)", labelpad=2)
ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

# Box energy --------------------------------
subplot(223)
plot(
    t[1:] * unit_time_in_cgs / Gyr_in_cgs,
    (box_energy[1:] - box_energy[0]) * unit_energy_in_cgs,
    linewidth=0.5,
    color="k",
    ms=0.5,
    label="swift",
)
xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
ylabel("Change in total energy of gas particles (erg)", labelpad=2)
ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

# Box mass --------------------------------
subplot(224)
plot(
    t[1:] * unit_time_in_cgs / Gyr_in_cgs,
    (box_mass[1:] - box_mass[0]) * unit_mass_in_cgs / Msun_in_cgs,
    linewidth=0.5,
    color="k",
    marker="*",
    ms=0.5,
    label="gas",
)
plot(
    t[1:] * unit_time_in_cgs / Gyr_in_cgs,
    (box_star_mass[1:] - box_star_mass[n_snapshots - 1])
    * unit_mass_in_cgs
    / Msun_in_cgs,
    linewidth=0.5,
    color="r",
    marker="*",
    ms=0.5,
    label="stars",
)
plot(
    t[1:] * unit_time_in_cgs / Gyr_in_cgs,
    (box_star_mass[1:] - box_star_mass[n_snapshots - 1] + box_mass[1:] - box_mass[0])
    * unit_mass_in_cgs
    / Msun_in_cgs,
    linewidth=0.5,
    color="g",
    marker="*",
    ms=0.5,
    label="total",
)
xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
ylabel("Change in total gas particle mass (Msun)", labelpad=2)
ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
legend()

savefig("box_evolution.png", dpi=200)
