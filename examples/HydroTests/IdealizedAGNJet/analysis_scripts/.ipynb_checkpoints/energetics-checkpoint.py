import pylab
from pylab import *
from scipy import stats
import h5py as h5
import matplotlib
from matplotlib.colors import LogNorm
import swiftsimio as sw
import glob
import re
import os

matplotlib.use("Agg")

# Plot parameters
params = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 9,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.figsize": (3.15, 3.15),
    "figure.subplot.left": 0.01,
    "figure.subplot.right": 0.92,
    "figure.subplot.bottom": 0.01,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0,
    "figure.subplot.hspace": 0,
    "lines.markersize": 6.0,
    "lines.linewidth": 2.0,
    "axes.facecolor": "black",
    "axes.edgecolor": "white",
}
rcParams.update(params)

path_to_file = "../"
snapshot_base = "output_"

jet_and_lobe_kinetic_energies = []
jet_and_lobe_thermal_energies = []
ambient_medium_kinetic_energies = []
ambient_medium_thermal_energies = []

Msun = 1.98848e33  # g
kpc = 3.08566e21  # cm
mu = 0.6  #
m_H = 1.67262e-24  # g
k_B = 1.38065e-16

snapshots = range(22)
for i in snapshots:
    # Read in the data
    if snapshots[i] < 10:
        data = sw.load(
            path_to_file + snapshot_base + "000" + str(snapshots[i]) + ".hdf5"
        )
    if snapshots[i] >= 10 and snapshots[i] < 100:
        data = sw.load(
            path_to_file + snapshot_base + "00" + str(snapshots[i]) + ".hdf5"
        )
    if snapshots[i] >= 100 and snapshots[i] < 1000:
        data = sw.load(path_to_file + snapshot_base + "0" + str(snapshots[i]) + ".hdf5")
    if snapshots[i] >= 1000:
        data = sw.load(path_to_file + snapshot_base + "" + str(snapshots[i]) + ".hdf5")

    # Properties we need
    gas_velocities = data.gas.velocities.value
    gas_temperatures = data.gas.temperatures.value
    gas_masses = data.gas.masses.value
    gas_ids = data.gas.particle_ids.value

    # Compute the (added) internal energy using the initial value
    if snapshots[i] == 0:
        gas_temperature_0 = gas_temperatures[-1]
    gas_thermal_energies = (
        (gas_temperatures - gas_temperature_0)
        * 1.5
        * k_B
        * gas_masses
        * 10 ** 10
        * Msun
        / (mu * m_H)
    )  # in cgs

    gas_kinetic_energies = (
        0.5
        * gas_masses
        * 10 ** 10
        * Msun
        * (
            gas_velocities[:, 0] ** 2
            + gas_velocities[:, 1] ** 2
            + gas_velocities[:, 2] ** 2
        )
        * 1e10
    )  # in cgs

    # Select all kicked particles (using their ids), or all particles hotter than some threshold value (see phase space diagram)
    jet_and_lobe_selection = ((gas_ids > 0e7) & (gas_ids < 2e7)) | (
        gas_temperatures > 1e8
    )

    # Select the other particles as ambient medium, but exclude particles still in the cone
    ambient_medium_selection = np.invert(jet_and_lobe_selection)  # | (gas_ids > 1e7)

    jet_and_lobe_kinetic_energies.append(
        np.sum(gas_kinetic_energies[jet_and_lobe_selection])
    )
    jet_and_lobe_thermal_energies.append(
        np.sum(gas_thermal_energies[jet_and_lobe_selection])
    )
    ambient_medium_kinetic_energies.append(
        np.sum(gas_kinetic_energies[ambient_medium_selection])
    )
    ambient_medium_thermal_energies.append(
        np.sum(gas_thermal_energies[ambient_medium_selection])
    )

# Create plot
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.style.use("classic")
fig = plt.figure(figsize=(8, 6))

times = np.array(snapshots[1:-1]) * 0.005
times_in_Myr = times * 1000

# Plot jet and lobe kinetic and thermal energies
plt.plot(
    times_in_Myr,
    jet_and_lobe_kinetic_energies[1:-1],
    color="blue",
    linestyle="-",
    linewidth=1.5,
    label="Jet & lobe kinetic",
)
plt.plot(
    times_in_Myr,
    jet_and_lobe_thermal_energies[1:-1],
    color="red",
    linestyle="-",
    linewidth=1.5,
    label="Jet & lobe thermal",
)

# Plot total jet and lobe energy
jet_and_lobe_total_energies = np.array(jet_and_lobe_kinetic_energies) + np.array(
    jet_and_lobe_thermal_energies
)
plt.plot(
    times_in_Myr,
    jet_and_lobe_total_energies[1:-1],
    color="black",
    linestyle="-",
    linewidth=2.5,
    label="Jet & lobe total",
)

# Plot ambient medium kinetic and thermal energies
plt.plot(
    times_in_Myr,
    ambient_medium_kinetic_energies[1:-1],
    color="blue",
    linestyle="--",
    linewidth=1.5,
    label="Ambient medium kinetic",
)
plt.plot(
    times_in_Myr,
    ambient_medium_thermal_energies[1:-1],
    color="red",
    linestyle="--",
    linewidth=1.5,
    label="Ambient medium thermal",
)

# Plot total jet and lobe energy
ambient_medium_total_energies = np.array(ambient_medium_kinetic_energies) + np.array(
    ambient_medium_thermal_energies
)
plt.plot(
    times_in_Myr,
    ambient_medium_total_energies[1:-1],
    color="black",
    linestyle="--",
    linewidth=2.5,
    label="Ambient medium total",
)

# Aggregate total
total_injected_energies = jet_and_lobe_total_energies + ambient_medium_total_energies
plt.plot(
    times_in_Myr,
    total_injected_energies[1:-1],
    color="black",
    linestyle="-",
    linewidth=4,
    label="Total injected (sum over simulation)",
)

# Sanity check - how much are we injecting versus time? This should agree perfectly with the above
P_jet = 1.56e8  # in internal units
P_jet *= 6.444262e36  # jet power in cgs units
plt.plot(
    times_in_Myr,
    P_jet * times_in_Myr * 3.1555e13,
    color="grey",
    linestyle="--",
    linewidth=4,
    label="Total injected (analytically)",
)

plt.yscale("log")

# Cosmetics
plt.xlabel(r"$t$ $[\mathrm{Myr}]$", fontsize=22)
plt.ylabel(r"$\Delta E$ $[\mathrm{erg}]$", fontsize=22)
plt.legend(loc="lower right", prop={"size": 13.5}, ncol=1, columnspacing=0.1)
plt.axis([0, 105, 3e57, 5e60])
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=8,
    width=1.2,
    which="major",
    labelsize=16,
)
plt.tick_params(
    axis="x",
    direction="in",
    bottom=True,
    top=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.tick_params(
    axis="y",
    direction="in",
    left=True,
    right=True,
    length=4,
    width=0.9,
    which="minor",
    labelsize=16,
)
plt.minorticks_on()
plt.setp(plt.gca().spines.values(), linewidth=1.35)

# Save figure
plt.savefig("energetics.png", bbox_inches="tight", pad_inches=0.1)
