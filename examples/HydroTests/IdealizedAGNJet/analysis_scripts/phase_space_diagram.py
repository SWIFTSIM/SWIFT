import pylab
from pylab import *
from scipy import stats
import h5py as h5
from sphviewer.tools import QuickView
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

snapshots = np.array([5, 10, 15, 20])
ages = [25, 50, 75, 100]

# Get the initial values of temperature and density.
data = sw.load(path_to_file + snapshot_base + "000" + str(0) + ".hdf5")
gas_temperature_0 = data.gas.temperatures.value[-1]
gas_density_0 = data.gas.densities.value[-1]

# Convert density to number density
Msun = 1.98848e33  # g
kpc = 3.08566e21  # cm
mu = 0.6  #
m_H = 1.67262e-24  # g
gas_density_0 *= (10 ** 10 * Msun / kpc ** 3) / (mu * m_H)

# Make plot
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.style.use("classic")
fig = plt.figure(figsize=(12, 12))
gs = gridspec.GridSpec(2, 2, hspace=0, wspace=0)

for i in range(4):
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
    gas_temperatures = data.gas.temperatures.value
    gas_ids = data.gas.particle_ids.value

    # Load density and convert to number density
    gas_densities = data.gas.densities.value * (10 ** 10 * Msun / kpc ** 3) / (mu * m_H)

    # Make a selection for jet and lobe particles only
    selection_jet_only = (gas_ids > 1e7) & (gas_ids < 2e7)

    # Selection for everything else
    selection_ambient_medium = np.invert(selection_jet_only)

    plt.subplot(gs[i])

    # Plot ambient medium and jet/lobe separately
    plt.scatter(
        gas_densities[selection_ambient_medium],
        gas_temperatures[selection_ambient_medium],
        marker="o",
        color="black",
        s=2,
        alpha=0.2,
        label="Ambient medium",
    )
    plt.scatter(
        gas_densities[selection_jet_only],
        gas_temperatures[selection_jet_only],
        marker="o",
        color="red",
        s=2,
        alpha=0.2,
        label="Directly kicked",
    )

    # Mark the initial temperature/density
    plt.scatter(
        gas_density_0,
        gas_temperature_0,
        marker="x",
        color="orange",
        s=100,
        linewidth=2,
        label="Initial state",
    )

    # The line of constant pressure (pressure equilibrium)
    plt.plot(
        [1e-10, 1e10],
        (gas_density_0 / np.array([1e-10, 1e10])) * gas_temperature_0,
        color="blue",
        linestyle="--",
        linewidth=1,
        label="Pressure equilibrium",
    )

    # Arbitrary division line between lobe and ambient medium, everything hotter than this value is lobe
    T_crit = 1e8  # K
    plt.plot(
        [1e-10, 1e10],
        [T_crit, T_crit],
        color="blue",
        linestyle=":",
        linewidth=1,
        label="Division between lobe" + str("\n") + " and ambient medium",
    )

    handles, labels = plt.gca().get_legend_handles_labels()
    order1 = [0, 1, 2]
    order2 = [3, 4]

    plt.yscale("log")
    plt.xscale("log")

    # Cosmetics
    if i == 0 or i == 2:
        plt.ylabel(r"$T$ $[\mathrm{K}]$", fontsize=20)
    plt.xlabel(r"$n_\mathrm{H}$ $[\mathrm{cm}^{-3}]$", fontsize=18)

    if i == 0:
        plt.legend(
            [handles[idx] for idx in order1],
            [labels[idx] for idx in order1],
            loc="upper right",
            prop={"size": 12.5},
            ncol=1,
            columnspacing=0.1,
        )
    if i == 1:
        plt.legend(
            [handles[idx] for idx in order2],
            [labels[idx] for idx in order2],
            loc="upper right",
            prop={"size": 12.5},
            ncol=1,
            columnspacing=0.1,
        )

    if i == 0:
        plt.axis([1e-4, 6e-2, 1.0001e6, 5e9])
    else:
        plt.axis([1e-4, 6e-2, 1.0001e6, 5e9])
    if i == 1 or i == 3:
        plt.yticks([1e6, 1e7, 1e8, 1e9], ["", "", "", ""])

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
plt.savefig("phase_space_diagram.png", bbox_inches="tight", pad_inches=0.1)
