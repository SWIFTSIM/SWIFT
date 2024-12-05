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

# Coefficient c1 in eqn 32 (Kaiser and Alexander 1997)
#   - opening_angle_deg: half-opening angle of the jet
#   - gamma: adiabatic index (here the same for jet, lobe and ambient medium)
#   - beta: slope of the density profile of the ambient medium
def c1(theta, gama, beta):
    A = (
        1.084 / theta * (1 + 0.5 * theta ** 2 * 1.084 ** 2)
    )  # approximate aspect ratio, assuming theta is small enough
    return (
        A ** 4
        / (18 * np.pi)
        * (5 - beta) ** 3
        * (gama ** 2 - 1)
        / (9 * (gama + (gama - 1) * A ** 2 * 0.5) - 4 - beta)
    ) ** (1 / (5 - beta))


# Length of the lobe (eqn 31 in Kaiser and Alexander 1997)
#   - opening_angle_deg: half-opening angle of the jet
#   - gamma: adiabatic index (here the same for jet, lobe and ambient medium)
#   - beta: slope of the density profile of the ambient medium
#   - rho_0: central density of the ambient medium
#   - P_j: jet power
#   - t: time
def D_jet(theta, gama, beta, rho_0, P_j, t):
    return c1(theta, gama, beta) * (P_j / rho_0 * t ** 3) ** (1 / 5)


# Make plot
snapshots = np.array([5, 10, 15, 20])
ages = [25, 50, 75, 100]

r_size = 250  # kpc, size of image along z-direction (y axis on image)
scale = 50  # kpc, bar scale to put on plot
slice_factor = 0.025  # thickness of slice to display
mintemp = 6.5  # minimum of color bar scale for temperature
maxtemp = 9.5  # maximum of color bar scale for temperature
aspect_ratio = 0.5  # aspect ratio of images, x:y

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(23, 10.5))
gs = gridspec.GridSpec(1, 4, hspace=0, wspace=0)

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

    # Box size
    boxsize = data.metadata.boxsize.value

    # Properties we need
    gas_positions = data.gas.coordinates.value
    gas_temperatures = data.gas.temperatures.value
    gas_masses = (
        data.gas.masses.value * 1e10
    )  # quickview complains later if this value is too small
    gas_smoothing_lengths = data.gas.smoothing_lengths.value
    gas_ids = data.gas.particle_ids.value

    # Recenter the positions
    gas_positions[:, 0] -= boxsize[0] / 2.0
    gas_positions[:, 1] -= boxsize[1] / 2.0
    gas_positions[:, 2] -= boxsize[2] / 2.0

    # Apply selection to exclude particles in the reservoir
    selection = gas_ids > 1e7

    # Also exclude particles outside the slice
    selection = (
        selection
        & (gas_positions[:, 1] < slice_factor * r_size)
        & (gas_positions[:, 1] > -slice_factor * r_size)
    )

    # Also exclude particles outside visualisation box
    selection = (
        selection
        & (gas_positions[:, 0] < aspect_ratio * r_size)
        & (gas_positions[:, 0] > -aspect_ratio * r_size)
    )
    selection = (
        selection & (gas_positions[:, 2] < r_size) & (gas_positions[:, 2] > -r_size)
    )

    gas_masses = gas_masses[selection]
    gas_temperatures = gas_temperatures[selection]
    gas_smoothing_lengths = gas_smoothing_lengths[selection]
    gas_positions = gas_positions[selection]

    # Surface density
    qv_gas = QuickView(
        gas_positions,
        hsml=gas_smoothing_lengths,
        mass=gas_masses,
        plot=False,
        logscale=True,
        r="infinity",
        p=0,
        t=90,
        extent=[-r_size, r_size, -r_size, r_size],
        x=0,
        y=0,
        z=0,
        xsize=1500,
        ysize=1500,
    )
    img_gas = qv_gas.get_image()
    ext_gas = qv_gas.get_extent()

    # Surface density times temperature
    qv_gas2 = QuickView(
        gas_positions,
        hsml=gas_smoothing_lengths,
        mass=gas_masses * gas_temperatures,
        plot=False,
        logscale=True,
        r="infinity",
        p=0,
        t=90,
        extent=[-r_size, r_size, -r_size, r_size],
        x=0,
        y=0,
        z=0,
        xsize=1500,
        ysize=1500,
    )
    img_gas2 = qv_gas2.get_image()

    # Projection-averaged temperature (logged)
    img_gas_temp = img_gas2 - img_gas

    plt.subplot(gs[i])
    # plot the temperature
    im = plt.imshow(
        img_gas_temp, vmin=mintemp, vmax=maxtemp, cmap="inferno", extent=ext_gas
    )

    # Some jet parameters for an analytical estimate
    opening_angle_in_degrees = 10
    theta = opening_angle_in_degrees / 180 * np.pi  # in radians
    lobe_aspect_ratio = 1.084 / theta * (1 + 0.5 * theta ** 2 * 1.084 ** 2)

    dens_0 = 1.47e-5  # density of the ambient medium
    P_jet = 1.56e8  # jet power
    t_jet = 0.005 * snapshots[i]  # jet age at current time

    # Get the lobe length and radius
    lobe_length = D_jet(theta, 5 / 3, 0, dens_0, P_jet, t_jet)
    lobe_radius = lobe_length / lobe_aspect_ratio

    # Plot analytical solution as ellipse
    t = np.linspace(0, 2 * np.pi, 100)
    plt.plot(
        lobe_radius * np.cos(t),
        lobe_length / 2 + lobe_length / 2 * np.sin(t),
        color="green",
        linewidth=2,
        linestyle="--",
    )

    # Add a scale
    plt.plot(
        np.linspace(-0.475 * r_size, -0.4 * r_size + scale, 2),
        [24 / 25 * r_size for x in np.linspace(-11 / 12.5 * r_size, -120, 2)],
        linewidth=4,
        color="white",
    )
    plt.text(
        0.45 * (-1.2 * r_size + scale),
        11.1 / 12.5 * r_size,
        str(scale) + " kpc",
        color="white",
        fontsize=28,
    )
    plt.text(
        1 / 12.5 * r_size,
        11.3 / 12.5 * r_size,
        "t=" + str(ages[i]) + " Myr",
        color="white",
        fontsize=28,
    )

    # Plot bounds
    plt.xlim(-aspect_ratio * r_size, aspect_ratio * r_size)
    plt.ylim(-r_size, r_size)

    # Remove ticks
    plt.xticks([])
    plt.yticks([])

# Add color bar
cax = fig.add_axes([0.93, 0.05, 0.025, 0.9])
cb = fig.colorbar(im, cax=cax, orientation="vertical")
cb.update_ticks()
cb.ax.set_ylabel(r"$\log_{10}T$ $[\mathrm{K}]$", rotation=270, fontsize=28)
cb.ax.get_yaxis().labelpad = 30
cb.ax.tick_params(labelsize=26)

# Save figure
plt.savefig("temperature_slices.png", bbox_inches="tight", pad_inches=0.1)
