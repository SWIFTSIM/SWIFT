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


### Compute the length/radius of the lobe
lobe_lengths = []
lobe_radii = []

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
    gas_coordinates = data.gas.coordinates.value
    gas_ids = data.gas.particle_ids.value

    # Recenter coordinates
    boxsize = data.metadata.boxsize.value  # box size
    gas_coordinates[:, 0] -= boxsize[0] / 2
    gas_coordinates[:, 1] -= boxsize[1] / 2
    gas_coordinates[:, 2] -= boxsize[2] / 2

    # Find length and radius of lobe as the mean position of the 10 particles launched into the jets that are farthest (z-direction and r-direction) from the origin.
    jet_selection = gas_ids < 2e7
    lobe_length = np.mean(
        np.sort(np.absolute(gas_coordinates[:, 2][jet_selection]))[::-1][0:10]
    )
    lobe_radius = np.mean(
        np.sort(
            np.sqrt(
                gas_coordinates[:, 0][jet_selection] ** 2
                + gas_coordinates[:, 1][jet_selection] ** 2
            )
        )[::-1][0:10]
    )
    lobe_lengths.append(lobe_length)
    lobe_radii.append(lobe_radius)

# Make the plot
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.style.use("classic")
fig = plt.figure(figsize=(8, 6))

times = np.array(snapshots[1:-1]) * 0.005
times_in_Myr = times * 1000

# Jet/lobe length
plt.plot(
    times_in_Myr,
    lobe_lengths[1:-1],
    color="black",
    linestyle="-",
    linewidth=3,
    label="Jet/lobe length",
)

# Jet/lobe length
plt.plot(
    times_in_Myr,
    lobe_radii[1:-1],
    color="black",
    linestyle="--",
    linewidth=3,
    label="Lobe radius",
)

# Plot analytical solution
opening_angle_in_degrees = 10
theta = opening_angle_in_degrees / 180 * np.pi  # in radians
lobe_aspect_ratio = 1.084 / theta * (1 + 0.5 * theta ** 2 * 1.084 ** 2)
dens_0 = 1.47e-5  # density of the ambient medium
P_jet = 1.56e8  # jet power

# Get the lobe length and radius
lobe_length_analytical = D_jet(theta, 5 / 3, 0, dens_0, P_jet, times)
lobe_radius_analytical = lobe_length_analytical / lobe_aspect_ratio

plt.plot(
    times_in_Myr,
    lobe_length_analytical,
    color="blue",
    linestyle="-",
    linewidth=1,
    label="Kaiser & Best (2007), analytical, $L\propto t^{0.6}$",
)
plt.plot(
    times_in_Myr, lobe_radius_analytical, color="blue", linestyle="--", linewidth=1
)


plt.yscale("log")
plt.xscale("log")

# Cosmetics
plt.xlabel(r"$t$ $[\mathrm{Myr}]$", fontsize=22)
plt.ylabel(r"$L$ $[\mathrm{kpc}]$", fontsize=22)
plt.legend(loc="lower right", prop={"size": 14.5}, ncol=1, columnspacing=0.1)
plt.axis([4.5, 110, 1e0, 3e2])
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
plt.savefig("lobe_dimensions.png", bbox_inches="tight", pad_inches=0.1)
