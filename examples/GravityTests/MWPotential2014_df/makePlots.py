###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Darwin Roduit (yves.revaz@.epfl.ch)
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
import numpy as np
import h5py
import matplotlib.pyplot as plt


#%%Functions


def get_positions_and_time(N_snapshots, N_part, output_dir, boxsize):
    xx = np.zeros((N_part, N_snapshots))
    yy = np.zeros((N_part, N_snapshots))
    zz = np.zeros((N_part, N_snapshots))
    time = np.zeros(N_snapshots)
    pot = np.zeros((N_part, N_snapshots))
    for i in range(0, N_snapshots):
        Data = h5py.File(output_dir + "/output_%04d.hdf5" % i, "r")
        header = Data["Header"]
        time[i] = header.attrs["Time"][0]
        particles = Data["PartType1"]
        positions = particles["Coordinates"]
        xx[:, i] = positions[:, 0] - boxsize / 2.0
        yy[:, i] = positions[:, 1] - boxsize / 2.0
        zz[:, i] = positions[:, 2] - boxsize / 2.0
        pot[:, i] = particles["Potentials"][:]
    return xx, yy, zz, time, pot


def plot_orbits(x, y, z, t, color, save_fig_name_suffix):
    # Plots the orbits
    fig, ax = plt.subplots(nrows=1, ncols=4, num=1, figsize=(12, 4.1))
    fig.suptitle("Orbits", fontsize=15)
    ax[0].clear()
    ax[1].clear()

    for i in range(0, N_part):
        ax[0].plot(x[i, :], y[i, :], color[i])

    ax[0].set_aspect("equal", "box")
    ax[0].set_xlim([-300, 300])
    ax[0].set_ylim([-300, 300])
    ax[0].set_ylabel("y (kpc)")
    ax[0].set_xlabel("x (kpc)")

    for i in range(0, N_part):
        ax[1].plot(x[i, :], z[i, :], col[i], label="SWIFT solution")

    ax[1].set_aspect("equal", "box")
    ax[1].set_xlim([-100, 100])
    ax[1].set_ylim([-100, 100])
    ax[1].set_ylabel("z (kpc)")
    ax[1].set_xlabel("x (kpc)")

    for i in range(0, N_part):
        ax[2].plot(y[i, :], z[i, :], col[i])

    ax[2].set_aspect("equal", "box")
    ax[2].set_xlim([-100, 100])
    ax[2].set_ylim([-100, 100])
    ax[2].set_ylabel("z (kpc)")
    ax[2].set_xlabel("y (kpc)")
    plt.tight_layout()

    for i in range(0, N_part):
        ax[3].plot(t, np.sqrt(x[i, :] ** 2 + y[i, :] ** 2 + z[i, :] ** 2), col[i])

    ax[3].set_aspect("auto", "box")
    ax[3].set_ylim([0, 100])
    ax[3].set_ylabel("r (kpc)")
    ax[3].set_xlabel("t (kpc)")
    plt.tight_layout()

    # add the reference orbit
    data = np.genfromtxt("orbit.csv", delimiter=",", skip_header=1)
    t = data[:, 0]
    x = data[:, 1]
    y = data[:, 2]
    z = data[:, 3]

    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    ax[0].plot(x, y, "grey", alpha=0.5, lw=5)
    ax[1].plot(x, z, "grey", alpha=0.5, lw=5, label="pNbody solution")
    ax[2].plot(y, z, "grey", alpha=0.5, lw=5)
    ax[3].plot(t, r, "grey", alpha=0.5, lw=5)

    ax[1].legend()

    plt.savefig("orbits" + save_fig_name_suffix + ".png", bbox_inches="tight")
    plt.close()


#%%Plots the orbits, the deviation from the circular orbit and the deviation from the original precomputed data
# Notice that in this examples, the ouputs are set in suitable units in the parameters files.

# General parameters
N_snapshots = 1001
N_part = 1
boxsize = 1000.0  # kpc
col = ["b", "r", "c", "y", "k"]

# First type of units (kpc)
output_dir = "output_1"
save_fig_name_suffix = "_simulation_kpc"
x_1, y_1, z_1, time_1, pot_1 = get_positions_and_time(
    N_snapshots, N_part, output_dir, boxsize
)
plot_orbits(x_1, y_1, z_1, time_1, col, save_fig_name_suffix)


# Second type of units (Mpc) (no need for units conversion to kpc, this is already done by swift in the snapshots)
output_dir = "output_2"
save_fig_name_suffix = "_simulation_Mpc"
x_2, y_2, z_2, time_2, pot_2 = get_positions_and_time(
    N_snapshots, N_part, output_dir, boxsize
)
plot_orbits(x_2, y_2, z_2, time_2, col, save_fig_name_suffix)
