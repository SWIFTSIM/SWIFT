###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Darwin Roduit (darwin.roduit@epfl.ch)
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


def plot_orbits(x, y, z, color, save_fig_name_suffix):
    # Plots the orbits
    fig, ax = plt.subplots(nrows=1, ncols=3, num=1, figsize=(12, 4.1))
    fig.suptitle("Orbits", fontsize=15)
    ax[0].clear()
    ax[1].clear()

    for i in range(0, N_part):
        ax[0].plot(x[i, :], y[i, :], color[i])

    ax[0].set_aspect("equal", "box")
    ax[0].set_xlim([-35, 35])
    ax[0].set_ylim([-35, 35])
    ax[0].set_ylabel("y (kpc)")
    ax[0].set_xlabel("x (kpc)")

    for i in range(0, N_part):
        ax[1].plot(x[i, :], z[i, :], col[i])

    ax[1].set_aspect("equal", "box")
    ax[1].set_xlim([-35, 35])
    ax[1].set_ylim([-35, 35])
    ax[1].set_ylabel("z (kpc)")
    ax[1].set_xlabel("x (kpc)")
    ax[1].legend(
        [
            "Particule 1, $R = 5$ kpc",
            "Particule 2, $R = 5$ kpc",
            "Particule 3, $R = 30$ kpc",
        ]
    )

    for i in range(0, N_part):
        ax[2].plot(y[i, :], z[i, :], col[i])

    ax[2].set_aspect("equal", "box")
    ax[2].set_xlim([-35, 35])
    ax[2].set_ylim([-35, 35])
    ax[2].set_ylabel("z (kpc)")
    ax[2].set_xlabel("y (kpc)")
    plt.tight_layout()
    plt.savefig("circular_orbits" + save_fig_name_suffix + ".png", bbox_inches="tight")
    plt.close()


def plot_deviation_from_circular_orbits(x, y, z, time, color, save_fig_name_suffix):
    # Plot of the deviation from circular orbit
    R_1_th = 5.0  # kpc
    R_2_th = 5.0  # kpc
    R_3_th = 30.0  # kpc

    fig2, ax2 = plt.subplots(nrows=1, ncols=2, num=2, figsize=(12, 4.5))
    fig2.suptitle("Deviation from circular orbit", fontsize=15)
    ax2[0].clear()
    ax2[1].clear()

    # Gather the x,y and z components of each particule into one array
    pos_1 = np.array([x[0, :], y[0, :], z[0, :]])
    pos_2 = np.array([x[1, :], y[1, :], z[1, :]])
    pos_3 = np.array([x[2, :], y[2, :], z[2, :]])

    # Compute the radii
    r_1 = np.linalg.norm(pos_1, axis=0)
    error_1 = np.abs(r_1 - R_1_th) / R_1_th * 100
    r_2 = np.linalg.norm(pos_2, axis=0)
    error_2 = np.abs(r_2 - R_2_th) / R_2_th * 100
    r_3 = np.linalg.norm(pos_3, axis=0)
    error_3 = np.abs(r_3 - R_3_th) / R_3_th * 100

    ax2[0].plot(time, error_1, color[0])
    ax2[1].plot(time, error_2, color[1])
    ax2[1].plot(time, error_3, color[2])
    ax2[0].set_ylabel("Deviation (\%)")
    ax2[0].set_xlabel("Time (Gyr)")
    ax2[1].set_ylabel("Deviation (\%)")
    ax2[1].set_xlabel("Time (Gyr)")
    ax2[0].legend(["Particule 1, $R = 5$ kpc"])
    ax2[1].legend(["Particule 2, $R = 5$ kpc", "Particule 3, $R = 30$ kpc"])

    plt.tight_layout()
    plt.savefig("deviation" + save_fig_name_suffix + ".png", bbox_inches="tight")
    plt.close()
    return r_1, r_2, r_3


def plot_deviation_from_original_data(r_1, r_2, r_3, time, color, save_fig_name_suffix):
    """Make a comparison with the obtained data and ours to check nothing is broken."""
    filename = "original_radii.txt"
    r_1_original, r_2_original, r_3_original = np.loadtxt(filename)

    # Plots the deviation wrt the original data
    fig3, ax3 = plt.subplots(nrows=1, ncols=3, num=3, figsize=(12, 4.3))
    fig3.suptitle("Deviation from the original data", fontsize=15)
    ax3[0].clear()
    ax3[1].clear()
    ax3[2].clear()

    error_1 = np.abs(r_1 - r_1_original) / r_1_original * 100
    error_2 = np.abs(r_2 - r_2_original) / r_2_original * 100
    error_3 = np.abs(r_3 - r_3_original) / r_3_original * 100

    ax3[0].plot(time, error_1, col[0])
    ax3[1].plot(time, error_2, col[1])
    ax3[2].plot(time, error_3, col[2])
    ax3[0].set_ylabel("Deviation (\%)")
    ax3[0].set_xlabel("Time (Gyr)")
    ax3[1].set_ylabel("Deviation (\%)")
    ax3[1].set_xlabel("Time (Gyr)")
    ax3[2].set_ylabel("Deviation (\%)")
    ax3[2].set_xlabel("Time (Gyr)")
    ax3[0].legend(["Particule 1, $R = 5$ kpc"])
    ax3[1].legend(["Particule 2, $R = 5$ kpc"])
    ax3[2].legend(["Particule 3, $R = 30$ kpc"])
    plt.tight_layout()
    plt.savefig(
        "deviation_from_original_data" + save_fig_name_suffix + ".png",
        bbox_inches="tight",
    )
    plt.close()


#%%Plots the orbits, the deviation from the circular orbit and the deviation from the original precomputed data
# Notice that in this examples, the ouputs are set in suitable units in the parameters files.

# General parameters
N_snapshots = 201
N_part = 3
boxsize = 1000.0  # kpc
col = ["b", "r", "c", "y", "k"]

# First type of units (kpc)
output_dir = "output_1"
save_fig_name_suffix = "_simulation_kpc"
x_1, y_1, z_1, time_1, pot_1 = get_positions_and_time(
    N_snapshots, N_part, output_dir, boxsize
)
plot_orbits(x_1, y_1, z_1, col, save_fig_name_suffix)
r_11, r_21, r_31 = plot_deviation_from_circular_orbits(
    x_1, y_1, z_1, time_1, col, save_fig_name_suffix
)
plot_deviation_from_original_data(r_11, r_21, r_31, time_1, col, save_fig_name_suffix)

# Second type of units (Mpc) (no need for units conversion to kpc, this is already done by swift in the snapshots)
output_dir = "output_2"
save_fig_name_suffix = "_simulation_Mpc"
x_2, y_2, z_2, time_2, pot_2 = get_positions_and_time(
    N_snapshots, N_part, output_dir, boxsize
)
plot_orbits(x_2, y_2, z_2, col, save_fig_name_suffix)
r_12, r_22, r_32 = plot_deviation_from_circular_orbits(
    x_2, y_2, z_2, time_2, col, save_fig_name_suffix
)
# plot_deviation_from_original_data(r_12, r_22, r_32, time_2, col, save_fig_name_suffix) #does not make sense since the original data are in kpc, not in Mpc

#%%Saves our data to be the reference ones (precomputed)
# Uncomment only if corrections of the precomputed data must occur !
# Original data :  If some corrections occur in the potential default parameters, allows to correct
# the data
# filename = "original_radii.txt"
# np.savetxt(filename, (r_1, r_2, r_3))
