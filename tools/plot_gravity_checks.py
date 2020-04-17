#!/usr/bin/env python

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

params = {
    "axes.labelsize": 14,
    "axes.titlesize": 18,
    "font.size": 11,
    "legend.fontsize": 12,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "text.usetex": True,
    "figure.figsize": (12, 9),
    "figure.subplot.left": 0.06,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.06,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.14,
    "figure.subplot.hspace": 0.14,
    "lines.markersize": 6,
    "lines.linewidth": 3.0
}
plt.rcParams.update(params)
plt.rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})

min_error = 1e-7
max_error = 3e-1
num_bins = 64

# Construct the bins
bin_edges = np.linspace(np.log10(min_error), np.log10(max_error), num_bins + 1)
bin_size = (np.log10(max_error) - np.log10(min_error)) / num_bins
bins = 0.5 * (bin_edges[1:] + bin_edges[:-1])
bin_edges = 10 ** bin_edges
bins = 10 ** bins

# Colours
cols = ["#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677"]

# Time-step to plot
step = int(sys.argv[1])
periodic = int(sys.argv[2])

# Find the files for the different expansion orders
order_list = glob.glob("gravity_checks_swift_step%.4d_order*.dat" % step)
num_order = len(order_list)

# Get the multipole orders
order = np.zeros(num_order)
for i in range(num_order):
    order[i] = int(order_list[i][35])
order = sorted(order)
order_list = sorted(order_list)

# Read the exact accelerations first
if periodic:
    data = np.loadtxt("gravity_checks_exact_periodic_step%.4d.dat" % step)
else:
    data = np.loadtxt("gravity_checks_exact_step%.4d.dat" % step)
exact_ids = data[:, 0]
exact_pos = data[:, 1:4]
exact_a = data[:, 4:7]
exact_pot = data[:, 7]
# Sort stuff
sort_index = np.argsort(exact_ids)
exact_ids = exact_ids[sort_index]
exact_pos = exact_pos[sort_index, :]
exact_a = exact_a[sort_index, :]
exact_pot = exact_pot[sort_index]
exact_a_norm = np.sqrt(exact_a[:, 0] ** 2 + exact_a[:, 1] ** 2 + exact_a[:, 2] ** 2)

print("Number of particles tested:", np.size(exact_ids))

# Start the plot
plt.figure()

ax1 = plt.subplot(231)
ax2 = plt.subplot(232)
ax3 = plt.subplot(233)
ax4 = plt.subplot(234)
ax5 = plt.subplot(235)

count = 0

# Get the Gadget-2 data if existing
if periodic:
    gadget2_file_list = glob.glob("forcetest_gadget2_periodic.txt")
else:
    gadget2_file_list = glob.glob("forcetest_gadget2.txt")
if len(gadget2_file_list) != 0:

    gadget2_data = np.loadtxt(gadget2_file_list[0])
    gadget2_ids = gadget2_data[:, 0]
    gadget2_pos = gadget2_data[:, 1:4]
    gadget2_a_exact = gadget2_data[:, 4:7]
    gadget2_a_grav = gadget2_data[:, 7:10]

    # Sort stuff
    sort_index = np.argsort(gadget2_ids)
    gadget2_ids = gadget2_ids[sort_index]
    gadget2_pos = gadget2_pos[sort_index, :]
    gadget2_a_exact = gadget2_a_exact[sort_index, :]
    gadget2_exact_a_norm = np.sqrt(
        gadget2_a_exact[:, 0] ** 2
        + gadget2_a_exact[:, 1] ** 2
        + gadget2_a_exact[:, 2] ** 2
    )
    gadget2_a_grav = gadget2_a_grav[sort_index, :]

    # Cross-checks
    if not np.array_equal(exact_ids, gadget2_ids):
        print("Comparing different IDs !")

    if np.max(np.abs(exact_pos - gadget2_pos) / np.abs(gadget2_pos)) > 1e-6:
        print("Comparing different positions ! max difference:")
        index = np.argmax(
            exact_pos[:, 0] ** 2
            + exact_pos[:, 1] ** 2
            + exact_pos[:, 2] ** 2
            - gadget2_pos[:, 0] ** 2
            - gadget2_pos[:, 1] ** 2
            - gadget2_pos[:, 2] ** 2
        )
        print(
            "Gadget2 (id=%d):" % gadget2_ids[index],
            gadget2_pos[index, :],
            "exact (id=%d):" % exact_ids[index],
            exact_pos[index, :],
            "\n",
        )

    diff = np.abs(exact_a_norm - gadget2_exact_a_norm) / np.abs(gadget2_exact_a_norm)
    max_diff = np.max(diff)
    if max_diff > 2e-6:
        print("Comparing different exact accelerations !")
        print(
            "Median=",
            np.median(diff),
            "Mean=",
            np.mean(diff),
            "99%=",
            np.percentile(diff, 99),
        )
        print("max difference ( relative diff =", max_diff, "):")
        # index = np.argmax(exact_a[:,0]**2 + exact_a[:,1]**2 + exact_a[:,2]**2 - gadget2_a_exact[:,0]**2 - gadget2_a_exact[:,1]**2 - gadget2_a_exact[:,2]**2)
        index = np.argmax(diff)
        print(
            "a_exact --- Gadget2:",
            gadget2_a_exact[index, :],
            "exact:",
            exact_a[index, :],
        )
        print(
            "pos ---     Gadget2: (id=%d):" % gadget2_ids[index],
            gadget2_pos[index, :],
            "exact (id=%d):" % gadget2_ids[index],
            gadget2_pos[index, :],
            "\n",
        )

    # Compute the error norm
    diff = gadget2_a_exact - gadget2_a_grav

    norm_diff = np.sqrt(diff[:, 0] ** 2 + diff[:, 1] ** 2 + diff[:, 2] ** 2)
    norm_a = np.sqrt(
        gadget2_a_exact[:, 0] ** 2
        + gadget2_a_exact[:, 1] ** 2
        + gadget2_a_exact[:, 2] ** 2
    )

    norm_error = norm_diff / norm_a
    error_x = abs(diff[:, 0]) / norm_a
    error_y = abs(diff[:, 1]) / norm_a
    error_z = abs(diff[:, 2]) / norm_a

    # Bin the error
    norm_error_hist, _ = np.histogram(norm_error, bins=bin_edges, density=False) / (
        np.size(norm_error) * bin_size
    )
    error_x_hist, _ = np.histogram(error_x, bins=bin_edges, density=False) / (
        np.size(norm_error) * bin_size
    )
    error_y_hist, _ = np.histogram(error_y, bins=bin_edges, density=False) / (
        np.size(norm_error) * bin_size
    )
    error_z_hist, _ = np.histogram(error_z, bins=bin_edges, density=False) / (
        np.size(norm_error) * bin_size
    )

    norm_median = np.median(norm_error)
    median_x = np.median(error_x)
    median_y = np.median(error_y)
    median_z = np.median(error_z)

    norm_per99 = np.percentile(norm_error, 99)
    per99_x = np.percentile(error_x, 99)
    per99_y = np.percentile(error_y, 99)
    per99_z = np.percentile(error_z, 99)

    norm_max = np.max(norm_error)
    max_x = np.max(error_x)
    max_y = np.max(error_y)
    max_z = np.max(error_z)

    print("Gadget-2 ---- ")
    print("Norm: median= %f 99%%= %f max= %f" % (norm_median, norm_per99, norm_max))
    print("X   : median= %f 99%%= %f max= %f" % (median_x, per99_x, max_x))
    print("Y   : median= %f 99%%= %f max= %f" % (median_y, per99_y, max_y))
    print("Z   : median= %f 99%%= %f max= %f" % (median_z, per99_z, max_z))
    print("")

    ax4.text(
        min_error * 1.5,
        1.55,
        "$50\\%%\\rightarrow%.5f~~ 99\\%%\\rightarrow%.5f$" % (norm_median, norm_per99),
        ha="left",
        va="top",
        alpha=0.8,
    )
    ax4.semilogx(bins, norm_error_hist, "k--", label="Gadget-2", alpha=0.8)

    ax1.semilogx(bins, error_x_hist, "k--", label="Gadget-2", alpha=0.8)
    ax1.text(
        min_error * 1.5,
        1.55,
        "$50\\%%\\rightarrow%.5f~~ 99\\%%\\rightarrow%.5f$" % (median_x, per99_x),
        ha="left",
        va="top",
        alpha=0.8,
    )

    ax2.semilogx(bins, error_y_hist, "k--", label="Gadget-2", alpha=0.8)
    ax2.text(
        min_error * 1.5,
        1.55,
        "$50\\%%\\rightarrow%.5f~~ 99\\%%\\rightarrow%.5f$" % (median_y, per99_y),
        ha="left",
        va="top",
        alpha=0.8,
    )

    ax3.semilogx(bins, error_z_hist, "k--", label="Gadget-2", alpha=0.8)
    ax3.text(
        min_error * 1.5,
        1.55,
        "$50\\%%\\rightarrow%.5f~~ 99\\%%\\rightarrow%.5f$" % (median_z, per99_z),
        ha="left",
        va="top",
        alpha=0.8,
    )

    count += 1


# Plot the different histograms
for i in range(num_order):
    data = np.loadtxt(order_list[i])
    ids = data[:, 0]
    pos = data[:, 1:4]
    a_grav = data[:, 4:7]
    pot = data[:, 7]

    # Sort stuff
    sort_index = np.argsort(ids)
    ids = ids[sort_index]
    pos = pos[sort_index, :]
    a_grav = a_grav[sort_index, :]
    pot = pot[sort_index]

    # Cross-checks
    if not np.array_equal(exact_ids, ids):
        print("Comparing different IDs !")

    if np.max(np.abs(exact_pos - pos) / np.abs(pos)) > 1e-6:
        print("Comparing different positions ! max difference:")
        index = np.argmax(
            exact_pos[:, 0] ** 2
            + exact_pos[:, 1] ** 2
            + exact_pos[:, 2] ** 2
            - pos[:, 0] ** 2
            - pos[:, 1] ** 2
            - pos[:, 2] ** 2
        )
        print(
            "SWIFT (id=%d):" % ids[index],
            pos[index, :],
            "exact (id=%d):" % exact_ids[index],
            exact_pos[index, :],
            "\n",
        )

    # Compute the error norm
    diff = exact_a - a_grav
    diff_pot = exact_pot - pot

    # Correct for different normalization of potential
    print("Difference in normalization of potential:", np.mean(diff_pot), end=" ")
    print(
        "std_dev=",
        np.std(diff_pot),
        "99-percentile:",
        np.percentile(diff_pot, 99) - np.median(diff_pot),
        "1-percentile:",
        np.median(diff_pot) - np.percentile(diff_pot, 1),
    )

    exact_pot -= np.mean(diff_pot)
    diff_pot = exact_pot - pot

    norm_diff = np.sqrt(diff[:, 0] ** 2 + diff[:, 1] ** 2 + diff[:, 2] ** 2)

    norm_error = norm_diff / exact_a_norm
    error_x = abs(diff[:, 0]) / exact_a_norm
    error_y = abs(diff[:, 1]) / exact_a_norm
    error_z = abs(diff[:, 2]) / exact_a_norm
    error_pot = abs(diff_pot) / abs(exact_pot)

    # Bin the error
    norm_error_hist, _ = np.histogram(norm_error, bins=bin_edges, density=False) / (
        np.size(norm_error) * bin_size
    )
    error_x_hist, _ = np.histogram(error_x, bins=bin_edges, density=False) / (
        np.size(norm_error) * bin_size
    )
    error_y_hist, _ = np.histogram(error_y, bins=bin_edges, density=False) / (
        np.size(norm_error) * bin_size
    )
    error_z_hist, _ = np.histogram(error_z, bins=bin_edges, density=False) / (
        np.size(norm_error) * bin_size
    )
    error_pot_hist, _ = np.histogram(error_pot, bins=bin_edges, density=False) / (
        np.size(norm_error) * bin_size
    )

    norm_median = np.median(norm_error)
    median_x = np.median(error_x)
    median_y = np.median(error_y)
    median_z = np.median(error_z)
    median_pot = np.median(error_pot)

    norm_per99 = np.percentile(norm_error, 99)
    per99_x = np.percentile(error_x, 99)
    per99_y = np.percentile(error_y, 99)
    per99_z = np.percentile(error_z, 99)
    per99_pot = np.percentile(error_pot, 99)

    norm_max = np.max(norm_error)
    max_x = np.max(error_x)
    max_y = np.max(error_y)
    max_z = np.max(error_z)
    max_pot = np.max(error_pot)

    print("Order %d ---- " % order[i])
    print("Norm: median= %f 99%%= %f max= %f" % (norm_median, norm_per99, norm_max))
    print("X   : median= %f 99%%= %f max= %f" % (median_x, per99_x, max_x))
    print("Y   : median= %f 99%%= %f max= %f" % (median_y, per99_y, max_y))
    print("Z   : median= %f 99%%= %f max= %f" % (median_z, per99_z, max_z))
    print("Pot : median= %f 99%%= %f max= %f" % (median_pot, per99_pot, max_pot))
    print("")

    ax1.semilogx(
        bins, error_x_hist, color=cols[i], label="SWIFT m-poles order %d" % order[i]
    )
    ax1.text(
        min_error * 1.5,
        1.5 - count / 10.0,
        "$50\\%%\\rightarrow%.5f~~ 99\\%%\\rightarrow%.5f$" % (median_x, per99_x),
        ha="left",
        va="top",
        color=cols[i],
    )
    ax2.semilogx(
        bins, error_y_hist, color=cols[i], label="SWIFT m-poles order %d" % order[i]
    )
    ax2.text(
        min_error * 1.5,
        1.5 - count / 10.0,
        "$50\\%%\\rightarrow%.5f~~ 99\\%%\\rightarrow%.5f$" % (median_y, per99_y),
        ha="left",
        va="top",
        color=cols[i],
    )
    ax3.semilogx(
        bins, error_z_hist, color=cols[i], label="SWIFT m-poles order %d" % order[i]
    )
    ax3.text(
        min_error * 1.5,
        1.5 - count / 10.0,
        "$50\\%%\\rightarrow%.5f~~ 99\\%%\\rightarrow%.5f$" % (median_z, per99_z),
        ha="left",
        va="top",
        color=cols[i],
    )
    ax4.semilogx(
        bins, norm_error_hist, color=cols[i], label="SWIFT m-poles order %d" % order[i]
    )
    ax4.text(
        min_error * 1.5,
        1.5 - count / 10.0,
        "$50\\%%\\rightarrow%.5f~~ 99\\%%\\rightarrow%.5f$" % (norm_median, norm_per99),
        ha="left",
        va="top",
        color=cols[i],
    )
    ax5.semilogx(
        bins, error_pot_hist, color=cols[i], label="SWIFT m-poles order %d" % order[i]
    )
    ax5.text(
        min_error * 1.5,
        1.5 - count / 10.0,
        "$50\\%%\\rightarrow%.5f~~ 99\\%%\\rightarrow%.5f$" % (median_pot, per99_pot),
        ha="left",
        va="top",
        color=cols[i],
    )

    count += 1


ax1.set_xlabel("$\delta a_x/|\overrightarrow{a}_{exact}|$")
# plt.ylabel("Density")
ax1.set_xlim(min_error, max_error)
ax1.set_ylim(0, 1.75)
# plt.legend(loc="center left")

ax2.set_xlabel("$\delta a_y/|\overrightarrow{a}_{exact}|$")
# plt.ylabel("Density")
ax2.set_xlim(min_error, max_error)
ax2.set_ylim(0, 1.75)
# plt.legend(loc="center left")

ax3.set_xlabel("$\delta a_z/|\overrightarrow{a}_{exact}|$")
# plt.ylabel("Density")
ax3.set_xlim(min_error, max_error)
ax3.set_ylim(0, 1.75)

ax4.set_xlabel("$|\delta \overrightarrow{a}|/|\overrightarrow{a}_{exact}|$")
# plt.ylabel("Density")
ax4.set_xlim(min_error, max_error)
ax4.set_ylim(0, 2.5)
ax4.legend(loc="upper left", fontsize=8)

ax5.set_xlabel("$\delta \phi/\phi_{exact}$")
# plt.ylabel("Density")
ax5.set_xlim(min_error, max_error)
ax5.set_ylim(0, 1.75)
# plt.legend(loc="center left")


plt.savefig("gravity_checks_step%.4d.png" % step, dpi=200)
plt.savefig("gravity_checks_step%.4d.pdf" % step, dpi=200)
