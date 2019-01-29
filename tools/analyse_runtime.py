#!/usr/bin/env python

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

import re
import sys
import matplotlib

matplotlib.use("Agg")
from pylab import *

# Plot parameters
params = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": True,
    "figure.figsize": (6.45, 6.45),
    "figure.subplot.left": 0.06,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.06,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.21,
    "figure.subplot.hspace": 0.13,
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
    "text.latex.unicode": True,
}
rcParams.update(params)

threshold = 0.008

num_files = len(sys.argv) - 1

labels = [
    ["Gpart assignment", 1],
    ["Mesh comunication", 1],
    ["Forward Fourier transform", 1],
    ["Green function", 1],
    ["Backwards Fourier transform", 1],
    ["engine_recompute_displacement_constraint:", 1],
    ["engine_exchange_top_multipoles:", 1],
    ["updating particle counts", 1],
    ["engine_estimate_nr_tasks:", 1],
    ["Making gravity tasks", 1],
    ["Making hydro tasks", 1],
    ["Splitting tasks", 1],
    ["Counting and linking tasks", 1],
    ["Setting super-pointers", 1],
    ["Making extra hydroloop tasks", 1],
    ["Making extra starsloop tasks", 1],
    ["Linking gravity tasks", 1],
    ["Creating send tasks", 1],
    ["Exchanging cell tags", 1],
    ["Creating recv tasks", 1],
    ["Counting number of foreign particles", 1],
    ["Recursively linking foreign arrays", 1],
    ["Setting unlocks", 1],
    ["Ranking the tasks", 1],
    ["scheduler_reweight:", 1],
    ["space_list_useful_top_level_cells:", 1],
    ["space_rebuild:", 1],
    ["engine_drift_all:", 0],
    ["engine_unskip:", 0],
    ["engine_collect_end_of_step:", 0],
    ["engine_launch:", 0],
    ["writing particle properties", 0],
    ["engine_repartition:", 0],
    ["engine_exchange_cells:", 1],
    ["Dumping restart files", 0],
    ["engine_print_stats:", 0],
    ["engine_marktasks:", 1],
    ["Reading initial conditions", 0],
    ["engine_print_task_counts:", 0],
    ["engine_drift_top_multipoles:", 0],
    ["Communicating rebuild flag", 0],
    ["engine_split:", 0],
    ["space_init", 0],
    ["engine_init", 0],
    ["engine_repartition_trigger:", 0],
    ["VR Collecting top-level cell info", 0],
    ["VR Collecting particle info", 0],
    ["VR Invokation of velociraptor", 0],
    ["VR Copying group information back", 0]
]
times = np.zeros(len(labels))
counts = np.zeros(len(labels))

cols = [
    "0.5",
    "#332288",
    "#88CCEE",
    "#44AA99",
    "#117733",
    "#999933",
    "#DDCC77",
    "#CC6677",
    "#882255",
    "#AA4499",
    "#661100",
    "#6699CC",
    "#AA4466",
    "#4477AA",
]

total_time = 0
lastline = ""

for i in range(num_files):

    filename = sys.argv[i + 1]
    print("Analysing %s" % filename)

    # Open stdout file
    file = open(filename, "r")

    # Search the different phrases
    for line in file:

        # Loop over the possbile labels
        for i in range(len(labels)):

            # Extract the different blocks
            if re.search("%s took" % labels[i][0], line):
                counts[i] += 1.0
                times[i] += float(
                    re.findall(r"[+-]?((\d+\.?\d*)|(\.\d+))", line)[-1][0]
                )

        # Find the last line with meaningful output (avoid crash report, batch system stuff....)
        if re.findall(r"\[[0-9]{4}\][ ]\[*", line) or re.findall(
            r"^\[[0-9]*[.][0-9]+\][ ]", line
        ):
            lastline = line

    # Total run time
    total_time += float(re.findall(r"[+-]?(\[[0-9]\])?(\[[0-9]*[.][0-9]*\])+", lastline)[0][1][1:-1])

# Conver to seconds
times /= 1000.0

# Total time
total_measured_time = np.sum(times)
print("\nTotal measured time: %.3f s" % total_measured_time)

print("Total time: %f  s\n" % total_time)

# Ratios
time_ratios = times / total_time

# Better looking labels
for i in range(len(labels)):
    labels[i][0] = labels[i][0].replace("_", " ")
    labels[i][0] = labels[i][0].replace(":", "")
    labels[i][0] = labels[i][0].title()

times = np.array(times)
time_ratios = np.array(time_ratios)

# Sort in order of importance
order = np.argsort(-times)
times = times[order]
counts = counts[order]
time_ratios = time_ratios[order]
labels = [labels[i] for i in order]

# Keep only the important components
important_times = [0.0]
important_ratios = [0.0]
important_is_rebuild = [0]
important_labels = ["Others (all below %.1f\%%)" % (threshold * 100)]
need_print = True
print("Time spent in the different code sections:")
for i in range(len(labels)):
    if time_ratios[i] > threshold:
        important_times.append(times[i])
        important_ratios.append(time_ratios[i])
        important_is_rebuild.append(labels[i][1])
        important_labels.append(labels[i][0])
    else:
        if need_print:
            print("Elements in 'Other' category (<%.1f%%):" % (threshold * 100))
            need_print = False
        important_times[0] += times[i]
        important_ratios[0] += time_ratios[i]

    print(" - '%-40s' (%5d calls, time: %.4fs): %.4f%%" % (labels[i][0], counts[i], times[i], time_ratios[i] * 100))

# Anything unaccounted for?
print(
    "\nUnaccounted for: %.4f%%\n"
    % (100 * (total_time - total_measured_time) / total_time)
)

important_ratios = np.array(important_ratios)
important_is_rebuild = np.array(important_is_rebuild)


figure()

def func(pct):
    return "$%4.2f\\%%$" % pct


pie, _, _ = pie(
    important_ratios,
    explode=important_is_rebuild * 0.2,
    autopct=lambda pct: func(pct),
    textprops=dict(color="0.1", fontsize=14),
    labeldistance=0.7,
    pctdistance=0.85,
    startangle=-15,
    colors=cols,
)
legend(pie, important_labels, title="SWIFT operations", loc="upper left")

savefig("time_pie.pdf", dpi=150)
