#!/usr/bin/env python

################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
from timed_functions import labels

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
    "figure.figsize": (12.45, 6.45),
    "figure.subplot.left": 0.06,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.06,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.0,
    "figure.subplot.hspace": 0.0,
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
    "hatch.linewidth": 4,
}
rcParams.update(params)

threshold = 0.008

num_files = len(sys.argv) - 1

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

tasks = [
    "dead time",
    "drift",
    "sorts",
    "resort",
    "hydro",
    "gravity",
    "feedback",
    "black holes",
    "cooling",
    "star formation",
    "limiter",
    "sync",
    "time integration",
    "mpi",
    "pack",
    "fof",
    "others",
    "sink",
    "neutrino",
    "RT",
    "CSDS",
    "total",
]

times_tasks = np.zeros(len(tasks))
counts_tasks = np.zeros(len(tasks))

# Zoom vs background breakdown (only populated if the zoom report produced by
# zoom_scheduler_report_task_times() is present in the log).
zoom_names = ["zoom", "background", "mixed (zoom + background)"]
times_zoom_tasks = np.zeros((3, len(tasks)))
counts_zoom_tasks = np.zeros((3, len(tasks)))

total_time = 0
lastline = ""

for i in range(num_files):

    # First analyse the code sections

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
    total_time += float(
        re.findall(r"[+-]?(\[[0-9]\])?(\[[0-9]*[.][0-9]*\])+", lastline)[0][1][1:-1]
    )

    file.close()

    # Now analyse the task times

    file = open(filename, "r")

    # Track which zoom report block we are in (-1 = outside the zoom report).
    zoom_block = -1

    # Search the different phrases
    for line in file:

        # Handle lines from the zoom-specific task time report (if present).
        if "zoom_scheduler_report_task_times:" in line:

            # Header lines switch the current block (check mixed first: its
            # header also mentions zoom and background).
            if "mixed" in line and "task categories" in line:
                zoom_block = 2
                continue
            if "zoom task categories" in line:
                zoom_block = 0
                continue
            if "background task categories" in line:
                zoom_block = 1
                continue
            if "Fraction of total" in line:
                zoom_block = -1
                continue

            # Category lines: accumulate into the current block.
            if zoom_block >= 0:
                for i in range(len(tasks)):
                    if re.search("%s" % tasks[i], line):
                        counts_zoom_tasks[zoom_block][i] += 1.0
                        times_zoom_tasks[zoom_block][i] += float(
                            re.findall(r":[ ]*[-+]?\d*\.\d+|\d+ ms", line)[0][1:]
                        )
            continue

        # Loop over the possbile labels
        for i in range(len(tasks)):

            # Extract the different blocks
            # Match only the standard report; zoom_scheduler_report_task_times
            # lines contain this string as a substring and must be excluded.
            if re.search(r"(?<!zoom_)scheduler_report_task_times: \*\*\*  ", line):
                if re.search("%s" % tasks[i], line):
                    counts_tasks[i] += 1.0
                    times_tasks[i] += float(
                        re.findall(r":[ ]*[-+]?\d*\.\d+|\d+ ms", line)[0][1:]
                    )

    file.close()

# Conver to seconds
times /= 1000.0
times_tasks /= 1000.0
times_zoom_tasks /= 1000.0

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
times_tasks = np.array(times_tasks)

times_tasks_ratios = times_tasks / times_tasks[-1]
times_tasks_ratios = np.array(times_tasks_ratios)

# Sort in order of importance
order = np.argsort(-times)
times = times[order]
counts = counts[order]
time_ratios = time_ratios[order]
labels = [labels[i] for i in order]

# Remove the regexp escapes to make the labels prettier
for i in range(len(labels)):
    labels[i][0] = labels[i][0].replace("\\", "")

# Keep only the important components
important_times = [0.0]
important_ratios = [0.0]
important_is_rebuild = [0]
important_is_fof = [0]
important_is_VR = [0]
important_is_mesh = [0]
important_labels = ["Others (all below %.1f\\%%)" % (threshold * 100)]
need_print = True
print("Time spent in the different code sections:")
for i in range(len(labels)):
    if time_ratios[i] > threshold:
        important_times.append(times[i])
        important_ratios.append(time_ratios[i])
        important_is_rebuild.append(labels[i][1] == 1)
        important_is_fof.append(labels[i][1] == 2)
        important_is_VR.append(labels[i][1] == 3)
        important_is_mesh.append(labels[i][1] == 4)
        important_labels.append(labels[i][0])
    else:
        if need_print:
            print("Elements in 'Other' category (<%.1f%%):" % (threshold * 100))
            need_print = False
        important_times[0] += times[i]
        important_ratios[0] += time_ratios[i]

    print(
        " - '%-40s' (%5d calls, time: %.4fs): %.4f%%"
        % (labels[i][0], counts[i], times[i], time_ratios[i] * 100)
    )

# Anything unaccounted for?
print(
    "\nUnaccounted for: %.4f%%\n"
    % (100 * (total_time - total_measured_time) / total_time)
)

important_ratios = np.array(important_ratios)
important_is_rebuild = np.array(important_is_rebuild)

print("Time spent in the different task categories (i.e. inside engine_launch()):")

for i in range(len(tasks) - 1):
    print(
        " - '%-40s' (%5d calls): %.4f%%"
        % (tasks[i], counts_tasks[i], 100.0 * times_tasks_ratios[i])
    )
print("")

# Print the zoom vs background breakdown, if the report was found in the logs.
if counts_zoom_tasks.sum() > 0:
    for b in range(3):
        total_b = np.sum(times_zoom_tasks[b])
        if total_b == 0.0:
            continue
        print(
            "Time spent in %s task categories (i.e. inside engine_launch()):"
            % zoom_names[b]
        )
        for i in range(len(tasks) - 1):
            print(
                " - '%-40s' (%5d calls): %8.4f%% within type, %8.4f%% of total"
                % (
                    tasks[i],
                    counts_zoom_tasks[b][i],
                    100.0 * times_zoom_tasks[b][i] / total_b,
                    100.0 * times_zoom_tasks[b][i] / times_tasks[-1],
                )
            )
        print("")

figure()

# Main code sections
subplot(121)


def func(pct):
    return "$%4.2f\\%%$" % pct


code_pie, _, _ = pie(
    important_ratios,
    explode=important_is_rebuild * 0.2,
    autopct=lambda pct: func(pct),
    textprops=dict(color="0.1", fontsize=14),
    labeldistance=0.7,
    pctdistance=0.85,
    startangle=-15,
    colors=cols,
)

# Use hashing for the FOF and VR wedges
for i in range(len(code_pie)):
    if important_is_fof[i]:
        code_pie[i].set_hatch("+")
        code_pie[i].set_edgecolor(code_pie[i].get_facecolor())
        code_pie[i].set_fill(False)
for i in range(len(code_pie)):
    if important_is_VR[i]:
        code_pie[i].set_hatch("*")
        code_pie[i].set_edgecolor(code_pie[i].get_facecolor())
        code_pie[i].set_fill(False)
for i in range(len(code_pie)):
    if important_is_mesh[i]:
        code_pie[i].set_hatch(".")
        code_pie[i].set_edgecolor(code_pie[i].get_facecolor())
        code_pie[i].set_fill(False)

legend(code_pie, important_labels, title="SWIFT operations", loc="upper left")

# Tasks
subplot(122)

tasks_pie, _, _ = pie(
    times_tasks_ratios[:-1],
    autopct=lambda pct: func(pct),
    textprops=dict(color="0.1", fontsize=14),
    labeldistance=0.7,
    pctdistance=0.85,
    startangle=-15,
    colors=cols,
)

legend(tasks_pie, tasks, title="SWIFT task categories", loc="upper left")

savefig("time_pie.pdf", dpi=150)

# Zoom vs background pie, if the report was found in the logs. Combining the
# cell types in one pie keeps every wedge relative to the total task time.
if counts_zoom_tasks.sum() > 0 and times_tasks[-1] > 0:
    figure()

    pie_times = []
    pie_labels = []
    pie_colors = []
    for b in range(3):
        for i in range(len(tasks) - 1):
            if times_zoom_tasks[b][i] > 0:
                pie_times.append(times_zoom_tasks[b][i])
                pie_labels.append("%s: %s" % (zoom_names[b], tasks[i]))
                pie_colors.append(cols[i % len(cols)])

    # Include time not assigned to a cell type so the pie sums to the total.
    other_time = times_tasks[-1] - np.sum(pie_times)
    if other_time > 0:
        pie_times.append(other_time)
        pie_labels.append("dead/unclassified")
        pie_colors.append("0.5")

    zoom_pie, _, _ = pie(
        pie_times,
        autopct=lambda pct: func(pct),
        textprops=dict(color="0.1", fontsize=14),
        labeldistance=0.7,
        pctdistance=0.85,
        startangle=-15,
        colors=pie_colors,
    )
    legend(
        zoom_pie,
        pie_labels,
        title="SWIFT task categories by cell type",
        loc="upper left",
    )

    savefig("time_pie_zoom.pdf", dpi=150)
