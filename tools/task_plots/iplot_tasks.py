#!/usr/bin/env python
"""
Interactive plot of a task dump.

Usage:
    iplot_tasks.py [options] input.dat

where input.dat is a thread info file for a step.  Use the '-y interval' flag
of the swift or swift_mpi commands to create these (these will need to be
built with the --enable-task-debugging configure option).

The task plot can be scrolled and zoomed using the standard matplotlib
controls, the type of task at a point can be queried by a mouse click
(unless the --motion option is in effect when a continuous readout is
shown) the task type and tic/toc range are reported in the terminal.

Requires the tkinter module.

This file is part of SWIFT.

Copyright (C) 2019 Peter W. Draper (p.w.draper@durham.ac.uk)
All Rights Reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.backends.backend_tkagg as tkagg
from matplotlib.figure import Figure
import tkinter as tk
import matplotlib.collections as collections
import matplotlib.ticker as plticker
import pylab as pl
import sys
import argparse

#  Handle the command line.
parser = argparse.ArgumentParser(description="Plot task graphs")

parser.add_argument("input", help="Thread data file (-y output)")
parser.add_argument(
    "-m",
    "--motion",
    dest="motion",
    help="Track mouse motion, otherwise clicks (def: clicks)",
    default=False,
    action="store_true",
)
parser.add_argument(
    "-l",
    "--limit",
    dest="limit",
    help="Upper time limit in millisecs (def: depends on data)",
    default=0,
    type=float,
)
parser.add_argument(
    "--height",
    dest="height",
    help="Height of plot in inches (def: 4)",
    default=4.0,
    type=float,
)
parser.add_argument(
    "--width",
    dest="width",
    help="Width of plot in inches (def: 16)",
    default=16.0,
    type=float,
)
parser.add_argument(
    "-v",
    "--verbose",
    dest="verbose",
    help="Show colour assignments and other details (def: False)",
    default=False,
    action="store_true",
)
parser.add_argument(
    "-r",
    "--rank",
    dest="rank",
    help="The rank to plot, if MPI in effect",
    default=0,
    type=int,
)

args = parser.parse_args()
infile = args.input
delta_t = args.limit
rank = args.rank

#  Basic plot configuration.
PLOT_PARAMS = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.figsize": (args.width, args.height),
    "figure.subplot.left": 0.03,
    "figure.subplot.right": 0.995,
    "figure.subplot.bottom": 0.1,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.0,
    "figure.subplot.hspace": 0.0,
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
}
pl.rcParams.update(PLOT_PARAMS)

#  Tasks and subtypes. Indexed as in tasks.h.
TASKTYPES = [
    "none",
    "sort",
    "self",
    "pair",
    "sub_self",
    "sub_pair",
    "init_grav",
    "init_grav_out",
    "ghost_in",
    "ghost",
    "ghost_out",
    "extra_ghost",
    "drift_part",
    "drift_spart",
    "drift_sink",
    "drift_bpart",
    "drift_gpart",
    "drift_gpart_out",
    "hydro_end_force",
    "kick1",
    "kick2",
    "timestep",
    "timestep_limiter",
    "timestep_sync",
    "send",
    "recv",
    "grav_long_range",
    "grav_mm",
    "grav_down_in",
    "grav_down",
    "grav_end_force",
    "cooling",
    "cooling_in",
    "cooling_out",
    "star_formation",
    "star_formation_in",
    "star_formation_out",
    "logger",
    "stars_in",
    "stars_out",
    "stars_ghost_in",
    "stars_ghost",
    "stars_ghost_out",
    "stars_sort",
    "stars_resort",
    "bh_in",
    "bh_out",
    "bh_ghost",
    "bh_swallow_ghost1",
    "bh_swallow_ghost2",
    "bh_swallow_ghost3",
    "fof_self",
    "fof_pair",
    "rt_in",
    "rt_out",
    "sink_formation",
    "rt_ghost1",
    "count",
]

SUBTYPES = [
    "none",
    "density",
    "gradient",
    "force",
    "limiter",
    "grav",
    "external_grav",
    "tend_part",
    "tend_gpart",
    "tend_spart",
    "tend_sink",
    "tend_bpart",
    "xv",
    "rho",
    "part_swallow",
    "bpart_merger",
    "gpart",
    "multipole",
    "spart",
    "stars_density",
    "stars_feedback",
    "sf_counts",
    "bpart_rho",
    "bpart_swallow",
    "bpart_feedback",
    "bh_density",
    "bh_swallow",
    "do_swallow",
    "bh_feedback",
    "sink",
    "rt_inject",
    "sink_compute_formation",
    "count",
]

#  Task/subtypes of interest.
FULLTYPES = [
    "self/limiter",
    "self/force",
    "self/gradient",
    "self/density",
    "self/grav",
    "sub_self/limiter",
    "sub_self/force",
    "sub_self/gradient",
    "sub_self/density",
    "pair/limiter",
    "pair/force",
    "pair/gradient",
    "pair/density",
    "pair/grav",
    "sub_pair/limiter",
    "sub_pair/force",
    "sub_pair/gradient",
    "sub_pair/density",
    "recv/xv",
    "send/xv",
    "recv/rho",
    "send/rho",
    "recv/tend_part",
    "send/tend_part",
    "recv/tend_gpart",
    "send/tend_gpart",
    "recv/tend_spart",
    "send/tend_spart",
    "recv/tend_bpart",
    "send/tend_bpart",
    "recv/gpart",
    "send/gpart",
    "recv/spart",
    "send/spart",
    "send/sf_counts",
    "recv/sf_counts",
    "recv/bpart",
    "send/bpart",
    "self/stars_density",
    "pair/stars_density",
    "sub_self/stars_density",
    "sub_pair/stars_density",
    "self/stars_feedback",
    "pair/stars_feedback",
    "sub_self/stars_feedback",
    "sub_pair/stars_feedback",
    "self/bh_density",
    "pair/bh_density",
    "sub_self/bh_density",
    "sub_pair/bh_density",
    "self/bh_swallow",
    "pair/bh_swallow",
    "sub_self/bh_swallow",
    "sub_pair/bh_swallow",
    "self/do_swallow",
    "pair/do_swallow",
    "sub_self/do_swallow",
    "sub_pair/do_swallow",
    "self/bh_feedback",
    "pair/bh_feedback",
    "sub_self/bh_feedback",
    "sub_pair/bh_feedback",
]

#  A number of colours for the various types. Recycled when there are
#  more task types than colours...
colours = [
    "cyan",
    "lightgray",
    "darkblue",
    "yellow",
    "tan",
    "dodgerblue",
    "sienna",
    "aquamarine",
    "bisque",
    "blue",
    "green",
    "lightgreen",
    "brown",
    "purple",
    "moccasin",
    "olivedrab",
    "chartreuse",
    "olive",
    "darkgreen",
    "green",
    "mediumseagreen",
    "mediumaquamarine",
    "darkslategrey",
    "mediumturquoise",
    "black",
    "cadetblue",
    "skyblue",
    "red",
    "slategray",
    "gold",
    "slateblue",
    "blueviolet",
    "mediumorchid",
    "firebrick",
    "magenta",
    "hotpink",
    "pink",
    "orange",
    "lightgreen",
]
maxcolours = len(colours)

#  Set colours of task/subtype.
TASKCOLOURS = {}
ncolours = 0
for task in TASKTYPES:
    TASKCOLOURS[task] = colours[ncolours]
    ncolours = (ncolours + 1) % maxcolours

SUBCOLOURS = {}
for task in FULLTYPES:
    SUBCOLOURS[task] = colours[ncolours]
    ncolours = (ncolours + 1) % maxcolours

for task in SUBTYPES:
    SUBCOLOURS[task] = colours[ncolours]
    ncolours = (ncolours + 1) % maxcolours

#  For fiddling with colours...
if args.verbose:
    print("#Selected colours:")
    for task in sorted(TASKCOLOURS.keys()):
        print(("# " + task + ": " + TASKCOLOURS[task]))
    for task in sorted(SUBCOLOURS.keys()):
        print(("# " + task + ": " + SUBCOLOURS[task]))

#  Read input.
data = pl.loadtxt(infile)

#  Do we have an MPI file?
full_step = data[0, :]
if full_step.size == 15:
    print("# MPI mode")
    mpimode = True
    ranks = list(range(int(max(data[:, 0])) + 1))
    print(("# Number of ranks:", len(ranks)))
    rankcol = 0
    threadscol = 1
    taskcol = 2
    subtaskcol = 3
    ticcol = 5
    toccol = 6
    qticcol = 13
    unlockscol = 14
else:
    print("# non MPI mode")
    mpimode = False
    rankcol = -1
    threadscol = 0
    taskcol = 1
    subtaskcol = 2
    ticcol = 4
    toccol = 5
    qticcol = 11
    unlockscol = 12

#  Get CPU_CLOCK to convert ticks into milliseconds.
CPU_CLOCK = float(full_step[-1]) / 1000.0
if args.verbose:
    print(("# CPU frequency:", CPU_CLOCK * 1000.0))

nthread = int(max(data[:, threadscol])) + 1
print(("# Number of threads:", nthread))

# Avoid start and end times of zero.
sdata = data[data[:, ticcol] != 0]
sdata = sdata[sdata[:, toccol] != 0]

# Calculate the data range, if not given.
delta_t = delta_t * CPU_CLOCK
if delta_t == 0:
    if mpimode:
        data = sdata[sdata[:, rankcol] == rank]
        full_step = data[0, :]

    tic_step = int(full_step[ticcol])
    toc_step = int(full_step[toccol])
    dt = toc_step - tic_step
    if dt > delta_t:
        delta_t = dt
    print(("# Data range: ", delta_t / CPU_CLOCK, "ms"))

# Once more doing the real gather and plots this time.
if mpimode:
    data = sdata[sdata[:, rankcol] == rank]
    full_step = data[0, :]
tic_step = int(full_step[ticcol])
toc_step = int(full_step[toccol])
print(("# Min tic = ", tic_step))
data = data[1:, :]

# Exit if no data.
if data.size == 0:
    print(("# Rank ", rank, " has no tasks"))
    sys.exit(1)

start_t = float(tic_step)
data[:, ticcol] -= start_t
data[:, toccol] -= start_t
data[:, qticcol] -= start_t
end_t = (toc_step - start_t) / CPU_CLOCK

tasks = {}
tasks[-1] = []
for i in range(nthread):
    tasks[i] = []

num_lines = pl.shape(data)[0]
for line in range(num_lines):
    thread = int(data[line, threadscol])
    tasktype = TASKTYPES[int(data[line, taskcol])]

    #  Only look at MPI?
    #if not "send" in tasktype and not "recv" in tasktype:
    #    continue

    tasks[thread].append({})
    subtype = SUBTYPES[int(data[line, subtaskcol])]
    tasks[thread][-1]["type"] = tasktype
    tasks[thread][-1]["subtype"] = subtype
    tic = float(data[line, ticcol]) / CPU_CLOCK
    toc = float(data[line, toccol]) / CPU_CLOCK
    qtic = tic - (float(data[line, qticcol]) / CPU_CLOCK)
    unlocks = int(data[line, unlockscol])
    tasks[thread][-1]["tic"] = tic
    tasks[thread][-1]["toc"] = toc
    tasks[thread][-1]["qtic"] = qtic
    tasks[thread][-1]["unlocks"] = unlocks
    if "fof" in tasktype:
        tasks[thread][-1]["colour"] = TASKCOLOURS[tasktype]
    elif ("self" in tasktype) or ("pair" in tasktype) or ("recv" in tasktype) or ("send" in tasktype):
        fulltype = tasktype + "/" + subtype
        if fulltype in SUBCOLOURS:
            tasks[thread][-1]["colour"] = SUBCOLOURS[fulltype]
        else:
            tasks[thread][-1]["colour"] = SUBCOLOURS[subtype]
    else:
        tasks[thread][-1]["colour"] = TASKCOLOURS[tasktype]

# Do the plotting.
fig = Figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlim(-delta_t * 0.01 / CPU_CLOCK, delta_t * 1.01 / CPU_CLOCK)
ax.set_ylim(0.5, nthread + 1.0)

ltics = []
lqtics = []
ltocs = []
llabels = []
lunlocks = []
for i in range(nthread):

    #  Collect ranges and colours into arrays. Also indexed lists for lookup tables.
    tictocs = []
    colours = []
    tics = []
    tocs = []
    qtics = []
    unlocks = []
    labels = []
    for task in tasks[i]:
        tictocs.append((task["tic"], task["toc"] - task["tic"]))
        colours.append(task["colour"])

        tics.append(task["tic"])
        tocs.append(task["toc"])
        qtics.append(task["qtic"])
        unlocks.append(task["unlocks"])
        labels.append(task["type"] + "/" + task["subtype"])

    #  Add to look up tables.
    ltics.append(tics)
    ltocs.append(tocs)
    lqtics.append(qtics)
    lunlocks.append(unlocks)
    llabels.append(labels)

    #  Now plot.
    ax.broken_barh(tictocs, [i + 0.55, 0.9], facecolors=colours, linewidth=0)

# Start and end of time-step
ax.plot([0, 0], [0, nthread + 1], "k--", linewidth=1)
ax.plot([end_t, end_t], [0, nthread + 1], "k--", linewidth=1)

# Labels.
ax.set_xlabel("Wall clock time [ms]")
ax.set_ylabel("Thread ID")

loc = plticker.MultipleLocator(base=1)
ax.yaxis.set_major_locator(loc)
ax.grid(True, which="major", axis="y", linestyle="-")

class Container:
    def __init__(self,  window, figure, motion, nthread, ltics, ltocs,
                 llabels, lqtics, lunlocks):
        self.window = window
        self.figure = figure
        self.motion = motion
        self.nthread = nthread
        self.ltics = ltics
        self.ltocs = ltocs
        self.llabels = llabels
        self.lqtics = lqtics
        self.lunlocks = lunlocks

    def plot(self):
        canvas = tkagg.FigureCanvasTkAgg(self.figure, master=self.window)
        wcanvas = canvas.get_tk_widget()
        wcanvas.config(width=1000, height=300)
        wcanvas.pack(side=tk.TOP, expand=True, fill=tk.BOTH)

        toolbar = tkagg.NavigationToolbar2Tk(canvas, self.window)
        toolbar.update()
        self.output = tk.StringVar()
        label = tk.Label(self.window, textvariable=self.output, bg="white", fg="red", bd=2)
        label.pack(side=tk.RIGHT, expand=True, fill=tk.X)
        wcanvas.pack(side=tk.TOP, expand=True, fill=tk.BOTH)

        canvas.draw()

        # Print task type using mouse clicks or motion.
        if self.motion:
            fig.canvas.mpl_connect('motion_notify_event', self.onclick)
        else:
            fig.canvas.mpl_connect('button_press_event', self.onclick)

    def onclick(self, event):
        # Find thread, then scan for bounded task.
        try:
            thread = int(round(event.ydata)) - 1
            if thread >= 0 and thread < self.nthread:
                tics = self.ltics[thread]
                tocs = self.ltocs[thread]
                labels = self.llabels[thread]
                qtics = self.lqtics[thread]
                unlocks = self.lunlocks[thread]
                for i in range(len(tics)):
                    if event.xdata > tics[i] and event.xdata < tocs[i]:
                        tic = "{0:.3f}".format(tics[i])
                        toc = "{0:.3f}".format(tocs[i])
                        qtic = "{0:.3f}".format(qtics[i])
                        unlock = "{0:d}".format(unlocks[i])
                        outstr = "task =  " + labels[i] + ",  tic/toc =  " + tic + " / " + toc + " qtic = " + qtic + " unlocks = " + unlock
                        self.output.set(outstr)
                        break
        except TypeError:
            #  Ignore out of bounds.
            print("out of bounds")
            pass

    def quit(self):
        self.window.destroy()

window = tk.Tk()
window.protocol("WM_DELETE_WINDOW", window.quit)
container = Container(window, fig, args.motion, nthread, ltics, ltocs,
                      llabels, lqtics, lunlocks)
container.plot()
window.mainloop()
