#!/usr/bin/env python3
"""
Interactive plot of a threadpool dump.

Usage:
    iplot_threadpool.py [options] input.dat

where input.dat is a threadpool info file for a step.  Use the '-Y interval'
flag of the swift or swift_mpi commands to create these (these will need to be
built with the --enable-threadpool-debugging configure option).

The plot can be scrolled and zoomed using the standard matplotlib
controls, the type of task at a point can be queried by a mouse click
(unless the --motion option is in effect when a continuous readout is
shown) the task type and tic/toc range are reported in the terminal.

Requires the tkinter module.

This file is part of SWIFT.

Copyright (C) 2022 Peter W. Draper (p.w.draper@durham.ac.uk)
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

matplotlib.use("TkAgg")
import numpy as np
import matplotlib.backends.backend_tkagg as tkagg
from matplotlib.figure import Figure
import tkinter as tk
import math
import matplotlib.collections as collections
import matplotlib.ticker as plticker
import pylab as pl
import sys
import argparse

#  Handle the command line.
parser = argparse.ArgumentParser(description="Plot threadpool function graphs")

parser.add_argument("input", help="Threadpool data file (-Y output)")
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

args = parser.parse_args()
infile = args.input
delta_t = args.limit

#  Basic plot configuration.
PLOT_PARAMS = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.figsize": (args.width, args.height),
    "figure.subplot.left": 0.03,
    "figure.subplot.right": 0.995,
    "figure.subplot.bottom": 0.09,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.0,
    "figure.subplot.hspace": 0.0,
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
}
pl.rcParams.update(PLOT_PARAMS)

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

#  Read header. First two lines.
with open(infile) as infid:
    head = [next(infid) for x in range(2)]
header = head[1][2:].strip()
header = eval(header)
nthread = int(header["num_threads"]) + 1
CPU_CLOCK = float(header["cpufreq"]) / 1000.0
print("Number of threads: ", nthread)
if args.verbose:
    print("CPU frequency:", CPU_CLOCK * 1000.0)

#  Read input.
data = pl.genfromtxt(infile, dtype=None, delimiter=" ", encoding=None)

#  Mixed types, so need to separate.
tics = []
tocs = []
funcs = []
threads = []
chunks = []
for i in data:
    if i[0] != "#":
        funcs.append(i[0].replace("_mapper", ""))
        if i[1] < 0:
            threads.append(nthread - 1)
        else:
            threads.append(i[1])
        chunks.append(i[2])
        tics.append(i[3])
        tocs.append(i[4])
tics = pl.array(tics)
tocs = pl.array(tocs)
funcs = pl.array(funcs)
threads = pl.array(threads)
chunks = pl.array(chunks)

#  Recover the start and end time
mintic_step = min(tics)
tic_step = mintic_step
toc_step = max(tocs)
print("# Min tic = ", mintic_step)

#  Calculate the time range, if not given.
delta_t = delta_t * CPU_CLOCK
if delta_t == 0:
    dt = toc_step - tic_step
    if dt > delta_t:
        delta_t = dt
    print("Data range: ", delta_t / CPU_CLOCK, "ms")

#  Once more doing the real gather and plots this time.
start_t = float(tic_step)
tics -= tic_step
tocs -= tic_step
end_t = (toc_step - start_t) / CPU_CLOCK

#  Get all "task" names and assign colours.
TASKTYPES = pl.unique(funcs)
print(TASKTYPES)

#  Set colours of tasks.
TASKCOLOURS = {}
ncolours = 0
for task in TASKTYPES:
    TASKCOLOURS[task] = colours[ncolours]
    ncolours = (ncolours + 1) % maxcolours

#  For fiddling with colours...
if args.verbose:
    print("#Selected colours:")
    for task in sorted(TASKCOLOURS.keys()):
        print("# " + task + ": " + TASKCOLOURS[task])

tasks = {}
tasks[-1] = []
for i in range(nthread):
    tasks[i] = []

for i in range(len(threads)):
    thread = threads[i]
    tasks[thread].append({})
    tasks[thread][-1]["type"] = funcs[i]
    tic = tics[i] / CPU_CLOCK
    toc = tocs[i] / CPU_CLOCK
    tasks[thread][-1]["tic"] = tic
    tasks[thread][-1]["toc"] = toc
    tasks[thread][-1]["colour"] = TASKCOLOURS[funcs[i]]

# Do the plotting.
fig = Figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlim(-delta_t * 0.01 / CPU_CLOCK, delta_t * 1.01 / CPU_CLOCK)
ax.set_ylim(0.5, nthread + 1.0)

ltics = []
ltocs = []
llabels = []
for i in range(nthread):

    #  Collect ranges and colours into arrays. Also indexed lists for lookup tables.
    tictocs = []
    colours = []
    tics = []
    tocs = []
    labels = []
    for task in tasks[i]:
        tictocs.append((task["tic"], task["toc"] - task["tic"]))
        colours.append(task["colour"])

        tics.append(task["tic"])
        tocs.append(task["toc"])
        labels.append(task["type"])

    #  Add to look up tables.
    ltics.append(tics)
    ltocs.append(tocs)
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
    def __init__(self, window, figure, motion, nthread, ltics, ltocs, llabels):
        self.window = window
        self.figure = figure
        self.motion = motion
        self.nthread = nthread
        self.ltics = ltics
        self.ltocs = ltocs
        self.llabels = llabels

    def plot(self):
        canvas = tkagg.FigureCanvasTkAgg(self.figure, master=self.window)
        wcanvas = canvas.get_tk_widget()
        wcanvas.config(width=1000, height=300)
        wcanvas.pack(side=tk.TOP, expand=True, fill=tk.BOTH)

        toolbar = tkagg.NavigationToolbar2Tk(canvas, self.window)
        toolbar.update()
        self.output = tk.StringVar()
        label = tk.Label(
            self.window, textvariable=self.output, bg="white", fg="red", bd=2
        )
        label.pack(side=tk.RIGHT, expand=True, fill=tk.X)
        wcanvas.pack(side=tk.TOP, expand=True, fill=tk.BOTH)

        canvas.draw()

        # Print task type using mouse clicks or motion.
        if self.motion:
            fig.canvas.mpl_connect("motion_notify_event", self.onclick)
        else:
            fig.canvas.mpl_connect("button_press_event", self.onclick)

        # Space bar to dump all tasks. Use with caution...
        fig.canvas.mpl_connect("key_press_event", self.dump)

    def dump(self, event):
        #  Dump all tasks to the console sorted by tic.
        xlow = float(event.inaxes.viewLim.x0)
        xhigh = float(event.inaxes.viewLim.x1)

        if event.key == " ":
            dumps = {}
            for thread in range(nthread):
                tics = self.ltics[thread]
                tocs = self.ltocs[thread]
                labels = self.llabels[thread]
                for i in range(len(tics)):
                    if (tics[i] > xlow and tics[i] < xhigh) or (
                        tocs[i] > xlow and tocs[i] < xhigh
                    ):
                        tic = "{0:.3f}".format(tics[i])
                        toc = "{0:.3f}".format(tocs[i])
                        dumps[tics[i]] = (
                            labels[i] + ",  tic/toc =  " + tic + " / " + toc
                        )
            print("")
            print("Tasks in time range: " + str(xlow) + " -> " + str(xhigh))
            for key in sorted(dumps):
                print(dumps[key])
            print("")

    def onclick(self, event):
        # Find thread, then scan for bounded task.
        try:
            thread = int(round(event.ydata)) - 1
            outstr = "none"
            if thread >= 0 and thread < self.nthread:
                tics = self.ltics[thread]
                tocs = self.ltocs[thread]
                labels = self.llabels[thread]
                for i in range(len(tics)):
                    if event.xdata > tics[i] and event.xdata < tocs[i]:
                        tic = "{0:.3f}".format(tics[i])
                        toc = "{0:.3f}".format(tocs[i])
                        outstr = labels[i] + ",  tic/toc =  " + tic + " / " + toc
                        break
            self.output.set(outstr)
        except TypeError:
            #  Out of bounds clears field.
            self.output.set("")
            pass

    def quit(self):
        self.window.destroy()


window = tk.Tk()
window.protocol("WM_DELETE_WINDOW", window.quit)
container = Container(window, fig, args.motion, nthread, ltics, ltocs, llabels)
container.plot()
window.mainloop()

sys.exit(0)
