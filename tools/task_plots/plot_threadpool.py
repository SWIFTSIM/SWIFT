#!/usr/bin/env python3
"""
Usage:
    plot_threadpool.py [options] input.dat output.png

where input.dat is a threadpool info file for a step.  Use the '-Y interval'
flag of the swift command to create these. The output plot will be called
'output.png'. The --limit option can be used to produce plots with the same
time span and the --expand option to expand each thread line into '*expand'
lines, so that adjacent tasks of the same type can be distinguished. Other
options can be seen using the --help flag.

This file is part of SWIFT.
Copyright (c) 2015 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
                   Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
                   Matthieu Schaller (schaller@strw.leidenuniv.nl)
          (c) 2017 Peter W. Draper (p.w.draper@durham.ac.uk)

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

matplotlib.use("Agg")
import matplotlib.collections as collections
import matplotlib.ticker as plticker
import pylab as pl
import sys
import argparse

#  Handle the command line.
parser = argparse.ArgumentParser(description="Plot threadpool function graphs")

parser.add_argument("input", help="Threadpool data file (-Y output)")
parser.add_argument("outpng", help="Name for output graphic file (PNG)")
parser.add_argument(
    "-l",
    "--limit",
    dest="limit",
    help="Upper time limit in millisecs (def: depends on data)",
    default=0,
    type=float,
)
parser.add_argument(
    "-e",
    "--expand",
    dest="expand",
    help="Thread expansion factor (def: 1)",
    default=1,
    type=int,
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
    "--nolegend",
    dest="nolegend",
    help="Whether to show the legend (def: False)",
    default=False,
    action="store_true",
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
    "-m",
    "--mintic",
    dest="mintic",
    help="Value of the smallest tic (def: least in input file)",
    default=-1,
    type=int,
)

args = parser.parse_args()
infile = args.input
outpng = args.outpng
delta_t = args.limit
expand = args.expand
mintic = args.mintic

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
if mintic > 0:
    tic_step = mintic

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

#  Set colours of task/subtype.
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
    for task in sorted(SUBCOLOURS.keys()):
        print("# " + task + ": " + SUBCOLOURS[task])

tasks = {}
tasks[-1] = []
for i in range(nthread * expand):
    tasks[i] = []

#  Counters for each thread when expanding.
ecounter = []
for i in range(nthread):
    ecounter.append(0)

for i in range(len(threads)):
    thread = threads[i]

    # Expand to cover extra lines if expanding.
    ethread = thread * expand + (ecounter[thread] % expand)
    ecounter[thread] = ecounter[thread] + 1
    thread = ethread

    tasks[thread].append({})
    tasks[thread][-1]["type"] = funcs[i]
    tic = tics[i] / CPU_CLOCK
    toc = tocs[i] / CPU_CLOCK
    tasks[thread][-1]["tic"] = tic
    tasks[thread][-1]["toc"] = toc
    tasks[thread][-1]["colour"] = TASKCOLOURS[funcs[i]]

# Use expanded threads from now on.
nthread = nthread * expand

typesseen = []
fig = pl.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlim(-delta_t * 0.01 / CPU_CLOCK, delta_t * 1.01 / CPU_CLOCK)
ax.set_ylim(0, nthread)

# Fake thread is used to colour the whole range, do that first.
tictocs = []
colours = []
j = 0
for task in tasks[nthread - expand]:
    tictocs.append((task["tic"], task["toc"] - task["tic"]))
    colours.append(task["colour"])
ax.broken_barh(tictocs, [0, (nthread - 1)], facecolors=colours, linewidth=0, alpha=0.15)

# And we don't plot the fake thread.
nthread = nthread - expand
for i in range(nthread):

    #  Collect ranges and colours into arrays.
    tictocs = []
    colours = []
    j = 0
    for task in tasks[i]:
        tictocs.append((task["tic"], task["toc"] - task["tic"]))
        colours.append(task["colour"])

        #  Legend support, collections don't add to this.
        qtask = task["type"]
        if qtask not in typesseen:
            pl.plot([], [], color=task["colour"], label=qtask)
            typesseen.append(qtask)

    #  Now plot.
    ax.broken_barh(tictocs, [i + 0.05, 0.90], facecolors=colours, linewidth=0)

#  Legend and room for it.
nrow = len(typesseen) / 5
if not args.nolegend:
    ax.fill_between([0, 0], nthread + 0.5, nthread + nrow + 0.5, facecolor="white")
    ax.set_ylim(0, nthread + 0.5)
    ax.legend(
        loc=1, shadow=True, bbox_to_anchor=(0.0, 1.05, 1.0, 0.2), mode="expand", ncol=5
    )
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height * 0.8])

# Start and end of time-step
real_start_t = (mintic_step - tic_step) / CPU_CLOCK
ax.plot([real_start_t, real_start_t], [0, nthread + nrow + 1], "k--", linewidth=1)

ax.plot([end_t, end_t], [0, nthread + nrow + 1], "k--", linewidth=1)

ax.set_xlabel("Wall clock time [ms]", labelpad=0.0)
if expand == 1:
    ax.set_ylabel("Thread ID", labelpad=0)
else:
    ax.set_ylabel("Thread ID * " + str(expand), labelpad=0)
ax.set_yticks(pl.array(list(range(nthread))), True)

loc = plticker.MultipleLocator(base=expand)
ax.yaxis.set_major_locator(loc)
ax.grid(True, which="major", axis="y", linestyle="-")

pl.show()
pl.savefig(outpng)
print("Graphics done, output written to", outpng)

sys.exit(0)
