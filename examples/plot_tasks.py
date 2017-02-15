#!/usr/bin/env python
"""
Usage:
    plot_tasks.py input.dat output.png [time-range-ms]

where input.dat is a thread info file for a step.  Use the '-y interval'
flag of the swift MPI commands to create these. The output plot will be
called 'output.png'. Use the time-range-ms in millisecs to produce
plots with the same time span.

This file is part of SWIFT.
Copyright (c) 2015 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
                   Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
                   Matthieu Schaller (matthieu.schaller@durham.ac.uk)
          (c) 2016 Peter W. Draper (p.w.draper@durham.ac.uk)

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
import matplotlib.collections as collections
matplotlib.use('Agg')
import pylab as pl
import numpy as np
import sys

#  Basic plot configuration.
PLOT_PARAMS = {"axes.labelsize": 10,
               "axes.titlesize": 10,
               "font.size": 12,
               "legend.fontsize": 12,
               "xtick.labelsize": 10,
               "ytick.labelsize": 10,
               "figure.figsize" : (16., 4.),
               "figure.subplot.left" : 0.03,
               "figure.subplot.right" : 0.995,
               "figure.subplot.bottom" : 0.1,
               "figure.subplot.top" : 0.99,
               "figure.subplot.wspace" : 0.,
               "figure.subplot.hspace" : 0.,
               "lines.markersize" : 6,
               "lines.linewidth" : 3.
               }
pl.rcParams.update(PLOT_PARAMS)

#  Tasks and subtypes. Indexed as in tasks.h.
TASKTYPES = ["none", "sort", "self", "pair", "sub_self", "sub_pair",
             "init", "ghost", "extra_ghost", "drift", "kick1", "kick2",
             "timestep", "send", "recv", "grav_gather_m", "grav_fft",
             "grav_mm", "grav_up", "cooling", "sourceterms", "count"]
SUBTYPES = ["none", "density", "gradient", "force", "grav", "external_grav",
            "tend", "xv", "rho", "gpart", "count"]

#  Task/subtypes of interest.
FULLTYPES = ["self/force", "self/density", "sub_self/force",
             "sub_self/density", "pair/force", "pair/density", "sub_pair/force",
             "sub_pair/density", "recv/xv", "send/xv", "recv/rho", "send/rho",
             "recv/tend", "send/tend"]

#  Get a number of colours for the various types.
colours = ["black", "gray", "rosybrown", "firebrick", "red", "darksalmon",
           "sienna", "sandybrown", "bisque", "tan", "moccasin", "gold", "darkkhaki",
           "lightgoldenrodyellow", "olivedrab", "chartreuse", "darksage", "lightgreen",
           "green", "mediumseagreen", "mediumaquamarine", "mediumturquoise", "darkslategrey",
           "cyan", "cadetblue", "skyblue", "dodgerblue", "slategray", "darkblue",
           "slateblue", "blueviolet", "mediumorchid", "purple", "magenta", "hotpink",
           "pink"]
maxcolours = len(colours)

#  Set colours of task/subtype.
TASKCOLOURS = {}
ncolours = 0
for task in TASKTYPES:
    TASKCOLOURS[task] = colours[ncolours]
    ncolours = (ncolours + 1) % maxcolours

SUBCOLOURS = {}
for task in SUBTYPES:
    SUBCOLOURS[task] = colours[ncolours]
    ncolours = (ncolours + 1) % maxcolours

for task in FULLTYPES:
    SUBCOLOURS[task] = colours[ncolours]
    ncolours = (ncolours + 1) % maxcolours

#  Show docs if help is requested.
if len( sys.argv ) == 2 and ( sys.argv[1][0:2] == "-h" or sys.argv[1][0:3] == "--h" ):
    from pydoc import help
    help( "__main__" )
    sys.exit( 0 )

#  Handle command-line.
if len( sys.argv ) != 3 and len( sys.argv ) != 4:
    print "Usage: ", sys.argv[0], "input.dat output.png [time-range-ms]"
    sys.exit(1)

infile = sys.argv[1]
outpng = sys.argv[2]
delta_t = 0
if len( sys.argv ) == 4:
    delta_t = int(sys.argv[3])

#  Read input.
data = pl.loadtxt( infile )

nthread = int(max(data[:,0])) + 1
print "Number of threads:", nthread

# Recover the start and end time
full_step = data[0,:]
tic_step = int(full_step[4])
toc_step = int(full_step[5])
CPU_CLOCK = float(full_step[-1])
data = data[1:,:]

print "CPU frequency:", CPU_CLOCK

# Avoid start and end times of zero.
data = data[data[:,4] != 0]
data = data[data[:,5] != 0]

# Calculate the time range, if not given.
delta_t = delta_t * CPU_CLOCK / 1000
if delta_t == 0:
    dt = max(data[:,5]) - min(data[:,4])
    if dt > delta_t:
        delta_t = dt
    print "Data range: ", delta_t / CPU_CLOCK * 1000, "ms"

# Once more doing the real gather and plots this time.
start_t = tic_step 
data[:,4] -= start_t
data[:,5] -= start_t
end_t = (toc_step - start_t) / CPU_CLOCK * 1000

tasks = {}
tasks[-1] = []
for i in range(nthread):
    tasks[i] = []

num_lines = pl.size(data) / 10
for line in range(num_lines):
    thread = int(data[line,0])
    tasks[thread].append({})
    tasktype = TASKTYPES[int(data[line,1])]
    subtype = SUBTYPES[int(data[line,2])]
    tasks[thread][-1]["type"] = tasktype
    tasks[thread][-1]["subtype"] = subtype
    tic = int(data[line,4]) / CPU_CLOCK * 1000
    toc = int(data[line,5]) / CPU_CLOCK * 1000
    tasks[thread][-1]["tic"] = tic
    tasks[thread][-1]["toc"] = toc
    tasks[thread][-1]["t"] = (toc + tic)/ 2
    if "self" in tasktype or "pair" in tasktype:
        fulltype = tasktype + "/" + subtype
        if fulltype in SUBCOLOURS:
            tasks[thread][-1]["colour"] = SUBCOLOURS[fulltype]
        else:
            tasks[thread][-1]["colour"] = SUBCOLOURS[subtype]
    else:
        tasks[thread][-1]["colour"] = TASKCOLOURS[tasktype]
    
for thread in range(nthread):
    tasks[thread] = sorted(tasks[thread], key=lambda l: l["t"])
            
typesseen = []
fig = pl.figure()
ax = fig.add_subplot(1,1,1)
ax.set_xlim(-delta_t * 0.03 * 1000 / CPU_CLOCK, delta_t * 1.03 * 1000 / CPU_CLOCK)
ax.set_ylim(0, nthread)
tictoc = np.zeros(2)
for i in range(nthread):

    #  Collect ranges and colours into arrays.
    tictocs = np.zeros(len(tasks[i])*2)
    colours = np.empty(len(tasks[i])*2, dtype='object')
    coloursseen = []
    j = 0
    for task in tasks[i]:
        tictocs[j] = task["tic"]
        tictocs[j+1] = task["toc"]
        colours[j] = task["colour"]
        colours[j+1] = task["colour"]
        j = j + 2
        if task["colour"] not in coloursseen:
            coloursseen.append(task["colour"])

        #  Legend support, collections don't add to this.
        if task["subtype"] != "none":
            qtask = task["type"] + "/" + task["subtype"]
        else:
            qtask = task["type"]
        if qtask not in typesseen:
            pl.plot([], [], color=task["colour"], label=qtask)
            typesseen.append(qtask)

    #  Now plot each colour, faster to use a mask to select colour ranges.
    for colour in coloursseen:
        collection = collections.BrokenBarHCollection.span_where(tictocs, ymin=i+0.05, ymax=i+0.95,
                                                                 where=colours == colour,
                                                                 facecolor=colour,
                                                                 linewidths=0)
        ax.add_collection(collection)


#  Legend and room for it.
nrow = len(typesseen) / 5
if len(typesseen) * 5 < nrow:
    nrow = nrow + 1
ax.fill_between([0, 0], nthread+0.5, nthread + nrow + 0.5, facecolor="white")
ax.set_ylim(0, nthread + nrow + 1)
ax.legend(loc=1, shadow=True, mode="expand", ncol=5)

# Start and end of time-step
ax.plot([0, 0], [0, nthread + nrow + 1], 'k--', linewidth=1)
ax.plot([end_t, end_t], [0, nthread + nrow + 1], 'k--', linewidth=1)

ax.set_xlabel("Wall clock time [ms]")
ax.set_ylabel("Thread ID" )
ax.set_yticks(pl.array(range(nthread)), True)

pl.show()
pl.savefig(outpng)
print "Graphics done, output written to", outpng

sys.exit(0)
