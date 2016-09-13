#!/usr/bin/env python
"""
Usage:
    plot_tasks_MPI.py input.dat png-output-prefix [time-range-ms]

where input.dat is a thread info file for a step of an MPI run.  Use the '-y
interval' flag of the swift MPI commands to create these. The output plots
will be called 'png-output-prefix<mpi-rank>.png', i.e. one each for all the
threads in each MPI rank. Use the time-range-ms in millisecs to produce
plots with the same time span.

See the command 'process_plot_tasks_MPI' to efficiently wrap this command to
process a number of thread info files and create an HTML file to view them.

This file is part of SWIFT.

Copyright (C) 2015 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
                   Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
                   Matthieu Schaller (matthieu.schaller@durham.ac.uk)
                   Peter W. Draper (p.w.draper@durham.ac.uk)
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
import matplotlib.collections as collections
matplotlib.use("Agg")
import pylab as pl
import numpy as np
import sys
#import warnings
#warnings.simplefilter("error")

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
TASKTYPES = ["none", "sort", "self", "pair", "sub_self", "sub_pair", "init",
             "ghost", "extra_ghost", "kick", "kick_fixdt", "send", "recv",
             "grav_gather_m", "grav_fft", "grav_mm", "grav_up",
             "grav_external", "cooling", "count"]

TASKCOLOURS = {"none": "black",
               "sort": "lightblue",
               "self": "greenyellow",
               "pair": "navy",
               "sub_self": "greenyellow",
               "sub_pair": "navy",
               "init": "indigo",
               "ghost": "cyan",
               "extra_ghost": "cyan",
               "kick": "green",
               "kick_fixdt": "green",
               "send": "yellow",
               "recv": "magenta",
               "grav_gather_m": "mediumorchid",
               "grav_fft": "mediumnightblue",
               "grav_mm": "mediumturquoise",
               "grav_up": "mediumvioletred",
               "grav_external": "darkred",
               "cooling": "darkblue",
               "count": "powerblue"}

SUBTYPES = ["none", "density", "gradient", "force", "grav", "tend", "count"]

SUBCOLOURS = {"none": "black",
              "density": "red",
              "gradient": "powerblue",
              "force": "blue",
              "grav": "indigo",
              "tend": "grey",
              "count": "black"}

#  Show docs if help is requested.
if len( sys.argv ) == 2 and ( sys.argv[1][0:2] == "-h" or sys.argv[1][0:3] == "--h" ):
    from pydoc import help
    help( "__main__" )
    sys.exit( 0 )

#  Handle command-line.
if len( sys.argv ) != 3 and len( sys.argv ) != 4:
    print "Usage: ", sys.argv[0], "input.dat png-output-prefix [time-range-ms]"
    sys.exit(1)


infile = sys.argv[1]
outbase = sys.argv[2]
delta_t = 0
if len( sys.argv ) == 4:
    delta_t = int(sys.argv[3])

#  Read input.
data = pl.loadtxt( infile )

# Recover the start and end time
full_step = data[0,:]
tic_step = int(full_step[5])
toc_step = int(full_step[6])
CPU_CLOCK = float(full_step[-1])

print "CPU frequency:", CPU_CLOCK


nranks = int(max(data[:,0])) + 1
print "Number of ranks:", nranks
nthread = int(max(data[:,1])) + 1
print "Number of threads:", nthread

# Avoid start and end times of zero.
sdata = data[data[:,5] != 0]
sdata = sdata[sdata[:,6] != 0]

# Each rank can have different clock (compute node), but we want to use the
# same delta times range for comparisons, so we suck it up and take the hit of
# precalculating this, unless the user knows better.
delta_t = delta_t * CPU_CLOCK / 1000
if delta_t == 0:
    for rank in range(nranks):
        data = sdata[sdata[:,0] == rank]
        dt = max(data[:,6]) - min(data[:,5])
        if dt > delta_t:
            delta_t = dt
    print "Data range: ", delta_t / CPU_CLOCK * 1000, "ms"


# Once more doing the real gather and plots this time.
for rank in range(nranks):
    data = sdata[sdata[:,0] == rank]

    full_step = data[0,:]
    tic_step = int(full_step[5])
    toc_step = int(full_step[6])
    data = data[1:,:]
    typesseen = []

    #  Dummy image for ranks that have no tasks.
    if data.size == 0:
        print "rank ", rank, " has no tasks"
        fig = pl.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_xlim(-delta_t * 0.03 * 1000 / CPU_CLOCK, delta_t * 1.03 * 1000 / CPU_CLOCK)
        ax.set_ylim(0, nthread)
        start_t = tic_step
        end_t = (toc_step - start_t) / CPU_CLOCK * 1000
    else:

        start_t = tic_step
        data[:,5] -= start_t
        data[:,6] -= start_t
        end_t = (toc_step - start_t) / CPU_CLOCK * 1000

        tasks = {}
        tasks[-1] = []
        for i in range(nthread):
            tasks[i] = []

        num_lines = pl.shape(data)[0]
        for line in range(num_lines):
            thread = int(data[line,1])
            tasks[thread].append({})
            tasks[thread][-1]["type"] = TASKTYPES[int(data[line,2])]
            tasks[thread][-1]["subtype"] = SUBTYPES[int(data[line,3])]
            tic = int(data[line,5]) / CPU_CLOCK * 1000
            toc = int(data[line,6]) / CPU_CLOCK * 1000
            tasks[thread][-1]["tic"] = tic
            tasks[thread][-1]["toc"] = toc
            tasks[thread][-1]["t"] = (toc + tic)/ 2

        combtasks = {}
        combtasks[-1] = []
        for i in range(nthread):
            combtasks[i] = []

        for thread in range(nthread):
            tasks[thread] = sorted(tasks[thread], key=lambda l: l["t"])
            lasttype = ""
            types = []
            for task in tasks[thread]:
                if task["type"] not in types:
                    types.append(task["type"])
                if lasttype == "" or not lasttype == task["type"]:
                    combtasks[thread].append({})
                    combtasks[thread][-1]["type"] = task["type"]
                    combtasks[thread][-1]["subtype"] = task["subtype"]
                    combtasks[thread][-1]["tic"] = task["tic"]
                    combtasks[thread][-1]["toc"] = task["toc"]
                    if task["type"] == "self" or task["type"] == "pair" or task["type"] == "sub":
                        combtasks[thread][-1]["colour"] = SUBCOLOURS[task["subtype"]]
                    else:
                        combtasks[thread][-1]["colour"] = TASKCOLOURS[task["type"]]
                    lasttype = task["type"]
                else:
                    combtasks[thread][-1]["toc"] = task["toc"]

        fig = pl.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_xlim(-delta_t * 0.03 * 1000 / CPU_CLOCK, delta_t * 1.03 * 1000 / CPU_CLOCK)
        ax.set_ylim(0, nthread)
        tictoc = np.zeros(2)
        for i in range(nthread):

            #  Collect ranges and colours into arrays.
            tictocs = np.zeros(len(combtasks[i])*2)
            colours = np.empty(len(combtasks[i])*2, dtype='object')
            coloursseen = []
            j = 0
            for task in combtasks[i]:
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
                collection = collections.BrokenBarHCollection.span_where(tictocs,
                                                                         ymin=i+0.05,
                                                                         ymax=i+0.95,
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
    if data.size > 0:
        ax.legend(loc=1, shadow=True, mode="expand", ncol=5)

    # Start and end of time-step
    ax.plot([0, 0], [0, nthread + nrow + 1], 'k--', linewidth=1)
    ax.plot([end_t, end_t], [0, nthread + nrow + 1], 'k--', linewidth=1)

    ax.set_xlabel("Wall clock time [ms]")
    ax.set_ylabel("Thread ID for MPI Rank " + str(rank) )
    ax.set_yticks(pl.array(range(nthread)), True)

    pl.show()
    outpng = outbase + str(rank) + ".png"
    pl.savefig(outpng)
    print "Graphics done, output written to", outpng

sys.exit(0)
