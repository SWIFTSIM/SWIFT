#!/usr/bin/env python
###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2015 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 #                    Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 #                    Peter W. Draper (p.w.draper@durham.ac.uk)
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
 ##############################################################################


import matplotlib
matplotlib.use("Agg")
import pylab as pl
import numpy as np
import sys

#  Handle command-line.
if len( sys.argv ) != 4:
    print "Usage: ", sys.argv[0], "input.dat rank output.png"
    sys.exit(1)

infile = sys.argv[1]
rank = int(sys.argv[2])
outpng = sys.argv[3]

#  CPU ticks per second.
CPU_CLOCK = 2.7e9

params = {"axes.labelsize": 10,
          "axes.titlesize": 10,
          "font.size": 12,
          "legend.fontsize": 12,
          "xtick.labelsize": 10,
          "ytick.labelsize": 10,
          "figure.figsize" : (12., 4.),
          "figure.subplot.left" : 0.03,
          "figure.subplot.right" : 0.995,
          "figure.subplot.bottom" : 0.1,
          "figure.subplot.top" : 0.99,
          "figure.subplot.wspace" : 0.,
          "figure.subplot.hspace" : 0.,
          "lines.markersize" : 6,
          "lines.linewidth" : 3.
}
pl.rcParams.update(params)

types = {"0": "none",
         "1": "sort",
         "2": "self",
         "3": "pair",
         "4": "sub",
         "5": "ghost",
         "6": "kick1",
         "7": "kick2",
         "8": "send",
         "9": "recv",
         "10": "grav_pp",
         "11": "grav_mm",
         "12": "grav_up",
         "13": "grav_down",
         "14": "psort",
         "15": "split_cell",
         "16": "rewait",
         "17": "count"}

subtypes = {"0": "none",
            "1": "density",
            "2": "force",
            "3": "grav",
            "4": "count"}


# Assign colours for all types.
#colours = ["blue","green","red","cyan","magenta","yellow","black"]
#colours = colours + list(matplotlib.colors.cnames)
colours = list(matplotlib.colors.cnames)

index = 0
subtypecolours = {}
for key in subtypes:
    subtypecolours[subtypes[key]] = colours[index]
    print subtypes[key], " = ", colours[index]
    index = index + 1

taskcolours = {}
for key in types:
    taskcolours[types[key]] = colours[index]
    print types[key], " = ", colours[index]
    index = index + 1

#  Read input.
data = pl.loadtxt( infile )

nranks = int(max(data[:,0])) + 1
print "Number of ranks:", nranks
if rank < 0 or rank >= nranks:
    print "Cannot read rank: ", rank, ", out of range"
    sys.exit(1)

nthread = int(max(data[:,1])) + 1
print "Number of threads:", nthread

tasks = {}
tasks[-1] = []
for i in range(nthread):
    tasks[i] = []

# Pick data for our rank, and avoid start and end times of zero.
sdata = data[data[:,5] != 0]
sdata = sdata[sdata[:,6] != 0]
sdata = sdata[sdata[:,0] == rank]

start_t = min(sdata[:,5])
end_t = max(sdata[:,6])
sdata[:,5] -= start_t
sdata[:,6] -= start_t

#np.savetxt("dump.txt", sdata, fmt="%d")

num_lines = pl.size(sdata) / 10
for line in range(num_lines):
    thread = int(sdata[line,1])
    tasks[thread].append({})
    tasks[thread][-1]["type"] = types[str(int(sdata[line,2]))]
    tasks[thread][-1]["subtype"] = subtypes[str(int(sdata[line,3]))]
    tic = int(sdata[line,5]) / CPU_CLOCK * 1000
    toc = int(sdata[line,6]) / CPU_CLOCK * 1000
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
                combtasks[thread][-1]["colour"] = subtypecolours[task["subtype"]]
            else:
                combtasks[thread][-1]["colour"] = taskcolours[task["type"]]
            lasttype = task["type"]
        else:
            combtasks[thread][-1]["toc"] = task["toc"]

typesseen = []
for i in range(nthread):
    for task in combtasks[i]:
        pl.fill_between([task["tic"], task["toc"]], i+0.05, i+0.95,
                        facecolor=task["colour"], linewidths=0)
        if task["subtype"] != "none":
            qtask = task["type"] + "/" + task["subtype"]
        else:
            qtask = task["type"]
        if qtask not in typesseen:
             pl.plot([], [], color=task["colour"], label=qtask)
             typesseen.append(qtask)

#  Legend and room for it.
pl.fill_between([0, 0], nthread, nthread+len(typesseen)/3, facecolor="white")
pl.legend(loc=1, shadow=True, mode="expand", ncol=3)

pl.xlabel("Wall clock time [ms]")
pl.xlim(0, (end_t - start_t)*1.03 * 1000 / CPU_CLOCK)
pl.ylabel("Thread ID for MPI Rank " + str(rank) )
pl.yticks(pl.array(range(nthread)) + 0.5, pl.array(range(nthread)))

pl.show()
pl.savefig(outpng)
