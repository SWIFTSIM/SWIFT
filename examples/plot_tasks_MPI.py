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
#matplotlib.use('Agg')
import pylab as pl
import numpy as np

CPU_CLOCK = 2.7e9

params = {'axes.labelsize': 10,
          'axes.titlesize': 10,
          'font.size': 12,
          'legend.fontsize': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': False,
          'figure.figsize' : (10., 4.),
          'figure.subplot.left' : 0.03,
          'figure.subplot.right' : 0.995,
          'figure.subplot.bottom' : 0.1,
          'figure.subplot.top' : 0.99,
          'figure.subplot.wspace' : 0.,
          'figure.subplot.hspace' : 0.,
          'lines.markersize' : 6,
          'lines.linewidth' : 3.,
          'text.latex.unicode': True
}
pl.rcParams.update(params)
#pl.rc('font',**{'family':'sans-serif','sans-serif':['Times']})

types = {"0": "task_none",
         "1": "task_sort",
         "2": "task_self",
         "3": "task_pair",
         "4": "task_pair",
         "5": "task_sub",
         "6": "task_ghost",
         "7": "task_kick1",
         "8": "task_kick2",
         "9": "task_send",
         "10": "task_recv",
         "11": "task_grav_pp",
         "12": "task_grav_mm",
         "13": "task_grav_up",
         "14": "task_grav_down",
         "15": "task_psort",
         "16": "task_split_cell",
         "17": "task_rewait",
         "18": "task_count"}

subtypes = {"0": "subtask_none",
            "1": "subtask_density",
            "2": "subtask_force",
            "3": "subtask_grav",
            "4": "subtask_count"}

# Assign colours for all types.
colours = ["red","blue","green","yellow","cyan","magenta","black"]
colours = colours + list(matplotlib.colors.cnames)
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


data = pl.loadtxt("thread_info_MPI.dat")

nranks = int(max(data[:,0])) + 1
print "Number of ranks:", nranks

nthread = int(max(data[:,1])) + 1
print "Number of threads:", nthread

tasks = {}
tasks[-1] = []
for i in range(nthread*nranks):
    tasks[i] = []

# Start and end times, avoiding embarrassing zeros.
tics = data[:,5]
tics = tics[tics != 0]
start_t = tics.min()
end_t = max(data[:,6])
data[:,5] -= start_t
data[:,6] -= start_t

num_lines = pl.size(data) / 10
for line in range(num_lines):
    tic = int(data[line,5]) / CPU_CLOCK * 1000
    toc = int(data[line,6]) / CPU_CLOCK * 1000
    if tic > 0 and toc > 0:
        rank = int(data[line,0])
        thread = int(data[line,1])
        index = thread*nranks + rank
        tasks[index].append({})
        tasks[index][-1]["type"] = types[ str(int(data[line,2])) ]
        tasks[index][-1]["subtype"] = subtypes[str(int(data[line,3]))]
        tasks[index][-1]["tic"] = tic
        tasks[index][-1]["toc"] = toc
        tasks[index][-1]["t"] = (toc - tic)/ 2

combtasks = {}
combtasks[-1] = []
for i in range(nthread*nranks):
    combtasks[i] = []

for index in range(nthread*nranks):
    tasks[index] = sorted(tasks[index], key=lambda l: l["t"])
    lasttype = ""
    types = []
    for task in tasks[index]:
        if task["type"] not in types:
            types.append(task["type"])
        if lasttype == "" or not lasttype == task["type"]:
            combtasks[index].append({})
            combtasks[index][-1]["type"] = task["type"]
            combtasks[index][-1]["subtype"] = task["subtype"]
            combtasks[index][-1]["tic"] = task["tic"]
            combtasks[index][-1]["toc"] = task["toc"]
            if task["type"] == 'task_self' or task["type"] == 'task_pair' or task["type"] == 'task_sub':
                combtasks[index][-1]["colour"] = subtypecolours[task["subtype"]]
            else:
                combtasks[index][-1]["colour"] = taskcolours[task["type"]]
            lasttype = task["type"]
        else:
            combtasks[index][-1]["toc"] = task["toc"]

typesseen = []
for i in range(nthread*nranks):
    for task in combtasks[i]:
        pl.fill_between([task["tic"], task["toc"]], i+0.05, i+0.95,
                        facecolor=task["colour"])
        qtask = task["type"] + "/" + task["subtype"]
        if qtask not in typesseen:
            pl.plot([], [], color=task["colour"], label=qtask)
            typesseen.append(qtask)

#  Legend and room for it.
pl.fill_between([0, 0], nthread*nranks+0.05, nthread*nranks+0.95, facecolor="white")
pl.legend(loc=1, shadow=True, mode="expand", ncol=3)

pl.xlabel("Wall clock time [ms]")
pl.xlim(0, (end_t - start_t)*1.03 * 1000 / CPU_CLOCK)
pl.ylabel("Thread ID*MPI Rank")
pl.yticks(pl.array(range(nthread*nranks)) + 0.5, pl.array(range(nthread*nranks)))

pl.show()
pl.savefig("task_graph.png")
