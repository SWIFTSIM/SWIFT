###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2015 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 #                    Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
matplotlib.use('Agg')
import pylab as pl
import numpy as np
import sys

CPU_CLOCK = 2.7e9

params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'text.fontsize': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
'figure.figsize' : (12., 4.),
'figure.subplot.left'    : 0.03,
'figure.subplot.right'   : 0.995  ,
'figure.subplot.bottom'  : 0.1  ,
'figure.subplot.top'     : 0.99  ,
'figure.subplot.wspace'  : 0.  ,
'figure.subplot.hspace'  : 0.  ,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
pl.rcParams.update(params)
pl.rc('font',**{'family':'sans-serif','sans-serif':['Times']})



types = {"0": "task_none",
         "1": "task_sort",
         "2": "task_self",
         "3": "task_pair",
         "4": "task_sub",
         "5": "task_ghost",
         "6": "task_kick1",
         "7": "task_kick2",
         "8": "task_send",
         "9": "task_recv",
         "10": "task_grav_pp",
         "11": "task_grav_mm",
         "12": "task_grav_up",
         "13": "task_grav_down",
         "14": "task_psort",
         "15": "task_split_cell",
         "16": "task_rewait",
         "17": "task_count"}

subtypes = {"0": "subtask_none",
            "1": "subtask_density",
            "2": "subtask_force",
            "3": "subtask_grav",
            "4": "subtask_count"}

# Assign colours for all types.
colors = ["red","blue","green","yellow","cyan","magenta","black"]
colors = colors + list(matplotlib.colors.cnames)
index = 0

subtypecolors = {}
for key in subtypes:
    subtypecolors[subtypes[key]] = colors[index]
    print subtypes[key], " = ", colors[index]
    index = index + 1

taskcolors = {}
for key in types:
    taskcolors[types[key]] = colors[index]
    print types[key], " = ", colors[index]
    index = index + 1

data = pl.loadtxt("thread_info.dat")

nthread = int(max(data[:,0])) + 1
print "Number of threads:", nthread

tasks = {}
tasks[-1] = []
for i in range(nthread):
    tasks[i] = []


start_t = min(data[:,4])
end_t = max(data[:,5])
data[:,4] -= start_t
data[:,5] -= start_t
num_lines = pl.size(data) / 8
for line in range(num_lines):
    thread = int(data[line,0])
    tasks[thread].append({})
    tasks[thread][-1]["type"] = types[ str(int(data[line,1])) ]
    tasks[thread][-1]["subtype"] = subtypes[str(int(data[line,2]))]
    tasks[thread][-1]["tic"] = int(data[line,4]) / CPU_CLOCK * 1000
    tasks[thread][-1]["toc"] = int(data[line,5]) / CPU_CLOCK * 1000
    tasks[thread][-1]["t"] = (tasks[thread][-1]["toc"]+tasks[thread][-1]["tic"]) / 2

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
            print "seen: ", task["type"], "/", task["subtype"]
        if lasttype == "" or not lasttype == task["type"]:
            combtasks[thread].append({})
            combtasks[thread][-1]["type"] = task["type"]
            combtasks[thread][-1]["subtype"] = task["subtype"]
            combtasks[thread][-1]["tic"] = task["tic"]
            combtasks[thread][-1]["toc"] = task["toc"]
            if task["type"] == 'task_self' or task["type"] == 'task_pair' or task["type"] == 'task_sub':
                combtasks[thread][-1]["color"] = subtypecolors[task["subtype"]]
            else:
                combtasks[thread][-1]["color"] = taskcolors[task["type"]]
            lasttype = task["type"]
        else:
            combtasks[thread][-1]["toc"] = task["toc"]

for i in range(nthread):
    for task in combtasks[i]:
        pl.fill_between([task["tic"], task["toc"]], i+0.05, i+0.95, facecolor=task["color"])

pl.xlabel("Wall clock time [ms]")
pl.xlim(0, (end_t - start_t)*1.03 * 1000 / CPU_CLOCK)
pl.ylabel("Thread ID")
pl.yticks(pl.array(range(nthread)) + 0.5, pl.array(range(nthread)))
pl.savefig("task_graph.png")


