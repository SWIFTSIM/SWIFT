#!/usr/bin/env python
"""
Usage:
    analyse_pcounters.py [options] input.dat

where input.dat is a thread info file for a step (MPI or non-MPI). Use the
'-y interval' flag of the swift and swift_mpi commands to create these
(you will also need to configure with the --enable-task-debugging option).

The output is an analysis of the task timings, including deadtime per thread
and step, total amount of time spent for each task type, for the whole step
and per thread and the minimum and maximum times spent per task type.

This file is part of SWIFT.
Copyright (c) 2017 Peter W. Draper (p.w.draper@durham.ac.uk)
Copright (c) 2018 STFC (author email aidan.chalk@stfc.ac.uk)

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
import matplotlib.pyplot as plt
import pylab as pl
import sys
import argparse
from os import listdir
from os.path import isfile, join


#  Handle the command line.

onlyfiles = [f for f in listdir('.') if isfile(join('.', f))]
onlyoutputs = [f for f in onlyfiles if "page_info-step" in f]

results = []
results2 = []
for i in range(1,len(onlyoutputs)):
    infile = "page_info-step{}.dat".format(i)
    
    data = pl.loadtxt(infile)
    
#    print("# cells:", len(data))
    
    node0 =0
    node1 =0
    for cell in data:
        if cell[1] == 1:
            node0 = node0 + 1
        if cell[2] == 1:
            node1 = node1 + 1
    
#    print("# node0", node0)
#    print("# node1", node1)
    results.append(node1)
    results2.append(node0)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(range(len(results)),results, c='b', label='node1')
ax1.scatter(range(len(results2)),results2, c='r', label='node0')
plt.savefig("foo.png")
