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
import pylab as pl
import sys
import argparse

#  Handle the command line.
parser = argparse.ArgumentParser(description="Analyse page dumps for statistics on particle data")

parser.add_argument("input", help="page data file (-y output)")

args = parser.parse_args()
infile = args.input

data = pl.loadtxt(infile)

print("# cells:", len(data))

node0 =0
node1 =0
for cell in data:
    if cell[1] == 1:
        node0 = node0 + 1
    if cell[2] == 1:
        node1 = node1 + 1

print("# node0", node0)
print("# node1", node1)
