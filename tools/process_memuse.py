#!/usr/bin/env python
"""
Usage:
    process_memuse.py output.dat

Parse the output of a run of SWIFT to convert the memuse output dumps into a
timeseries of memory use. Also outputs use in memory per labelled type.

This file is part of SWIFT.
Copyright (c) 2019 Peter W. Draper (p.w.draper@durham.ac.uk)

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

import sys
from collections import OrderedDict

#  Command-line arguments.
if len(sys.argv) != 2:
    print "usage: ", sys.argv[0], " memuse_report<>.dat"
    sys.exit(1)

memuse = OrderedDict()
labels = {}
totalmem = 0
process_use = ""
peak = 0.0

with open(sys.argv[1]) as infile:
    print "# tic label MB"
    for line in infile:
        if line[0] == "#":
            if "# Current use:" in line:
                process_use = line[14:-1]
        else:
            tic , adr, rank, step, allocated, label, size = line.split()
            rank = int(rank)
            step = int(step)
            allocated = int(allocated)
            size = int(size)

            doprint = True
            if allocated == 1:
                #  Allocation.
                totalmem = totalmem + size
                if not adr in memuse:
                    memuse[adr] = [size]
                    labels[adr] = label
                else:
                    memuse[adr].append(size)
            else:
                #  Free, locate allocation.
                if adr in memuse:
                    allocs = memuse[adr]
                    totalmem = totalmem - allocs[0]
                    if len(allocs) > 1:
                        memuse[adr] = allocs[1:]
                    else:
                        del memuse[adr]
                else:
                    #  Unmatched free, skip for now.
                    doprint = False
            if doprint:
                if totalmem > peak:
                    peak = totalmem
                print tic, label, totalmem/(1048576.0)

totals = {}
for adr in labels:
    #  If any remaining allocations.
    if adr in memuse:
        if labels[adr] in totals:
            totals[labels[adr]] = totals[labels[adr]] + memuse[adr][0]
        else:
            totals[labels[adr]] = memuse[adr][0]

print "# Aligned memory use by label"
total = 0.0
for label in sorted(totals):
    mem = totals[label]/(1048576.0)
    total = total + mem
    print "## ", label, '{:.3f}'.format(mem)
print "# Total aligned memory still in use : ", '{:.3f}'.format(total), " (MB)"
print "# Peak aligned memory usage         : ", '{:.3f}'.format(peak/1048576.0), " (MB)"
if process_use != "":
    print "#"
    print "# Memory use by process (all/system):", process_use
sys.exit(0)
