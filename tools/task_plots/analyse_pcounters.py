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
parser = argparse.ArgumentParser(description="Analyse task dumps for performance counters")

parser.add_argument("input", help="Thread data file (-y output)")
parser.add_argument(
    "-v",
    "--verbose",
    dest="verbose",
    help="Verbose output (default: False)",
    default=False,
    action="store_true",
)
parser.add_argument(
    "-r",
    "--rank",
    dest="rank",
    help="Rank to process (default: all)",
    default="all",
    action="store",
)

args = parser.parse_args()
infile = args.input

#  Tasks and subtypes. Indexed as in tasks.h.
TASKTYPES = [
    "none",
    "sort",
    "self",
    "pair",
    "sub_self",
    "sub_pair",
    "init_grav",
    "init_grav_out",
    "ghost_in",
    "ghost",
    "ghost_out",
    "extra_ghost",
    "drift_part",
    "drift_gpart",
    "drift_gpart_out",
    "end_force",
    "kick1",
    "kick2",
    "timestep",
    "send",
    "recv",
    "grav_long_range",
    "grav_mm",
    "grav_down_in",
    "grav_down",
    "grav_mesh",
    "cooling",
    "star_formation",
    "sourceterms",
    "stars_ghost_in",
    "stars_ghost",
    "stars_ghost_out",
    "count",
]

SUBTYPES = [
    "none",
    "density",
    "gradient",
    "force",
    "grav",
    "external_grav",
    "tend",
    "xv",
    "rho",
    "gpart",
    "multipole",
    "spart",
    "stars_density",
    "count",
]

SIDS = [
    "(-1,-1,-1)",
    "(-1,-1, 0)",
    "(-1,-1, 1)",
    "(-1, 0,-1)",
    "(-1, 0, 0)",
    "(-1, 0, 1)",
    "(-1, 1,-1)",
    "(-1, 1, 0)",
    "(-1, 1, 1)",
    "( 0,-1,-1)",
    "( 0,-1, 0)",
    "( 0,-1, 1)",
    "( 0, 0,-1)",
]

#  Read input.
data = pl.loadtxt(infile)
full_step = data[0, :]

#  Do we have an MPI file?
full_step = data[0, :]
if full_step.size == 17:
    print("# MPI mode")
    mpimode = True
    nranks = int(max(data[:, 0])) + 1
    print("# Number of ranks:", nranks)
    rankcol = 0
    threadscol = 1
    taskcol = 2
    subtaskcol = 3
    ticcol = 5
    toccol = 6
    cpuidcol = 13
    local_hitscol = 14
    local_misscol = 15
    instr_retcol = 16
    updates = int(full_step[7])
    g_updates = int(full_step[8])
    s_updates = int(full_step[9])
    CPU_CLOCK = float(full_step[12]) / 1000.0
else:
    print("# non MPI mode")
    nranks = 1
    mpimode = False
    rankcol = -1
    threadscol = 0
    taskcol = 1
    subtaskcol = 2
    ticcol = 4
    toccol = 5
    cpuidcol = 11
    local_hitscol = 12
    local_misscol = 13
    instr_retcol = 14
    updates = int(full_step[6])
    g_updates = int(full_step[7])
    s_updates = int(full_step[8])
    CPU_CLOCK = float(full_step[10]) / 1000.0

#  Get the CPU clock to convert ticks into milliseconds.
if args.verbose:
    print("# CPU frequency:", CPU_CLOCK * 1000.0)
print("#   updates:", updates)
print("# g_updates:", g_updates)
print("# s_updates:", s_updates)

if mpimode:
    if args.rank == "all":
        ranks = list(range(nranks))
    else:
        ranks = [int(args.rank)]
        if ranks[0] >= nranks:
            print("Error: maximum rank is " + str(nranks - 1))
            sys.exit(1)
else:
    ranks = [1]

maxthread = int(max(data[:,threadscol])) + 1
print("# Maximum thread id:", maxthread)

#  Avoid start and end times of zero.
sdata = data[data[:, ticcol] != 0]
sdata = data[data[:, toccol] != 0]

for rank in ranks:
    if mpimode:
        print("# Rank", rank)
        data = sdata[sdata[:, rankcol] == rank]
        full_step = data[0, :]
    else:
        data = sdata

    #  Recover the start and end time
    tic_step = int(full_step[ticcol])
    toc_step = int(full_step[toccol])
    data = data[1:, :]

    #  Avoid start and end times of zero.
    data = data[data[:, ticcol] != 0]
    data = data[data[:, toccol] != 0]

    #  Calculate the time range.
    total_t = (toc_step - tic_step) / CPU_CLOCK
    print("# Data range: ", total_t, "ms")
    print()

    #  Correct times to relative values.
    start_t = float(tic_step)
    data[:, ticcol] -= start_t
    data[:, toccol] -= start_t
    end_t = (toc_step - start_t) / CPU_CLOCK

    tasks = {}
    tasks[-1] = []
    for i in range(maxthread):
        tasks[i] = []

    #  Gather into by thread data.
    num_lines = pl.shape(data)[0]
    for line in range(num_lines):
        thread = int(data[line, threadscol])
        tic = int(data[line, ticcol]) / CPU_CLOCK
        toc = int(data[line, toccol]) / CPU_CLOCK
        tasktype = int(data[line, taskcol])
        subtype = int(data[line, subtaskcol])
        cpuid = int(data[line,cpuidcol])
        local_hits = int(data[line,local_hitscol])
        local_miss = int(data[line,local_misscol])
        instr_ret = int(data[line,instr_retcol])

        tasks[thread].append([tic,toc,tasktype,subtype,cpuid,local_hits,local_miss,instr_ret])

    threadids = []
    for i in range(maxthread):
        tasks[i] = sorted(tasks[i], key=lambda task: task[0])
        threadids.append(i)

    print("# Task times:")
    print("# -----------")
    print("# {0:<17s}: {1:>7s} {2:>9s} {3:>9s} {4:>9s} {5:>15s} {6:>15s}".format(
            "type/subtype","count","minimum","maximum","%hits","min_inst","max_inst")
         )

    all_tasktimes={}
    all_localphits={}
    all_instructions={}
    for i in threadids:
        tasktimes = {}
        localphits = {}
        instructions = {}
        for task in tasks[i]:
            key = TASKTYPES[task[2]] + "/" + SUBTYPES[task[3]]
            dt = task[1] - task[0]
            if(task[5] > 0):
                phits = task[5] / (task[5]+task[6]) *100.0
            instruct = task[7]
            if not key in tasktimes:
                tasktimes[key] = []
                localphits[key] = []
                instructions[key] = []
            tasktimes[key].append(dt)
            if(task[5] > 0):
                localphits[key].append(phits)
            instructions[key].append(instruct)

            if not key in all_tasktimes:
                all_tasktimes[key] = []
                all_localphits[key] = []
                all_instructions[key] = []
            all_tasktimes[key].append(dt)
            all_localphits[key].append(phits)
            all_instructions[key].append(instruct)

        print("# Thread : ", i, " CPU: ", tasks[i][0][4])
        for key in sorted(tasktimes.keys()):
            taskmin = min(tasktimes[key])
            taskmax = max(tasktimes[key])
            instmin = min(instructions[key])
            instmax = max(instructions[key])
            pctsum = sum(localphits[key])
            print("{0:19s}: {1:7d} {2:9.4f} {3:9.4f} {4:9.5f} {5:15d} {6:15d}".format(
                   key,
                   len(tasktimes[key]),
                   taskmin,
                   taskmax,
                   pctsum / len(tasktimes[key]),
                   instmin,
                   instmax)
                 )
        print()
        print("# All threads : ")
        for key in sorted(all_tasktimes.keys()):
            taskmin = min(all_tasktimes[key])
            taskmax = max(all_tasktimes[key])
            instmin = min(all_instructions[key])
            instmax = max(all_instructions[key])
            pctsum = sum(all_localphits[key])
            print("{0:19s}: {1:7d} {2:9.4f} {3:9.4f} {4:9.5f} {5:15d} {6:15d}".format(
                   key,
                   len(all_tasktimes[key]),
                   taskmin,
                   taskmax,
                   pctsum / len(all_tasktimes[key]),
                   instmin,
                   instmax)
                 )



