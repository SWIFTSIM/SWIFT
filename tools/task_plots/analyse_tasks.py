#!/usr/bin/env python
"""
Usage:
    analyse_tasks.py [options] input.dat

where input.dat is a thread info file for a step (MPI or non-MPI). Use the
'-y interval' flag of the swift and swift_mpi commands to create these
(you will also need to configure with the --enable-task-debugging option).

The output is an analysis of the task timings, including deadtime per thread
and step, total amount of time spent for each task type, for the whole step
and per thread and the minimum and maximum times spent per task type.

This file is part of SWIFT.
Copyright (c) 2017 Peter W. Draper (p.w.draper@durham.ac.uk)

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
parser = argparse.ArgumentParser(description="Analyse task dumps")

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
    "drift_spart",
    "drift_gpart",
    "drift_gpart_out",
    "hydro_end_force",
    "kick1",
    "kick2",
    "timestep",
    "timestep_limiter",
    "send",
    "recv",
    "grav_long_range",
    "grav_mm",
    "grav_down_in",
    "grav_down",
    "grav_mesh",
    "grav_end_force",
    "cooling",
    "star_formation",
    "logger",
    "stars_in",
    "stars_out",
    "stars_ghost_in",
    "stars_ghost",
    "stars_ghost_out",
    "stars_sort",
    "fof_self",
    "fof_pair",
    "count",
]

SUBTYPES = [
    "none",
    "density",
    "gradient",
    "force",
    "limiter",
    "grav",
    "external_grav",
    "tend",
    "xv",
    "rho",
    "gpart",
    "multipole",
    "spart",
    "stars_density",
    "stars_feedback",
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
if full_step.size == 13:
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
    updates = int(full_step[7])
    g_updates = int(full_step[8])
    s_updates = int(full_step[9])
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
    updates = int(full_step[6])
    g_updates = int(full_step[7])
    s_updates = int(full_step[8])

#  Get the CPU clock to convert ticks into milliseconds.
CPU_CLOCK = float(full_step[-1]) / 1000.0
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

maxthread = int(max(data[:, threadscol])) + 1
print("# Maximum thread id:", maxthread)

#  Avoid start and end times of zero.
sdata = data[data[:, ticcol] != 0]
sdata = data[data[:, toccol] != 0]

#  Now we process the required ranks.
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
        sid = int(data[line, -1])

        tasks[thread].append([tic, toc, tasktype, subtype, sid])

    #  Sort by tic and gather used threads.
    threadids = []
    for i in range(maxthread):
        tasks[i] = sorted(tasks[i], key=lambda task: task[0])
        threadids.append(i)

    #  Times per task.
    print("# Task times:")
    print("# -----------")
    print(
        "# {0:<17s}: {1:>7s} {2:>9s} {3:>9s} {4:>9s} {5:>9s} {6:>9s}".format(
            "type/subtype", "count", "minimum", "maximum", "sum", "mean", "percent"
        )
    )

    alltasktimes = {}
    sidtimes = {}
    for i in threadids:
        tasktimes = {}
        for task in tasks[i]:
            key = TASKTYPES[task[2]] + "/" + SUBTYPES[task[3]]
            dt = task[1] - task[0]
            if not key in tasktimes:
                tasktimes[key] = []
            tasktimes[key].append(dt)

            if not key in alltasktimes:
                alltasktimes[key] = []
            alltasktimes[key].append(dt)

            my_sid = task[4]
            if my_sid > -1:
                if not my_sid in sidtimes:
                    sidtimes[my_sid] = []
                sidtimes[my_sid].append(dt)

        print("# Thread : ", i)
        for key in sorted(tasktimes.keys()):
            taskmin = min(tasktimes[key])
            taskmax = max(tasktimes[key])
            tasksum = sum(tasktimes[key])
            print(
                "{0:19s}: {1:7d} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.4f} {6:9.2f}".format(
                    key,
                    len(tasktimes[key]),
                    taskmin,
                    taskmax,
                    tasksum,
                    tasksum / len(tasktimes[key]),
                    tasksum / total_t * 100.0,
                )
            )
        print()

    print("# All threads : ")
    for key in sorted(alltasktimes.keys()):
        taskmin = min(alltasktimes[key])
        taskmax = max(alltasktimes[key])
        tasksum = sum(alltasktimes[key])
        print(
            "{0:18s}: {1:7d} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.4f} {6:9.2f}".format(
                key,
                len(alltasktimes[key]),
                taskmin,
                taskmax,
                tasksum,
                tasksum / len(alltasktimes[key]),
                tasksum / (len(threadids) * total_t) * 100.0,
            )
        )
    print()

    # For pairs, show stuff sorted by SID
    print("# By SID (all threads): ")
    print(
        "# {0:<17s}: {1:>7s} {2:>9s} {3:>9s} {4:>9s} {5:>9s} {6:>9s}".format(
            "Pair/Sub-pair SID", "count", "minimum", "maximum", "sum", "mean", "percent"
        )
    )

    for sid in range(0, 13):
        if sid in sidtimes:
            sidmin = min(sidtimes[sid])
            sidmax = max(sidtimes[sid])
            sidsum = sum(sidtimes[sid])
            sidcount = len(sidtimes[sid])
            sidmean = sidsum / sidcount
        else:
            sidmin = 0.0
            sidmax = 0.0
            sidsum = 0.0
            sidcount = 0
            sidmean = 0.0
        print(
            "{0:3d} {1:15s}: {2:7d} {3:9.4f} {4:9.4f} {5:9.4f} {6:9.4f} {7:9.2f}".format(
                sid,
                SIDS[sid],
                sidcount,
                sidmin,
                sidmax,
                sidsum,
                sidmean,
                sidsum / (len(threadids) * total_t) * 100.0,
            )
        )
    print()

    #  Dead times.
    print("# Times not in tasks (deadtimes)")
    print("# ------------------------------")
    print("# Time before first task:")
    print("# no.    : {0:>9s} {1:>9s}".format("value", "percent"))
    predeadtimes = []
    for i in threadids:
        if len(tasks[i]) > 0:
            predeadtime = tasks[i][0][0]
            print(
                "thread {0:2d}: {1:9.4f} {2:9.4f}".format(
                    i, predeadtime, predeadtime / total_t * 100.0
                )
            )
            predeadtimes.append(predeadtime)
        else:
            predeadtimes.append(0.0)

    predeadmin = min(predeadtimes)
    predeadmax = max(predeadtimes)
    predeadsum = sum(predeadtimes)
    print(
        "#        : {0:>9s} {1:>9s} {2:>9s} {3:>9s} {4:>9s} {5:>9s}".format(
            "count", "minimum", "maximum", "sum", "mean", "percent"
        )
    )
    print(
        "all      : {0:9d} {1:9.4f} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.2f}".format(
            len(predeadtimes),
            predeadmin,
            predeadmax,
            predeadsum,
            predeadsum / len(predeadtimes),
            predeadsum / (len(threadids) * total_t) * 100.0,
        )
    )
    print()

    print("# Time after last task:")
    print("# no.    : {0:>9s} {1:>9s}".format("value", "percent"))
    postdeadtimes = []
    for i in threadids:
        if len(tasks[i]) > 0:
            postdeadtime = total_t - tasks[i][-1][1]
            print(
                "thread {0:2d}: {1:9.4f} {2:9.4f}".format(
                    i, postdeadtime, postdeadtime / total_t * 100.0
                )
            )
            postdeadtimes.append(postdeadtime)
        else:
            postdeadtimes.append(0.0)

    postdeadmin = min(postdeadtimes)
    postdeadmax = max(postdeadtimes)
    postdeadsum = sum(postdeadtimes)
    print(
        "#        : {0:>9s} {1:>9s} {2:>9s} {3:>9s} {4:>9s} {5:>9s}".format(
            "count", "minimum", "maximum", "sum", "mean", "percent"
        )
    )
    print(
        "all      : {0:9d} {1:9.4f} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.2f}".format(
            len(postdeadtimes),
            postdeadmin,
            postdeadmax,
            postdeadsum,
            postdeadsum / len(postdeadtimes),
            postdeadsum / (len(threadids) * total_t) * 100.0,
        )
    )
    print()

    #  Time in engine, i.e. from first to last tasks.
    print("# Time between tasks (engine deadtime):")
    print(
        "# no.    : {0:>9s} {1:>9s} {2:>9s} {3:>9s} {4:>9s} {5:>9s}".format(
            "count", "minimum", "maximum", "sum", "mean", "percent"
        )
    )
    enginedeadtimes = []
    for i in threadids:
        deadtimes = []
        if len(tasks[i]) > 0:
            last = tasks[i][0][0]
        else:
            last = 0.0
        for task in tasks[i]:
            dt = task[0] - last
            deadtimes.append(dt)
            last = task[1]

        #  Drop first value, last value already gone.
        if len(deadtimes) > 1:
            deadtimes = deadtimes[1:]
        else:
            #  Only one or fewer tasks, so no deadtime by definition.
            deadtimes = [0.0]

        deadmin = min(deadtimes)
        deadmax = max(deadtimes)
        deadsum = sum(deadtimes)
        print(
            "thread {0:2d}: {1:9d} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.4f} {6:9.2f}".format(
                i,
                len(deadtimes),
                deadmin,
                deadmax,
                deadsum,
                deadsum / len(deadtimes),
                deadsum / total_t * 100.0,
            )
        )
        enginedeadtimes.extend(deadtimes)

    deadmin = min(enginedeadtimes)
    deadmax = max(enginedeadtimes)
    deadsum = sum(enginedeadtimes)
    print(
        "all      : {0:9d} {1:9.4f} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.2f}".format(
            len(enginedeadtimes),
            deadmin,
            deadmax,
            deadsum,
            deadsum / len(enginedeadtimes),
            deadsum / (len(threadids) * total_t) * 100.0,
        )
    )
    print()

    #  All times in step.
    print("# All deadtimes:")
    print(
        "# no.    : {0:>9s} {1:>9s} {2:>9s} {3:>9s} {4:>9s} {5:>9s}".format(
            "count", "minimum", "maximum", "sum", "mean", "percent"
        )
    )
    alldeadtimes = []
    for i in threadids:
        deadtimes = []
        last = 0
        for task in tasks[i]:
            dt = task[0] - last
            deadtimes.append(dt)
            last = task[1]
        dt = total_t - last
        deadtimes.append(dt)

        deadmin = min(deadtimes)
        deadmax = max(deadtimes)
        deadsum = sum(deadtimes)
        print(
            "thread {0:2d}: {1:9d} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.4f} {6:9.2f}".format(
                i,
                len(deadtimes),
                deadmin,
                deadmax,
                deadsum,
                deadsum / len(deadtimes),
                deadsum / total_t * 100.0,
            )
        )
        alldeadtimes.extend(deadtimes)

    deadmin = min(alldeadtimes)
    deadmax = max(alldeadtimes)
    deadsum = sum(alldeadtimes)
    print(
        "all      : {0:9d} {1:9.4f} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.2f}".format(
            len(alldeadtimes),
            deadmin,
            deadmax,
            deadsum,
            deadsum / len(alldeadtimes),
            deadsum / (len(threadids) * total_t) * 100.0,
        )
    )
    print()

sys.exit(0)
