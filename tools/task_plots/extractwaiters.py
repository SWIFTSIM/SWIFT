#!/usr/bin/env python
"""
Usage:
    extractwaiters.py [options] time input.dat

where input.dat is a thread info file for a step (MPI or non-MPI). Use the
'-y interval' flag of the swift and swift_mpi commands to create these
(you will also need to configure with the --enable-task-debugging option).

The output is a list of all the tasks that are waiting to run at the
given offset time (in ms, not ticks).

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

import pylab as pl
import sys
import argparse

#  Handle the command line.
parser = argparse.ArgumentParser(description="Extract instantaneous waiters from task dumps")

parser.add_argument("time",
                    type=float,
                    help="Time of interest (ms from start of step)")
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
parser.add_argument(
    "-t",
    "--todo",
    dest="todo",
    help="Extract all tasks still to run instead (default: False)",
    default=False,
    action="store_true",
)

args = parser.parse_args()
infile = args.input
ourtime = args.time

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
    "drift_bpart",
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
    "star_formation_in",
    "star_formation_out",
    "logger",
    "stars_in",
    "stars_out",
    "stars_ghost_in",
    "stars_ghost",
    "stars_ghost_out",
    "stars_sort",
    "stars_resort",
    "bh_in",
    "bh_out",
    "bh_ghost",
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
    "tend_part",
    "tend_gpart",
    "tend_spart",
    "tend_bpart",
    "xv",
    "rho",
    "gpart",
    "multipole",
    "spart",
    "stars_density",
    "stars_feedback",
    "sf_counts"
    "bpart",
    "bh_density",
    "bh_swallow",
    "do_swallow",
    "bh_feedback",
    "count",
]

#  Read input.
data = pl.loadtxt(infile)
full_step = data[0, :]

#  Do we have an MPI file?
full_step = data[0, :]
if full_step.size == 14:
    print "## MPI mode"
    mpimode = True
    nranks = int(max(data[:, 0])) + 1
    print "## Number of ranks:", nranks
    rankcol = 0
    threadscol = 1
    taskcol = 2
    subtaskcol = 3
    ticcol = 5
    toccol = 6
    qticcol = 7
    updates = int(full_step[7])
    g_updates = int(full_step[8])
    s_updates = int(full_step[9])
else:
    print "## non MPI mode"
    nranks = 1
    mpimode = False
    rankcol = -1
    threadscol = 0
    taskcol = 1
    subtaskcol = 2
    ticcol = 4
    toccol = 5
    qticcol = 6
    updates = int(full_step[6])
    g_updates = int(full_step[7])
    s_updates = int(full_step[8])

#  Get the CPU clock to convert ticks into milliseconds.
CPU_CLOCK = float(full_step[-1]) / 1000.0
if args.verbose:
    print "## CPU frequency:", CPU_CLOCK * 1000.0
print "##   updates:", updates
print "## g_updates:", g_updates
print "## s_updates:", s_updates

if mpimode:
    if args.rank == "all":
        ranks = list(range(nranks))
    else:
        ranks = [int(args.rank)]
        if ranks[0] >= nranks:
            print "Error: maximum rank is " + str(nranks - 1)
            sys.exit(1)
else:
    ranks = [1]

maxthread = int(max(data[:, threadscol])) + 1
print "## Maximum thread id:", maxthread

#  Avoid start and end times of zero.
sdata = data[data[:, ticcol] != 0]
sdata = data[data[:, toccol] != 0]

#  Now we process the required ranks.
for rank in ranks:
    if mpimode:
        print "## Rank", rank
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
    print "## Data range: ", total_t, "ms"
    print

    #  Correct times to relative values.
    start_t = float(tic_step)
    data[:, ticcol] -= start_t
    data[:, toccol] -= start_t
    data[:, qticcol] -= start_t
    end_t = (toc_step - start_t) / CPU_CLOCK

    #  Extract tasks that are queued, but not started or tasks that are still
    #  to do, that is run or complete after ourtime.
    print "# rank thread task subtask tic toc qtic time qtime"
    tasks = []
    num_lines = pl.shape(data)[0]
    for line in range(num_lines):
        thread = int(data[line, threadscol])
        tic = int(data[line, ticcol]) / CPU_CLOCK
        toc = int(data[line, toccol]) / CPU_CLOCK
        qtic = int(data[line, qticcol]) / CPU_CLOCK
        if (ourtime >= qtic and ourtime <= tic) or (args.todo and ourtime < toc):
            tasktype = TASKTYPES[int(data[line, taskcol])]
            subtype = SUBTYPES[int(data[line, subtaskcol])]
            print rank, thread, tasktype, subtype, tic, toc, qtic, (toc - tic), (tic - qtic)

sys.exit(0)
