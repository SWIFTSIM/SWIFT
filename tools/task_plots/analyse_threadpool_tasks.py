#!/usr/bin/env python
"""
Usage:
    analsyse_threadpool_tasks.py [options] input.dat

where input.dat is a threadpool dump for a step.  Use the '-Y interval' flag
of the swift command to create these.

The output is an analysis of the threadpool task timings, including deadtime
per thread and step, total amount of time spent for each task type, for the
whole step and per thread and the minimum and maximum times spent per task
type.

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

parser.add_argument("input", help="Threadpool data file (-y output)")
parser.add_argument(
    "-v",
    "--verbose",
    dest="verbose",
    help="Verbose output (default: False)",
    default=False,
    action="store_true",
)
parser.add_argument(
    "--html",
    dest="html",
    help="Use html titles and anchors in the output (default: False)",
    default=False,
    action="store_true",
)

args = parser.parse_args()
infile = args.input
with_html = args.html

#  Read header. First two lines.
with open(infile) as infid:
    head = [next(infid) for x in range(2)]
header = head[1][2:].strip()
header = eval(header)
nthread = int(header["num_threads"]) + 1
CPU_CLOCK = float(header["cpufreq"]) / 1000.0
print("Number of threads: ", nthread - 1)
if args.verbose:
    print("CPU frequency:", CPU_CLOCK * 1000.0)

#  Read input.
data = pl.genfromtxt(infile, dtype=None, delimiter=" ", encoding=None)

#  Mixed types, so need to separate.
tics = []
tocs = []
funcs = []
threads = []
chunks = []
for i in data:
    if i[0] != "#":
        funcs.append(i[0].replace("_mapper", ""))
        if i[1] < 0:
            threads.append(nthread - 1)
        else:
            threads.append(i[1])
        chunks.append(i[2])
        tics.append(i[3])
        tocs.append(i[4])
tics = pl.array(tics)
tocs = pl.array(tocs)
funcs = pl.array(funcs)
threads = pl.array(threads)
chunks = pl.array(chunks)

#  Recover the start and end time
tic_step = min(tics)
toc_step = max(tocs)

#  Calculate the time range.
total_t = (toc_step - tic_step) / CPU_CLOCK
print("# Data range: ", total_t, "ms")
print()

#  Correct times to relative millisecs.
start_t = float(tic_step)
tics = (tics - start_t) / CPU_CLOCK
tocs = (tocs - start_t) / CPU_CLOCK

tasks = {}
tasks[-1] = []
for i in range(nthread):
    tasks[i] = []

#  Gather into by thread data.
for i in range(len(tics)):
    tasks[threads[i]].append([tics[i], tocs[i], funcs[i]])

#  Don't actually process the fake thread.
nthread = nthread - 1

#  Sort by tic and gather used thread ids.
threadids = []
for i in range(nthread):
    if len(tasks[i]) > 0:
        tasks[i] = sorted(tasks[i], key=lambda task: task[0])
        threadids.append(i)

#  Times per task.
print("# Task times:")
print("# -----------")
print(
    "# {0:<31s}: {1:>7s} {2:>9s} {3:>9s} {4:>9s} {5:>9s} {6:>9s}".format(
        "type/subtype", "count", "minimum", "maximum", "sum", "mean", "percent"
    )
)
alltasktimes = {}
sidtimes = {}
for i in threadids:
    tasktimes = {}
    for task in tasks[i]:
        key = task[2]
        dt = task[1] - task[0]
        if not key in tasktimes:
            tasktimes[key] = []
        tasktimes[key].append(dt)

        if not key in alltasktimes:
            alltasktimes[key] = []
        alltasktimes[key].append(dt)

    if with_html:
        print('<div id="thread{}"></div>'.format(i))
    print("# Thread : ", i)
    for key in sorted(tasktimes.keys()):
        taskmin = min(tasktimes[key])
        taskmax = max(tasktimes[key])
        tasksum = sum(tasktimes[key])
        print(
            "{0:33s}: {1:7d} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.4f} {6:9.2f}".format(
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

if with_html:
    print('<div id="all"></div>')
print("# All threads : ")
for key in sorted(alltasktimes.keys()):
    taskmin = min(alltasktimes[key])
    taskmax = max(alltasktimes[key])
    tasksum = sum(alltasktimes[key])
    print(
        "{0:33s}: {1:7d} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.4f} {6:9.2f}".format(
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

#  Dead times.
print("# Times not in tasks (deadtimes)")
print("# ------------------------------")
if with_html:
    print('<div id="before"></div>')
print("# Time before first task:")
print("# no.    : {0:>9s} {1:>9s}".format("value", "percent"))
predeadtimes = []
for i in threadids:
    predeadtime = tasks[i][0][0]
    print(
        "thread {0:2d}: {1:9.4f} {2:9.4f}".format(
            i, predeadtime, predeadtime / total_t * 100.0
        )
    )
    predeadtimes.append(predeadtime)

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

if with_html:
    print('<div id="after"></div>')
print("# Time after last task:")
print("# no.    : {0:>9s} {1:>9s}".format("value", "percent"))
postdeadtimes = []
for i in threadids:
    postdeadtime = total_t - tasks[i][-1][1]
    print(
        "thread {0:2d}: {1:9.4f} {2:9.4f}".format(
            i, postdeadtime, postdeadtime / total_t * 100.0
        )
    )
    postdeadtimes.append(postdeadtime)

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

#  Time in threadpool, i.e. from first to last tasks.
if with_html:
    print('<div id="between"></div>')
print("# Time between tasks (threadpool deadtime):")
print(
    "# no.    : {0:>9s} {1:>9s} {2:>9s} {3:>9s} {4:>9s} {5:>9s}".format(
        "count", "minimum", "maximum", "sum", "mean", "percent"
    )
)
threadpooldeadtimes = []
for i in threadids:
    deadtimes = []
    last = tasks[i][0][0]
    for task in tasks[i]:
        dt = task[0] - last
        deadtimes.append(dt)
        last = task[1]

    #  Drop first value, last value already gone.
    if len(deadtimes) > 1:
        deadtimes = deadtimes[1:]
    else:
        #  Only one task, so no deadtime by definition.
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
    threadpooldeadtimes.extend(deadtimes)

deadmin = min(threadpooldeadtimes)
deadmax = max(threadpooldeadtimes)
deadsum = sum(threadpooldeadtimes)
print(
    "all      : {0:9d} {1:9.4f} {2:9.4f} {3:9.4f} {4:9.4f} {5:9.2f}".format(
        len(threadpooldeadtimes),
        deadmin,
        deadmax,
        deadsum,
        deadsum / len(threadpooldeadtimes),
        deadsum / (len(threadids) * total_t) * 100.0,
    )
)
print()

#  All times in step.
if with_html:
    print('<div id="dead"></div>')
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
