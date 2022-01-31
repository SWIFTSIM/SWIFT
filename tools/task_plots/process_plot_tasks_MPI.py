#!/usr/bin/env python3
"""
Usage:
    process_plot_tasks_MPI.py TIME-RANGE
       [--files FILES] [--weights] [--nproc NPROC]

where NPROC is the number of parallel processes to use to process task info
files and TIME-RANGE is the range (in ms) for the horizontal axis in the task
plots (a value of 0 means the task data is used to determine the range).

This script will process task info input files named
"thread_info_MPI-step<n>.dat" that can be produced by configuring SWIFT with
'--enable-task-debugging' and running SWIFT with the '-y INTERVAL' option. The
script produces the following output:
 - task plots ("step<n>r<r>.png") for each input file and rank
 - "index.html": an overview of all task plots
 - step pages ("step<n>r.html") with an overview of the task plots per rank for
   a single step. These can be opened by clicking on a task plot on the
   overview page.
 - step rank pages ("step<n>r<r>.html") with statistical information about the
   tasks for a specific step and a specific rank. These can be opened by
   clicking on a task plot on the step page for the same step.

By default, all thread_info_MPI*.dat files in the current working directory are
processed. The optional FILES argument allows more fine-grained control over
which files get included. Note that this script acts as a wrapper for
'plot_tasks.py' and 'analyse_tasks.py'; these scripts still expect some
files to be present in the current working directory.

The optional argument --weights will process the task info files in reverse
order of their size. This significantly improves load-balancing when using
a large number of processes for an inhomogeneous set of input files, but can
also lead to a large memory usage, since NPROC large files will be loaded into
memory simultaneously.

Note that this script is a Python version of an earlier bash script by Peter
Draper.

This file is part of SWIFT.

Copyright (C) 2015 Peter W. Draper (p.w.draper@durham.ac.uk)
          (C) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
All Rights Reserved.

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

import subprocess
import argparse
import glob
import re
import sys
import os
import numpy as np

# location of this script
# used to find other scripts we need to run
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

# parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument(
    "time_range", type=float, help="Time range to use for the horizontal axis (in ms)."
)
argparser.add_argument(
    "--files",
    "-f",
    default=None,
    nargs="+",
    help="Files to process (default: all thread_info_MPI-step*.dat files in the current directory).",
)
argparser.add_argument(
    "--weights",
    "-w",
    action="store_true",
    help="Use file sizes as weights to determine task order.",
)
argparser.add_argument(
    "--nproc", "-j", default=1, type=int, help="Number of parallel processes to use."
)
args = argparser.parse_args()

nproc = args.nproc
trange = args.time_range

# function used to extract the step counter from a thread_info_MPI-step*.dat file
def getcount(filename):
    return int(re.findall("\d+", filename)[0])


# function that parses a thread_info_MPI-step*.dat file and gets the number of ranks
def get_nrank(filename):
    # read the rank column from the file
    data = np.loadtxt(filename, usecols=[0], dtype=np.int32)
    # find the maximum
    # the number of ranks is this value + 1
    nrank = data.max() + 1
    return nrank


# sort the files based on the step number
# also do this if a list of files was provided
files = args.files
if files is None:
    files = sorted(glob.glob("thread_info_MPI-step*.dat"), key=getcount)
else:
    files = sorted(files, key=getcount)

# create a list of step numbers matching the files
nfile = len(files)
steps = np.zeros(nfile, dtype=np.int32)
for ifile in range(nfile):
    steps[ifile] = getcount(files[ifile])

nrank = get_nrank(files[0])

# create all the analysis commands we want to run without actually running them
cmds = []
weights = []
for ifile in range(nfile):
    outname = "step{0}r".format(steps[ifile])
    weight = os.path.getsize(files[ifile])

    cmd = "{0}/plot_tasks.py --expand 1 --limit {1} {2} {3}".format(
        script_dir, trange, files[ifile], outname
    )
    cmds.append(cmd)
    # plot_task commands are more expensive than their analyse_tasks
    # counterpart because they also need to write large image files
    weights.append(2 * weight)

    for irank in range(nrank):
        outname = "step{0}r{1}.stats".format(steps[ifile], irank)
        cmd = "{0}/analyse_tasks.py -r {1} --html {2} > {3}".format(
            script_dir, irank, files[ifile], outname
        )
        cmds.append(cmd)
        weights.append(weight)

if args.weights:
    # sort the commands according to their weight
    # long/expensive commands will be launched first, to achieve maximum overlap
    # with shorter commands
    weights = np.array(weights)
    order = np.argsort(weights)[::-1]
    cmds = np.array(cmds)[order]

# now run all commands in parallel using the requested number of processes

print("Done generating analysis commands, running them in parallel...")

# first, send off 'nproc' processes
icmd = 0
dfs = [None] * nproc
while icmd < len(cmds) and icmd < nproc:
    cmd = cmds[icmd]
    print("Starting {0} in slot {1}".format(cmd, icmd))
    dfs[icmd] = subprocess.Popen(cmd, shell=True)
    icmd += 1

# now keep spawning more processes until all commands have been submitted
while icmd < len(cmds):
    # loop over the running processes
    for iproc in range(nproc):
        # if a process finished, replace it with a new command
        if not dfs[iproc].poll() is None:
            print("Slot {0} finished".format(iproc))
            cmd = cmds[icmd]
            print("Starting {0} in slot {1}".format(cmd, iproc))
            dfs[iproc] = subprocess.Popen(cmd, shell=True)
            icmd += 1
            # we are out of commands, exit the for-loop (and the while-loop)
            if icmd == len(cmds):
                break

# now wait for the remaining processes
# we don't care that the processes will not finish in list order, since we
# have to wait for the slowest one anyway
for iproc in range(nproc):
    if not dfs[iproc] is None:
        dfs[iproc].wait()
        print("Slot {0} finished".format(iproc))

print("Done processing files. Creating web pages...")

# generate web pages
# note that we don't bother adding newlines, since those are ignored by the
# browser anyway
htmltag = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">'

index_file = open("index.html", "w")
index_file.write(htmltag)
index_file.write("<html><head><title>SWIFT task graphs</title></head><body>")
index_file.write("<h1>SWIFT task graphs</h1>")

for ifile in range(nfile):
    step = steps[ifile]
    index_file.write("<h2>Step {0}</h2>".format(step))
    index_file.write('<ul style="list-style-type:none"><li>')

    step_file = open("step{0}r.html".format(step), "w")
    step_file.write(htmltag)
    step_file.write("<html><body>")

    for irank in range(nrank):
        index_file.write('<a href="step{0}r.html">'.format(step))
        index_file.write(
            '<img src="step{0}r{1}.png" width="400px"/></a>'.format(step, irank)
        )

        step_file.write('<a href="step{0}r{1}.html">'.format(step, irank))
        step_file.write('<img src="step{0}r{1}.png"/></a>'.format(step, irank))

        with open("step{0}r{1}.html".format(step, irank), "w") as rank_file:
            rank_file.write(htmltag)
            rank_file.write(
                '<html></body><img src="step{0}r{1}.png"/>'.format(step, irank)
            )
            rank_file.write('<pre><nav>Jump to: <a href="#all">all threads</a> ')
            rank_file.write('<a href="#dead">dead times</a></nav>\n')
            with open("step{0}r{1}.stats".format(step, irank), "r") as stats_file:
                rank_file.write(stats_file.read())
            rank_file.write("</pre></body></html>")

    step_file.write("</body></html>")

    index_file.write("</li></ul>")

index_file.write("</body></html>")

print("Done.")
