#!/usr/bin/env python3
"""A module that defines a generic task parser class.

This module should not be used directly, but instead should be imported by other
modules which need to parse the thread_info* files.

This can be used by any of the analysis scripts which need to parse the
thread_info*files.

This file is part of SWIFT.

Copyright (C) 2015 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
                   Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
                   Matthieu Schaller (schaller@strw.leidenuniv.nl)
          (C) 2017 Peter W. Draper (p.w.draper@durham.ac.uk)
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
import sys
import numpy as np

# import hardcoded data
from swift_hardcoded_data import TASKTYPES, SUBTYPES


# Cell types.
CELLTYPES = [
    "Regular",
    "Zoom",
    "Buff",
    "Bkg",
]

# Cell subtypes.
CELLSUBTYPES = [
    "Regular",
    "Neighbour",
    "Void",
    # Empty will never appear
]


class Task:
    def __init__(
        self,
        rank,
        thread,
        type_int,
        subtype_int,
        tic,
        toc,
        ci_part_count,
        cj_part_count,
        ci_gpart_count,
        cj_gpart_count,
        ci_type,
        cj_type,
        ci_subtype,
        cj_subtype,
        ci_depth,
        cj_depth,
        min_dist,
        mpole_dist,
    ):
        self.rank = rank
        self.thread = thread
        self.type = TASKTYPES[type_int]
        self.subtype = SUBTYPES[subtype_int]
        self.task = (
            self.type + "/" + self.subtype
            if self.subtype != "none"
            else self.type
        )
        self.tic = tic
        self.toc = toc
        self.ci_part_count = ci_part_count
        self.cj_part_count = cj_part_count
        self.ci_gpart_count = ci_gpart_count
        self.cj_gpart_count = cj_gpart_count
        self.ci_type = CELLTYPES[ci_type]
        self.cj_type = CELLTYPES[cj_type]
        self.ci_subtype = CELLSUBTYPES[ci_subtype]
        self.cj_subtype = CELLSUBTYPES[cj_subtype]
        self.ci_depth = ci_depth
        self.cj_depth = cj_depth
        self.min_dist = min_dist
        self.mpole_dist = mpole_dist
        self.dt = toc - tic

    def __str__(self):
        return (
            "Rank: %d, Thread: %d, Task: %d, Subtask: %d, Tic: %d, Toc: %d, "
            "ci_type: %d, cj_type: %d, ci_subtype: %d, cj_subtype: %d, "
            "ci_depth: %d, cj_depth: %d"
            % (
                self.rank,
                self.thread,
                self.task,
                self.subtask,
                self.tic,
                self.toc,
                self.ci_type,
                self.cj_type,
                self.ci_subtype,
                self.cj_subtype,
                self.ci_depth,
                self.cj_depth,
            )
        )


class TaskParser:
    def __init__(
        self, filename, ranks=None, verbose=True, delta_t=0.0, mintic=-1
    ):
        # Define the filename
        self.filename = filename

        # Are we talking?
        self.verbose = verbose

        # Define the user set duration. (If not set this will be derived later)
        self.delta_t = delta_t

        # Ensure what we've been given is sensible
        if self.delta_t < 0.0:
            print("The time-range must be >=0!")
            sys.exit(1)

        # Define the user set start time. (If not set this will be extracted
        # from the data)
        self.start_t = mintic

        # Define the attribute to hold the end tick
        self.end_t = 0

        # Load the data
        self.data = None
        self.full_step = None  # header containing the full step information
        self._load_data()

        # Flag for MPI mode
        self.mpimode = "MPI" in filename

        # Define a list of all the ranks we have
        self.ranks = ranks

        # Define the column look up table and populate it
        self._col_look_up = {}
        self._define_columns()

        # Process the header
        self.cpu_clock = None
        self.nthread = None
        self._process_header(mintic)

        # Clean up the data
        self._clean_up_data()

        # How many tasks are there total?
        self.ntasks = self.data[:, 0].size

        # Parse the file populating the look up table
        self.tasks = None
        self.dt = None
        self._parse_tasks()

    def _load_data(self):
        self.data = np.loadtxt(self.filename)
        self.full_step = self.data[0, :]

    def _define_columns(self):
        # If we have been given a subset of ranks, parse them
        if self.ranks is not None:
            self.ranks = [int(item) for item in self.ranks.split(",")]

        #  Do we have an MPI file?
        if self.mpimode:
            print("# MPI mode")
            if self.ranks is None:
                self.ranks = list(range(int(max(self.data[:, 0])) + 1))
            print("# Number of ranks:", len(self.ranks))
            self._col_look_up["rank"] = 0
            self._col_look_up["threads"] = 1
            self._col_look_up["task"] = 2
            self._col_look_up["subtask"] = 3
            self._col_look_up["tic"] = 5
            self._col_look_up["toc"] = 6
            self._col_look_up["ci_part_count"] = 7
            self._col_look_up["cj_part_count"] = 8
            self._col_look_up["ci_gpart_count"] = 9
            self._col_look_up["cj_gpart_count"] = 10
            self._col_look_up["ci_type"] = 13
            self._col_look_up["cj_type"] = 14
            self._col_look_up["ci_subtype"] = 15
            self._col_look_up["cj_subtype"] = 16
            self._col_look_up["ci_depth"] = 17
            self._col_look_up["cj_depth"] = 18
            self._col_look_up["min_dist"] = 19
            self._col_look_up["mpole_dist"] = 20
        else:
            print("# non MPI mode")
            self.ranks = [0]
            self._col_look_up["rank"] = -1
            self._col_look_up["threads"] = 0
            self._col_look_up["task"] = 1
            self._col_look_up["subtask"] = 2
            self._col_look_up["tic"] = 4
            self._col_look_up["toc"] = 5
            self._col_look_up["ci_part_count"] = 6
            self._col_look_up["cj_part_count"] = 7
            self._col_look_up["ci_gpart_count"] = 8
            self._col_look_up["cj_gpart_count"] = 9
            self._col_look_up["ci_type"] = 11
            self._col_look_up["cj_type"] = 12
            self._col_look_up["ci_subtype"] = 13
            self._col_look_up["cj_subtype"] = 14
            self._col_look_up["ci_depth"] = 15
            self._col_look_up["cj_depth"] = 16
            self._col_look_up["min_dist"] = 17
            self._col_look_up["mpole_dist"] = 18

    def _extract_column(self, column):
        return self.data[:, self._col_look_up[column]]

    def _process_header(self, mintic):
        # Extract the CPU clock
        if self.mpimode:
            self.cpu_clock = float(self.full_step[12]) / 1000.0
        else:
            self.cpu_clock = float(self.full_step[10]) / 1000.0
        if self.verbose:
            print("# CPU frequency:", self.cpu_clock * 1000.0)

        # Count the number of threads
        self.nthread = int(max(self.data[:, self._col_look_up["threads"]])) + 1
        print("# Number of threads:", self.nthread)

        # Each rank can have different clocks (compute node), but we want to
        # use the same delta times range for comparisons, so we suck it up and
        # take the hit of precalculating this, unless the user knows better.
        self.delta_t = self.delta_t * self.cpu_clock
        if self.delta_t == 0:
            for rank in self.ranks:
                if self.mpimode:
                    data = self.data[self.task_ranks == rank]
                else:
                    data = self.data

                # Get a local version of the full step.
                full_step = data[0, :]

                #  Start and end times for this rank. Can be changed using the
                #  mintic option. This moves our zero time to other time.
                #  Useful for comparing to other plots.
                if mintic < 0:
                    tic_step = int(full_step[self._col_look_up["tic"]])
                else:
                    tic_step = mintic
                toc_step = int(full_step[self._col_look_up["toc"]])
                dt = toc_step - tic_step
                if dt > self.delta_t:
                    self.delta_t = dt
        print("# Data range: ", self.delta_t / self.cpu_clock, "ms")

        # Get the start tic
        if self.start_t < 0:
            self.start_t = np.inf
            for rank in self.ranks:
                if self.mpimode:
                    data = self.data[self.task_ranks == rank]
                else:
                    data = self.data

                # Get a local version of the full step.
                full_step = data[0, :]

                # Get the start tic for this rank
                tic_step = int(full_step[self._col_look_up["tic"]])
                if self.start_t > tic_step:
                    self.start_t = tic_step

        # Set the end toc
        self.end_t = self.start_t + self.delta_t

    def _clean_up_data(self):
        # Remove the header
        self.data = self.data[1:, :]

        # Remove start and end times of zero
        mask = np.logical_and(
            self.tics != 0,
            self.tocs != 0,
        )
        self.data = self.data[mask, :]

        # Make tics and tocs relative to start
        self.data[:, self._col_look_up["tic"]] -= self.start_t
        self.data[:, self._col_look_up["toc"]] -= self.start_t

        # Convert tics and tocs to ms
        self.data[:, self._col_look_up["tic"]] /= self.cpu_clock
        self.data[:, self._col_look_up["toc"]] /= self.cpu_clock

    def _parse_tasks(self):
        # Prepare the arrays we'll populate
        self.tasks = np.zeros(self.ntasks, dtype=object)
        self.task_labels = np.zeros(self.ntasks, dtype=object)
        self.dt = np.zeros(self.ntasks, dtype=np.float64)

        # Get local copies of the columns to avoid extracting every single loop
        task_ranks = self.task_ranks
        task_threads = self.task_threads
        task_types = self.task_types
        task_subtypes = self.task_subtypes
        tics = self.tics
        tocs = self.tocs
        ci_part_count = self.ci_part_count
        cj_part_count = self.cj_part_count
        ci_gpart_count = self.ci_gpart_count
        cj_gpart_count = self.cj_gpart_count
        ci_types = self.ci_types
        cj_types = self.cj_types
        ci_subtypes = self.ci_subtypes
        cj_subtypes = self.cj_subtypes
        ci_depths = self.ci_depths
        cj_depths = self.cj_depths
        min_dists = self.min_dists
        mpole_dists = self.mpole_dists

        # Loop over tasks creating task objects
        for i in range(self.ntasks):
            self.tasks[i] = Task(
                task_ranks[i],
                task_threads[i],
                task_types[i],
                task_subtypes[i],
                tics[i],
                tocs[i],
                ci_part_count[i],
                cj_part_count[i],
                ci_gpart_count[i],
                cj_gpart_count[i],
                ci_types[i],
                cj_types[i],
                ci_subtypes[i],
                cj_subtypes[i],
                ci_depths[i],
                cj_depths[i],
                min_dists[i],
                mpole_dists[i],
            )
            self.task_labels[i] = self.tasks[i].task
            self.dt[i] = self.tasks[i].dt

    def _get_tasks_with_mask(self, mask=None):
        if mask is not None:
            tasks = self.tasks[mask]
        else:
            tasks = self.tasks
        unique_tasks = np.unique(self.task_labels[mask])
        unique_count = len(unique_tasks)
        return tasks, unique_tasks, unique_count

    def get_tasks(self):
        return self._get_tasks_with_mask()

    def get_tasks_on_rank(self, rank):
        return self._get_tasks_with_mask(mask=self.task_ranks == rank)

    def get_tasks_at_depth(self, depth):
        return self._get_tasks_with_mask(
            mask=np.logical_or(
                self.ci_depths == depth,
                self.cj_depths == depth,
            )
        )

    def get_tasks_tictoc_by_thread(
        self,
        task_type=None,
        cell_type=None,
        cell_subtype=None,
        depth=None,
    ):
        tasks = self.tasks
        mask = np.ones(len(tasks), dtype=bool)
        if task_type is not None:
            mask = np.logical_and(mask, self.task_types == task_type)
        if depth is not None:
            mask = np.logical_and(
                mask,
                np.logical_or(
                    self.ci_depths == depth, self.cj_depths == depth
                ),
            )
        if cell_type is not None:
            mask = np.logical_and(
                mask,
                np.logical_or(
                    self.ci_types == cell_type, self.cj_types == cell_type
                ),
            )
        if cell_subtype is not None:
            mask = np.logical_and(
                mask,
                np.logical_or(
                    self.ci_subtypes == cell_subtype,
                    self.cj_subtypes == cell_subtype,
                ),
            )
        labels = self.task_labels[mask]
        tics = self.tics[mask]
        tocs = self.tocs[mask]
        threads = self.task_threads[mask]
        unique_labels = set(labels)
        return labels, tics, tocs, threads, len(unique_labels)

    @property
    def task_ranks(self):
        return np.int32(self._extract_column("rank"))

    @property
    def task_threads(self):
        return np.int32(self._extract_column("threads"))

    @property
    def task_types(self):
        return np.int32(self._extract_column("task"))

    @property
    def task_subtypes(self):
        return np.int32(self._extract_column("subtask"))

    @property
    def tics(self):
        return np.float64(self._extract_column("tic"))

    @property
    def tocs(self):
        return np.float64(self._extract_column("toc"))

    @property
    def ci_part_count(self):
        return np.int32(self._extract_column("ci_part_count"))

    @property
    def cj_part_count(self):
        return np.int32(self._extract_column("cj_part_count"))

    @property
    def ci_gpart_count(self):
        return np.int32(self._extract_column("ci_gpart_count"))

    @property
    def cj_gpart_count(self):
        return np.int32(self._extract_column("cj_gpart_count"))

    @property
    def ci_types(self):
        return np.int32(self._extract_column("ci_type"))

    @property
    def cj_types(self):
        return np.int32(self._extract_column("cj_type"))

    @property
    def ci_subtypes(self):
        return np.int32(self._extract_column("ci_subtype"))

    @property
    def cj_subtypes(self):
        return np.int32(self._extract_column("cj_subtype"))

    @property
    def ci_depths(self):
        return np.int32(self._extract_column("ci_depth"))

    @property
    def cj_depths(self):
        return np.int32(self._extract_column("cj_depth"))

    @property
    def min_dists(self):
        return np.float64(self._extract_column("min_dist"))

    @property
    def mpole_dists(self):
        return np.float64(self._extract_column("mpole_dist"))
