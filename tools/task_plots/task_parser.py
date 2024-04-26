"""A module that defines a generic task parser class.

This module should not be used directly, but instead should be imported by other
modules which need to parse the thread_info* files.

This can be used by any of the analysis scripts which need to parse the
thread_info*files.

This file is part of SWIFT.

Copyright (C) 2024 Will Roper (w.roper@sussex.ac.uk)
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
import random
import numpy as np

# import hardcoded data
from swift_hardcoded_data import TASKTYPES, SUBTYPES


# Set the random seed for reproducibility
random.seed(42)
np.random.seed(42)


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


def hsv_to_rgb(h, s, v):
    """
    Convert HSV color space to RGB color space.

    Args:
        h (float): Hue (0-360 degrees)
        s (float): Saturation (0-1)
        v (float): Value (0-1)

    Returns:
        tuple: Corresponding RGB values scaled between 0 and 255
    """
    c = v * s
    x = c * (1 - abs((h / 60) % 2 - 1))
    m = v - c

    if h < 60:
        r, g, b = c, x, 0
    elif h < 120:
        r, g, b = x, c, 0
    elif h < 180:
        r, g, b = 0, c, x
    elif h < 240:
        r, g, b = 0, x, c
    elif h < 300:
        r, g, b = x, 0, c
    else:
        r, g, b = c, 0, x

    r = (r + m) * 255
    g = (g + m) * 255
    b = (b + m) * 255

    return int(r), int(g), int(b)


def rgb_to_hex(r, g, b):
    """
    Convert RGB color space to Hex color code.

    Args:
        r (int): Red component (0-255)
        g (int): Green component (0-255)
        b (int): Blue component (0-255)

    Returns:
        str: Corresponding Hex color code
    """
    return f"#{r:02x}{g:02x}{b:02x}"


def generate_distinct_colors(n):
    """
    Generate n distinct hex colors using HSV color space.

    Args:
        n (int): Number of distinct colors to generate.

    Returns:
        list: List of hex color codes.
    """
    colors = []
    for i in range(n):
        h = (360 / n) * i  # Evenly space the hue values
        s = 1  # Set saturation to 100%
        v = 1  # Set value to 100%
        r, g, b = hsv_to_rgb(h, s, v)
        hex_color = rgb_to_hex(r, g, b)
        colors.append(hex_color)
    return colors


def assign_colors(labels):
    """
    Generate a dictionary mapping each label to a distinct hex color.

    Args:
        labels (list): A list of labels for which to generate distinct colors.

    Returns:
        dict: A dictionary with labels as keys and hex color codes as values.
    """
    num_labels = len(labels)
    colors = generate_distinct_colors(num_labels)
    random.shuffle(colors)
    return dict(zip(labels, colors))


class Task:
    """A class to hold the data for a single task.

    This class is used to hold the data for a single task. It is used by the
    TaskParser class to hold the data for each task in the file.

    Attributes:
        rank: The rank of the task.
        thread: The thread of the task.
        type: The type of the task.
        subtype: The subtype of the task.
        task: The task string.
        tic: The start time of the task.
        toc: The end time of the task.
        ci_part_count: The number of particles in the ci cell.
        cj_part_count: The number of particles in the cj cell.
        ci_gpart_count: The number of ghost particles in the ci cell.
        cj_gpart_count: The number of ghost particles in the cj cell.
        ci_type: The type of the ci cell.
        cj_type: The type of the cj cell.
        ci_subtype: The subtype of the ci cell.
        cj_subtype: The subtype of the cj cell.
        ci_depth: The depth of the ci cell.
        cj_depth: The depth of the cj cell.
        min_dist: The minimum distance.
        mpole_dist: The multipole distance.
        dt: The duration of the task.
    """

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
        """
        Initialise the Task object.

        Args:
            rank: The rank of the task.
            thread: The thread of the task.
            type_int: The type of the task.
            subtype_int: The subtype of the task.
            tic: The start time of the task.
            toc: The end time of the task.
            ci_part_count: The number of particles in the ci cell.
            cj_part_count: The number of particles in the cj cell.
            ci_gpart_count: The number of ghost particles in the ci cell.
            cj_gpart_count: The number of ghost particles in the cj cell.
            ci_type: The type of the ci cell.
            cj_type: The type of the cj cell.
            ci_subtype: The subtype of the ci cell.
            cj_subtype: The subtype of the cj cell.
            ci_depth: The depth of the ci cell.
            cj_depth: The depth of the cj cell.
            min_dist: The minimum distance.
            mpole_dist: The multipole distance.
        """
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
        """
        Return a string representation of the task.

        Returns:
            A string representation of the task.
        """
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
    """A class to parse the thread_info-step*.dat files.

    This will ingest a thread_info-step*.dat file and parse it to extract the
    data and convert it to a useful human readable format. The methods then
    provide an interface to access this data.

    Attributes:
        filename: The filename of the file to parse.
        name: The name of the data set.
        ranks: A list of ranks to parse.
        verbose: A flag to determine if the parser should print information.
        delta_t: The time range to parse.
        mintic: The start time to parse.
        filename: The filename of the file to parse.
        name: The name of the data set.
        verbose: A flag to determine if the parser should print information.
        delta_t: The time range to parse.
        mintic: The start time to parse.
        end_t: The end time of the data.
        data: The data from the file.
        full_step: The header of the file.
        mpimode: A flag to determine if the file is an MPI file.
        ranks: A list of ranks to parse.
        _col_look_up: A dictionary to look up the columns in the data.
        cpu_clock: The CPU clock frequency.
        nthread: The number of threads.
        tasks: An array of tasks.
        dt: An array of task durations.
    """

    def __init__(
        self, filename, name, ranks=None, verbose=True, delta_t=0.0, mintic=-1
    ):
        """
        Initialise the TaskParser object.

        Args:
            filename: The filename of the file to parse.
            name: The name of the data set.
            ranks: A list of ranks to parse.
            verbose: A flag to determine if the parser should print information.
            delta_t: The time range to parse.
            mintic: The start time to parse.
        """
        # Define the filename
        self.filename = filename

        # Define the name for this data set
        self.name = name

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
        """Load the data from the file."""
        self.data = np.loadtxt(self.filename)
        self.full_step = self.data[0, :]

    def _define_columns(self):
        """
        Define the columns of the data.

        This is needed since the data is stored in different columns depending
        on whether the file is an MPI file or not.

        This populates a look up dictionary which can be utilised to get data
        by label.
        """
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
        """
        Get a column without needing to know the index.

        Args:
            column: The string defining the column to extract.
        """
        return self.data[:, self._col_look_up[column]]

    def _process_header(self, mintic):
        """Process the header to extract metadata."""
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
        """
        Clean up the data.

        This method will remove the header, remove any zero start and end times
        and convert the tics and tocs to ms.
        """
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
        """Parse the tasks creating Task objects."""
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

    def get_mask(
        self,
        ci_type=None,
        cj_type=None,
        ci_subtype=None,
        cj_subtype=None,
        depth=None,
    ):
        """
        Get a mask for the data based on some filters.

        Args:
            ci_type: The type of the ci cell.
            cj_type: The type of the cj cell.
            ci_subtype: The subtype of the ci cell.
            cj_subtype: The subtype of the cj cell.
            depth: The depth of the cells.
        """
        mask = np.ones(len(self.task_labels), dtype=bool)
        if (ci_type is not None and cj_type is None) or (
            ci_type is None and cj_type is not None
        ):
            cell_type = ci_type if ci_type is not None else cj_type
            mask = np.logical_and(
                mask,
                np.logical_or(
                    self.ci_types == cell_type,
                    self.cj_types == cell_type,
                ),
            )
        if ci_type is not None and cj_type is not None:
            mask = np.logical_and(
                mask,
                np.logical_or(
                    np.logical_and(
                        self.ci_types == ci_type,
                        self.cj_types == cj_type,
                    ),
                    np.logical_and(
                        self.ci_types == cj_type,
                        self.cj_types == ci_type,
                    ),
                ),
            )
        if (ci_subtype is not None and cj_subtype is None) or (
            ci_subtype is None and cj_subtype is not None
        ):
            cell_subtype = ci_subtype if ci_subtype is not None else cj_subtype
            mask = np.logical_and(
                mask,
                np.logical_or(
                    self.ci_subtypes == cell_subtype,
                    self.cj_subtypes == cell_subtype,
                ),
            )
        if ci_subtype is not None and cj_subtype is not None:
            mask = np.logical_and(
                mask,
                np.logical_or(
                    np.logical_and(
                        self.ci_subtypes == ci_subtype,
                        self.cj_subtypes == cj_subtype,
                    ),
                    np.logical_and(
                        self.ci_subtypes == cj_subtype,
                        self.cj_subtypes == ci_subtype,
                    ),
                ),
            )
        if depth is not None:
            mask = np.logical_and(
                mask,
                np.logical_or(
                    self.ci_depths == depth, self.cj_depths == depth
                ),
            )

        return mask

    def _get_tasks_with_mask(self, mask=None):
        """
        Get a subset of the tasks.

        Args:
            mask: The mask to apply to the data.
        """
        if mask is not None:
            tasks = self.tasks[mask]
        else:
            tasks = self.tasks
        unique_tasks = np.unique(self.task_labels[mask])
        unique_count = len(unique_tasks)
        return tasks, unique_tasks, unique_count

    def get_tasks(self):
        """Get all tasks."""
        return self._get_tasks_with_mask()

    def get_tasks_on_rank(self, rank):
        """Get tasks filtered by ranks."""
        return self._get_tasks_with_mask(mask=self.task_ranks == rank)

    def get_tasks_at_depth(self, depth):
        """Get tasks filtered by depth."""
        return self._get_tasks_with_mask(
            mask=np.logical_or(
                self.ci_depths == depth,
                self.cj_depths == depth,
            )
        )

    def get_tasks_tictoc_by_thread(
        self,
        task=None,
        ci_type=None,
        cj_type=None,
        ci_subtype=None,
        cj_subtype=None,
        depth=None,
    ):
        """
        Get the tics and tocs split by thread.

        This will return the labels, tics, and tocs split by thread for a given
        set of tasks or all tasks if no filters are applied.

        This method will also generate a color for each label. These colors are
        evenly spaced in the HSV color space and then randomly assigned to each
        label.

        Args:
            task: The specific task to return.
            ci_type: The type of the ci cell.
            cj_type: The type of the cj cell.
            ci_subtype: The subtype of the ci cell.
            cj_subtype: The subtype of the cj cell.
            depth: The depth of the cells.
        """
        # Get a mask
        mask = self.get_mask(
            ci_type=ci_type,
            cj_type=cj_type,
            ci_subtype=ci_subtype,
            cj_subtype=cj_subtype,
            depth=depth,
        )

        # Combine the task into the mask
        if task is not None:
            mask = np.logical_and(mask, self.task_labels == task)

        # Extract the data
        _labels = self.task_labels[mask]
        _tics = self.tics[mask]
        _tocs = self.tocs[mask]
        threads = self.task_threads[mask]

        # Count the threads
        nthreads = threads.max() + 1

        # Create dictionaries split by thread
        labels = {tid: [] for tid in range(nthreads)}
        tics = {tid: [] for tid in range(nthreads)}
        tocs = {tid: [] for tid in range(nthreads)}

        # Populate the dictionaries
        for i in range(len(labels)):
            labels[threads[i]].append(_labels[i])
            tics[threads[i]].append(_tics[i])
            tocs[threads[i]].append(_tocs[i])

        return labels, tics, tocs, assign_colors(labels)

    @property
    def task_ranks(self):
        """Get the ranks of the tasks."""
        return np.int32(self._extract_column("rank"))

    @property
    def task_threads(self):
        """Get the threads of the tasks."""
        return np.int32(self._extract_column("threads"))

    @property
    def task_types(self):
        """Get the types of the tasks."""
        return np.int32(self._extract_column("task"))

    @property
    def task_subtypes(self):
        """Get the subtypes of the tasks."""
        return np.int32(self._extract_column("subtask"))

    @property
    def tics(self):
        """Get the tics of the tasks."""
        return np.float64(self._extract_column("tic"))

    @property
    def tocs(self):
        """Get the tocs of the tasks."""
        return np.float64(self._extract_column("toc"))

    @property
    def ci_part_count(self):
        """Get the number of hydro particles in the ci cell."""
        return np.int32(self._extract_column("ci_part_count"))

    @property
    def cj_part_count(self):
        """Get the number of hydro particles in the cj cell."""
        return np.int32(self._extract_column("cj_part_count"))

    @property
    def ci_gpart_count(self):
        """Get the number of gravity particles in the ci cell."""
        return np.int32(self._extract_column("ci_gpart_count"))

    @property
    def cj_gpart_count(self):
        """Get the number of gravity particles in the cj cell."""
        return np.int32(self._extract_column("cj_gpart_count"))

    @property
    def ci_types(self):
        """Get the types of the ci cells."""
        return np.int32(self._extract_column("ci_type"))

    @property
    def cj_types(self):
        """Get the types of the cj cells."""
        return np.int32(self._extract_column("cj_type"))

    @property
    def ci_subtypes(self):
        """Get the subtypes of the ci cells."""
        return np.int32(self._extract_column("ci_subtype"))

    @property
    def cj_subtypes(self):
        """Get the subtypes of the cj cells."""
        return np.int32(self._extract_column("cj_subtype"))

    @property
    def ci_depths(self):
        """Get the depths of the ci cells."""
        return np.int32(self._extract_column("ci_depth"))

    @property
    def cj_depths(self):
        """Get the depths of the cj cells."""
        return np.int32(self._extract_column("cj_depth"))

    @property
    def min_dists(self):
        """Get the minimum distances."""
        return np.float64(self._extract_column("min_dist"))

    @property
    def mpole_dists(self):
        """Get the multipole distances."""
        return np.float64(self._extract_column("mpole_dist"))
