#!/usr/bin/env python3
"""
Usage:

python analyse_zoom_tasks.py --files file1 file2 --labels label1 label2

Running this script will produce a series of plots that compare different
thread_info-step*.dat files produced by SWIFT when configured with
--enable-task-debugging. This will produce a series of diagnostic plots that
show the task counts, runtime and distance metrics for different types of tasks
in the simulation.

This script is specifically designed to work with zoom simulations, where the
tasking is complicated by having several different types of cell and complex
geometries.

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
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

from task_parser import TaskParser


def make_task_hist_time_split(runs, order_by_count=True, output=""):
    """
    Plot the runtime of each task type including cell and depth information.

    Args:
        runs: Dictionary of TaskParser objects.
        order_by_count: If True, order the tasks by count, otherwise order
                        by time.
        output: Output filepath.
    """
    fig = plt.figure(figsize=(12, 16))
    ax = fig.add_subplot(111)
    ax.set_xscale("log")
    ax.grid(True)

    # Combine all information into the labels
    labels_dict = {
        name: np.zeros(run.ntasks, dtype=object) for name, run in runs.items()
    }
    time_dict = {}
    for name, run in runs.items():
        for i in range(run.ntasks):
            task = run.task_labels[i]
            labels_dict[name][i] = f"{task}:{run.tasks[i].ci_type}"
            if run.tasks[i].ci_subtype != "Regular":
                labels_dict[name][i] += f"({run.tasks[i].ci_subtype})"
            if "pair" in task:
                labels_dict[name][i] += f"->{run.tasks[i].cj_type}"
                if run.tasks[i].cj_subtype != "Regular":
                    labels_dict[name][i] += f"({run.tasks[i].cj_subtype})"
            labels_dict[name][i] += f"@{run.tasks[i].ci_depth}"
            time_dict[labels_dict[name][i]] = (
                time_dict.get(labels_dict[name][i], 0) + run.tasks[i].dt
            )

    for i, (name, run) in enumerate(runs.items()):
        labels, counts = np.unique(labels_dict[name], return_counts=True)

        if order_by_count:
            # Sort the labels
            if i == 0:
                sinds = np.argsort(-counts)
            labels = labels[sinds]
            counts = counts[sinds]

            # Unpack the times
            times = np.array([time_dict[lab] for lab in labels])
        else:
            # Unpack the times
            times = np.array([time_dict[lab] for lab in labels])

            # Sort the labels by time
            if i == 0:
                sinds = np.argsort(-times)
            labels = labels[sinds]
            times = times[sinds]

        # Calculate positions for horizontal bars
        positions = np.arange(len(labels))

        # Compute the width between labels
        width = 0.8 / (len(runs) + 1)

        # Create horizontal bar plot
        ax.barh(
            positions + (i * width),
            times,
            height=0.75 / len(runs),
            label=name,
            alpha=0.7,
        )

    ax.set_yticks(np.arange(len(labels)) + 0.2)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()

    ax.set_xlabel("Time (ms)")

    # Place the legend at the bottom of the plot
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=3)

    # Define the filename
    filename = "task_time_comp_split"
    if order_by_count:
        filename += "_count_ordered"
    filename += ".png"

    fig.tight_layout()
    fig.savefig(f"{output}{filename}", bbox_inches="tight")

    plt.close(fig)


def make_task_hist_split(runs, output=""):
    """
    Plot the count of each task type including cell and depth information.

    Args:
        runs: Dictionary of TaskParser objects.
        output: Output filepath.
    """
    fig = plt.figure(figsize=(12, 16))
    ax = fig.add_subplot(111)
    ax.set_xscale("log")
    ax.grid(True)

    # Combine all information into the labels
    labels_dict = {
        name: np.zeros(run.ntasks, dtype=object) for name, run in runs.items()
    }
    for name, run in runs.items():
        for i in range(run.ntasks):
            task = run.task_labels[i]
            if (
                run.tasks[i].ci_subtype != "Regular"
                or run.tasks[i].cj_subtype != "Regular"
            ):
                print(
                    "Getting labels:",
                    run.tasks[i].ci_subtype,
                    run.tasks[i].cj_subtype,
                )
            labels_dict[name][i] = f"{task}:{run.tasks[i].ci_type}"
            # if run.tasks[i].ci_subtype != "Regular":
            labels_dict[name][i] += f"({run.tasks[i].ci_subtype})"
            if "pair" in task:
                labels_dict[name][i] += f"->{run.tasks[i].cj_type}"
                # if run.tasks[i].cj_subtype != "Regular":
                labels_dict[name][i] += f"({run.tasks[i].cj_subtype})"
            labels_dict[name][i] += f"@{run.tasks[i].ci_depth}"

    for i, (name, run) in enumerate(runs.items()):
        labels, counts = np.unique(labels_dict[name], return_counts=True)

        # Sort the labels
        if i == 0:
            sinds = np.argsort(-counts)
        labels = labels[sinds]
        counts = counts[sinds]

        for lab in labels:
            print(lab)

        # Calculate positions for horizontal bars
        positions = np.arange(len(labels))

        # Compute the width between labels
        width = 0.8 / (len(runs) + 1)

        # Create horizontal bar plot
        ax.barh(
            positions + (i * width),
            counts,
            height=0.75 / len(runs),
            label=name,
            alpha=0.7,
        )

    ax.set_yticks(np.arange(len(labels)) + 0.2)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()

    ax.set_xlabel("Count")

    # Place the legend at the bottom of the plot
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05), ncol=3)

    # Define the filename
    filename = "task_count_comp_split.png"

    fig.tight_layout()
    fig.savefig(f"{output}{filename}", bbox_inches="tight")

    plt.close(fig)


def make_task_hist(
    runs,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
    output="",
):
    """
    Plot the count of each task type.

    If any of the filter arguments are specifie then a subset of tasks will be
    plotted. If all are None then all tasks will be plotted.

    Note that ci and ci are symmetric, i.e. ci_type=1, cj_type=2 is the same
    as ci_type=2, cj_type=1.

    Args:
        runs: Dictionary of TaskParser objects.
        ci_type: Cell type of the first cell to filter on. (Optional)
        cj_type: Cell type of the second cell to filter on. (Optional)
        ci_subtype: Cell subtype of the first cell to filter on. (Optional)
        cj_subtype: Cell subtype of the second cell to filter on. (Optional)
        depth: Depth of the tasks to filter on. (Optional)
        output: Output filepath.
    """
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.set_xscale("log")
    ax.grid(True)

    # Combine all information into the labels
    nempty_runs = 0
    need_sort = True
    for i, (name, run) in enumerate(runs.items()):
        mask = run.get_mask(ci_type, cj_type, ci_subtype, cj_subtype, depth)

        # Check we have something to plot
        if mask.sum() == 0:
            nempty_runs += 1
            continue

        labels, counts = np.unique(run.task_labels[mask], return_counts=True)

        # Sort the labels and counts by counts in descending order
        if need_sort:
            sorted_indices = np.argsort(-counts)
            need_sort = False
        labels = labels[sorted_indices]
        counts = counts[sorted_indices]

        # Calculate positions for horizontal bars
        positions = np.arange(len(labels))

        # Compute the width between labels
        width = 0.8 / (len(runs) + 1)

        # Create horizontal bar plot
        ax.barh(
            positions + (i * width),
            counts,
            height=0.75 / len(runs),
            label=name,
            alpha=0.7,
        )

    # Exit if there's nothing to plot
    if nempty_runs == len(runs):
        print(
            f"Nothing to plot for ci_type={ci_type}, cj_type={cj_type} "
            f"ci_subtype={ci_subtype}, cj_subtype={cj_subtype}, depth={depth}"
        )
        return

    ax.set_yticks(np.arange(len(labels)) + 0.2)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()

    ax.set_xlabel("Count")

    # Place the legend at the bottom of the plot
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3)

    # Define the filename
    filename = "task_count_comp"
    if ci_type is not None and cj_type is not None:
        filename += f"_types{ci_type}-{cj_type}"
    if ci_subtype is not None and cj_subtype is not None:
        filename += f"_subtypes{ci_subtype}-{cj_subtype}"
    if ci_type is not None and cj_type is None:
        filename += f"_type{ci_type}"
    if ci_subtype is not None and cj_subtype is None:
        filename += f"_subtype{ci_subtype}"
    if ci_type is None and cj_type is not None:
        filename += f"_type{cj_type}"
    if ci_subtype is None and cj_subtype is not None:
        filename += f"_subtype{cj_subtype}"
    if depth is not None:
        filename += f"_depth{depth}"

    fig.tight_layout()
    fig.savefig(f"{output}{filename}", bbox_inches="tight")

    plt.close(fig)


def make_task_hist_time_weighted(
    runs,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
    output="",
):
    """
    Plot the runtime of each task type.

    If any of the filter arguments are specifie then a subset of tasks will be
    plotted. If all are None then all tasks will be plotted.

    Note that ci and ci are symmetric, i.e. ci_type=1, cj_type=2 is the same
    as ci_type=2, cj_type=1.

    Args:
        runs: Dictionary of TaskParser objects.
        ci_type: Cell type of the first cell to filter on. (Optional)
        cj_type: Cell type of the second cell to filter on. (Optional)
        ci_subtype: Cell subtype of the first cell to filter on. (Optional)
        cj_subtype: Cell subtype of the second cell to filter on. (Optional)
        depth: Depth of the tasks to filter on. (Optional)
        output: Output filepath.
    """
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.set_xscale("log")
    ax.grid(True)

    # Combine all information into the labels
    nempty_runs = 0
    need_sort = True
    for i, (name, run) in enumerate(runs.items()):
        mask = run.get_mask(ci_type, cj_type, ci_subtype, cj_subtype, depth)

        # Check we have something to plot
        if mask.sum() == 0:
            nempty_runs += 1
            continue

        # Loop over tasks collecting their runtime
        labels = np.unique(run.task_labels[mask])
        counts = np.array(
            [np.sum(run.dt[mask][run.task_labels[mask] == k]) for k in labels]
        )

        # Sort the labels and counts by counts in descending order
        if need_sort:
            sorted_indices = np.argsort(-counts)
            need_sort = False
        labels = labels[sorted_indices]
        counts = counts[sorted_indices]

        # Calculate positions for horizontal bars
        positions = np.arange(len(labels))

        # Compute the width between labels
        width = 0.8 / (len(runs) + 1)

        # Create horizontal bar plot
        ax.barh(
            positions + (i * width),
            counts,
            height=0.75 / len(runs),
            label=name,
            alpha=0.7,
        )

    # Exit if there's nothing to plot
    if nempty_runs == len(runs):
        print(
            f"Nothing to plot for ci_type={ci_type}, cj_type={cj_type} "
            f"ci_subtype={ci_subtype}, cj_subtype={cj_subtype}, depth={depth}"
        )
        return

    ax.set_yticks(np.arange(len(labels)) + 0.2)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()

    ax.set_xlabel("Time (ms)")

    # Place the legend at the bottom of the plot
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3)

    # Define the filename
    filename = "task_time_comp"
    if ci_type is not None and cj_type is not None:
        filename += f"_types{ci_type}-{cj_type}"
    if ci_subtype is not None and cj_subtype is not None:
        filename += f"_subtypes{ci_subtype}-{cj_subtype}"
    if ci_type is not None and cj_type is None:
        filename += f"_type{ci_type}"
    if ci_subtype is not None and cj_subtype is None:
        filename += f"_subtype{ci_subtype}"
    if ci_type is None and cj_type is not None:
        filename += f"_type{cj_type}"
    if ci_subtype is None and cj_subtype is not None:
        filename += f"_subtype{cj_subtype}"
    if depth is not None:
        filename += f"_depth{depth}"

    fig.tight_layout()
    fig.savefig(f"{output}{filename}", bbox_inches="tight")

    plt.close(fig)


def make_pair_mindist_plot(
    runs,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
    nbins=30,
    output="",
):
    """
    Histogram of the minimum distances between cells.

    This will histogram the output of sqrt(cell_min_dist2).

    If any of the filter arguments are specifie then a subset of tasks will be
    plotted. If all are None then all tasks will be plotted.

    Note that ci and ci are symmetric, i.e. ci_type=1, cj_type=2 is the same
    as ci_type=2, cj_type=1.

    Args:
        runs: Dictionary of TaskParser objects.
        ci_type: Cell type of the first cell to filter on. (Optional)
        cj_type: Cell type of the second cell to filter on. (Optional)
        ci_subtype: Cell subtype of the first cell to filter on. (Optional)
        cj_subtype: Cell subtype of the second cell to filter on. (Optional)
        depth: Depth of the tasks to filter on. (Optional)
        nbins: Number of bins in the histogram.
        output: Output filepath.
    """
    # Make the figure
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.set_yscale("log")
    ax.grid(True)

    # Collect the distances
    dists = {}
    for i, (name, run) in enumerate(runs.items()):
        mask = run.get_mask(ci_type, cj_type, ci_subtype, cj_subtype, depth)

        # Ensure we only have pair tasks (i.e. the string "pair" is in the
        # task label)
        mask = np.logical_and(
            mask, np.array(["pair" in t for t in run.task_labels])
        )

        # Get the distances
        dists[name] = run.min_dists[mask]

    # Collect together all the distances
    all_dists = np.concatenate(list(dists.values()))

    # Exit if there are no distances
    if all_dists.size == 0:
        print(
            f"No distances to plot for ci_type={ci_type}, cj_type={cj_type} "
            f"ci_subtype={ci_subtype}, cj_subtype={cj_subtype}, depth={depth}"
        )
        return

    # Construct the bins
    bins = np.linspace(all_dists.min(), all_dists.max(), nbins + 1)
    bin_cents = (bins[:-1] + bins[1:]) / 2

    # Compute histogram and plot
    for name in dists.keys():
        linestyle = "--" if "long_range" in name else "-"
        H, _ = np.histogram(dists[name], bins=bins)
        ax.plot(bin_cents, H, label=name, linestyle=linestyle)

    ax.set_xlabel("sqrt(cell_min_dist2) (U_L)")
    ax.set_ylabel("Count")

    # Place the legend at the bottom of the plot
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3)

    # Define the filename
    filename = "pair_min_dist_comp"
    if ci_type is not None and cj_type is not None:
        filename += f"_types{ci_type}-{cj_type}"
    if ci_subtype is not None and cj_subtype is not None:
        filename += f"_subtypes{ci_subtype}-{cj_subtype}"
    if ci_type is not None and cj_type is None:
        filename += f"_type{ci_type}"
    if ci_subtype is not None and cj_subtype is None:
        filename += f"_subtype{ci_subtype}"
    if ci_type is None and cj_type is not None:
        filename += f"_type{cj_type}"
    if ci_subtype is None and cj_subtype is not None:
        filename += f"_subtype{cj_subtype}"
    if depth is not None:
        filename += f"_depth{depth}"

    fig.tight_layout()
    fig.savefig(f"{output}{filename}", bbox_inches="tight")

    plt.close(fig)


def make_pair_mpoledist_plot(
    runs,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
    nbins=30,
    output="",
):
    """
    Histogram of the minimum distances between multipoles.

    This will histogram the distances between the centres of mass of the
    multipoles which have pair tasks.

    If any of the filter arguments are specifie then a subset of tasks will be
    plotted. If all are None then all tasks will be plotted.

    Note that ci and ci are symmetric, i.e. ci_type=1, cj_type=2 is the same
    as ci_type=2, cj_type=1.

    Args:
        runs: Dictionary of TaskParser objects.
        ci_type: Cell type of the first cell to filter on. (Optional)
        cj_type: Cell type of the second cell to filter on. (Optional)
        ci_subtype: Cell subtype of the first cell to filter on. (Optional)
        cj_subtype: Cell subtype of the second cell to filter on. (Optional)
        depth: Depth of the tasks to filter on. (Optional)
        nbins: Number of bins in the histogram.
        output: Output filepath.
    """
    # Make the figure
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.set_yscale("log")
    ax.grid(True)

    # Collect the distances
    dists = {}
    for i, (name, run) in enumerate(runs.items()):
        mask = run.get_mask(ci_type, cj_type, ci_subtype, cj_subtype, depth)

        # Ensure we only have pair tasks (i.e. the string "pair" is in the
        # task label)
        mask = np.logical_and(
            mask, np.array(["pair" in t for t in run.task_labels])
        )

        # Get the distances
        dists[name] = run.mpole_dists[mask]

    # Collect together all the distances
    all_dists = np.concatenate(list(dists.values()))

    # Exit if there are no distances
    if all_dists.size == 0:
        print(
            f"No distances to plot for ci_type={ci_type}, cj_type={cj_type} "
            f"ci_subtype={ci_subtype}, cj_subtype={cj_subtype}, depth={depth}"
        )
        return

    # Construct the bins
    bins = np.linspace(all_dists.min(), all_dists.max(), nbins + 1)
    bin_cents = (bins[:-1] + bins[1:]) / 2

    # Compute histogram and plot
    for name in dists.keys():
        linestyle = "--" if "long_range" in name else "-"
        H, _ = np.histogram(dists[name], bins=bins)
        ax.plot(bin_cents, H, label=name, linestyle=linestyle)

    ax.set_xlabel("Multipole CoM distance (U_L)")
    ax.set_ylabel("Count")

    # Place the legend at the bottom of the plot
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3)

    # Define the filename
    filename = "pair_mpole_dist_comp"
    if ci_type is not None and cj_type is not None:
        filename += f"_types{ci_type}-{cj_type}"
    if ci_subtype is not None and cj_subtype is not None:
        filename += f"_subtypes{ci_subtype}-{cj_subtype}"
    if ci_type is not None and cj_type is None:
        filename += f"_type{ci_type}"
    if ci_subtype is not None and cj_subtype is None:
        filename += f"_subtype{ci_subtype}"
    if ci_type is None and cj_type is not None:
        filename += f"_type{cj_type}"
    if ci_subtype is None and cj_subtype is not None:
        filename += f"_subtype{cj_subtype}"
    if depth is not None:
        filename += f"_depth{depth}"

    fig.tight_layout()
    fig.savefig(f"{output}{filename}", bbox_inches="tight")
    plt.close(fig)


def make_mindist_mpoledist_comp(
    run,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
    output="",
):
    """
    Plot a scatter comparing the different distances.

    This function will compare the minimum distance between cells and the
    distance between the multipoles.

    Unlike other functions this function acts on a single run at a time.

    Args:
        run: A TaskParser object.
        ci_type: Cell type of the first cell to filter on. (Optional)
        cj_type: Cell type of the second cell to filter on. (Optional)
        ci_subtype: Cell subtype of the first cell to filter on. (Optional)
        cj_subtype: Cell subtype of the second cell to filter on. (Optional)
        depth: Depth of the tasks to filter on. (Optional)
        output: Output filepath.
    """
    # Make the mask
    mask = run.get_mask(ci_type, cj_type, ci_subtype, cj_subtype, depth)

    # Ensure we only have pair tasks (i.e. the string "pair" is in the
    # task label)
    mask = np.logical_and(
        mask, np.array(["pair" in t for t in run.task_labels])
    )

    # Exit if there are no distances
    if np.sum(mask) == 0:
        print(
            f"No distances to plot for ci_type={ci_type}, cj_type={cj_type} "
            f"ci_subtype={ci_subtype}, cj_subtype={cj_subtype}, depth={depth}"
        )
        return

    # Get the distances
    min_dists = run.min_dists[mask]
    mpole_dists = run.mpole_dists[mask]

    # Set up the figure
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.loglog()
    ax.grid(True)
    ax.scatter(min_dists, mpole_dists, marker=".", color="grey", alpha=0.7)
    ax.set_xlabel("Minimum distance between cells (U_L)")
    ax.set_ylabel("Distance between multipoles (U_L)")

    # Define the filename
    filename = "min_dist_mpole_dist_comp"
    if ci_type is not None and cj_type is not None:
        filename += f"_types{ci_type}-{cj_type}"
    if ci_subtype is not None and cj_subtype is not None:
        filename += f"_subtypes{ci_subtype}-{cj_subtype}"
    if ci_type is not None and cj_type is None:
        filename += f"_type{ci_type}"
    if ci_subtype is not None and cj_subtype is None:
        filename += f"_subtype{ci_subtype}"
    if ci_type is None and cj_type is not None:
        filename += f"_type{cj_type}"
    if ci_subtype is None and cj_subtype is not None:
        filename += f"_subtype{cj_subtype}"
    if depth is not None:
        filename += f"_depth{depth}"

    fig.tight_layout()

    # Save the figure
    fig.savefig(f"{output}{run.name}_{filename}.png", bbox_inches="tight")
    plt.close(fig)


def make_task_plot(
    run,
    task_type=None,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
    output="",
):
    """
    Plot tasks as a function of time per thread.

    This will give the same output as plot_tasks.py but offers more flexibility
    in the filters that can be applied.

    Args:
        run: A TaskParser object.
        task_type: The type of task to plot to filter on. (Optional)
        ci_type: Cell type of the first cell to filter on. (Optional)
        cj_type: Cell type of the second cell to filter on. (Optional)
        ci_subtype: Cell subtype of the first cell to filter on. (Optional)
        cj_subtype: Cell subtype of the second cell to filter on. (Optional)
        depth: Depth of the tasks to filter on. (Optional)
    """
    # If we have nothing then exit and move on
    mask = run.get_mask(ci_type, cj_type, ci_subtype, cj_subtype, depth)
    ntasks_tot = mask.sum()
    if ntasks_tot == 0:
        print(
            f"No tasks to plot for task_type={task_type} ci_type={ci_type} "
            f"cj_type={cj_type} ci_subtype={ci_subtype} "
            f"cj_subtype={cj_subtype} depth={depth}"
        )
        return

    # Get the dictionaries containing the labels, tics, tocs and colors
    labels, tics, tocs, colors = run.get_tasks_tictoc_by_thread(
        task_type, ci_type, cj_type, ci_subtype, cj_subtype, depth
    )

    # Get the unique tasks
    unique_tasks, task_counts = np.unique(
        run.task_labels[mask], return_counts=True
    )
    ntasks = len(unique_tasks)

    # Set up the figure
    fig = plt.figure(figsize=(16, 0.15 * run.nthread))
    ax = fig.add_subplot(111)
    ax.grid(True, zorder=0)
    ax.set_xlim(-run.delta_t * 0.01, run.delta_t * 1.01)
    ax.set_ylim(0.5, run.nthread + 1.0)

    # Loop over threads plotting the tasks
    for i in labels.keys():
        # Collect the ranges and colors into lists
        _tictocs = []
        _colors = []

        # Loop over the tasks
        for j in range(len(labels[i])):
            _tictocs.append((tics[i][j], tocs[i][j] - tics[i][j]))
            _colors.append(colors[labels[i][j]])

        ax.broken_barh(
            _tictocs,
            [i + 0.55, 0.9],
            facecolors=_colors,
            linewidth=0,
            zorder=1,
        )

    # Create legend handles from sorted labels and colors
    sinds = np.argsort(-task_counts)
    sorted_labels = unique_tasks[sinds]
    sorted_colors = [colors[lab] for lab in sorted_labels]
    handles = [
        mpatches.Patch(color=color, label=label)
        for label, color in zip(sorted_labels, sorted_colors)
    ]

    # Legend and room for it.
    nrow = ntasks / 8
    ax.fill_between([0, 0], run.nthread, run.nthread + nrow, facecolor="white")
    ax.legend(
        handles=handles,
        loc="lower left",
        shadow=True,
        bbox_to_anchor=(0.0, 1.0, 1.0, 0.2),
        mode="expand",
        ncol=8,
    )

    # Start and end of time-step
    ax.plot([0, 0], [0, run.nthread + nrow + 1], "k--", linewidth=1)
    ax.plot(
        [run.delta_t, run.delta_t],
        [0, run.nthread + nrow + 1],
        "k--",
        linewidth=1,
        zorder=1,
    )

    ax.set_xlabel("Wall clock time [ms]")
    ax.set_ylabel("Thread ID")

    # Define the filename
    filename = "tasks"
    if ci_type is not None and cj_type is not None:
        filename += f"_types{ci_type}-{cj_type}"
    if ci_subtype is not None and cj_subtype is not None:
        filename += f"_subtypes{ci_subtype}-{cj_subtype}"
    if ci_type is not None and cj_type is None:
        filename += f"_type{ci_type}"
    if ci_subtype is not None and cj_subtype is None:
        filename += f"_subtype{ci_subtype}"
    if ci_type is None and cj_type is not None:
        filename += f"_type{cj_type}"
    if ci_subtype is None and cj_subtype is not None:
        filename += f"_subtype{cj_subtype}"
    if depth is not None:
        filename += f"_depth{depth}"
    if run.mpimode:
        filename += f"_rank{run.rank}"

    fig.tight_layout()

    # Save the figure
    fig.savefig(f"{output}{run.name}_{filename}.png", bbox_inches="tight")
    plt.close(fig)


def make_task_activity_plot(
    run,
    task_type=None,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
    xbins=1000,
    sort_threads=True,
    output="",
):
    """
    Histogram tasks as a function of time per thread.

    This will bin tasks into a grid, counting how many tasks ran in a
    particular bin. This is particularly useful for visualising runs containing
    numerous short lived tasks which will not necessarily appear in a normal
    task plot.

    Args:
        run: A TaskParser object.
        task_type: The type of task to plot to filter on. (Optional)
        ci_type: Cell type of the first cell to filter on. (Optional)
        cj_type: Cell type of the second cell to filter on. (Optional)
        ci_subtype: Cell subtype of the first cell to filter on. (Optional)
        cj_subtype: Cell subtype of the second cell to filter on. (Optional)
        depth: Depth of the tasks to filter on. (Optional)
        xbins: Number of bins in the x-axis. This defines the minimum time
               period that can be resolved.
        sort_threads: Sort the threads by the end of their tasks. (Optional)
    """
    # If we have nothing then exit and move on
    ntasks_tot = run.get_mask(
        ci_type, cj_type, ci_subtype, cj_subtype, depth
    ).sum()
    if ntasks_tot == 0:
        print(
            f"No tasks to plot for task_type={task_type} ci_type={ci_type} "
            f"cj_type={cj_type} ci_subtype={ci_subtype} "
            f"cj_subtype={cj_subtype} depth={depth}"
        )
        return

    # Get the dictionaries containing the labels, tics, tocs and colors
    labels, tics, tocs, colors = run.get_tasks_tictoc_by_thread(
        task_type, ci_type, cj_type, ci_subtype, cj_subtype, depth
    )

    # Exit if there are no tasks to plot
    if len(labels) == 0:
        print(
            f"No tasks to plot for task_type={task_type} ci_type={ci_type} "
            f"cj_type={cj_type} ci_subtype={ci_subtype} "
            f"cj_subtype={cj_subtype} depth={depth}"
        )
        return

    # Define the grid of bins. This has shape (run.nthread, xbins)
    grid = np.zeros((run.nthread, xbins))

    # Populate the grid thread by thread and task by task
    for i in labels.keys():
        for j in range(len(labels[i])):
            # Calculate the tic bin
            xtic_bin = int(tics[i][j] / run.delta_t * xbins)

            # Calculate the toc bin
            xtoc_bin = int(tocs[i][j] / run.delta_t * xbins)

            # Populate the bins
            for xbin in range(xtic_bin, xtoc_bin):
                grid[i, xbin] = 1

    # Sort the rows of the grid
    if sort_threads:
        counts = [np.sum(grid[i, :]) for i in range(run.nthread)]
        sinds = np.argsort(counts)[::-1]
        grid = grid[sinds, :]

    # Set up the figure with a main plot for the grid and a histogram on top
    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(
        2, 2, height_ratios=[2, 10], width_ratios=[20, 1], hspace=0.0
    )
    ax_grid = fig.add_subplot(gs[1, 0])
    ax_hist = fig.add_subplot(gs[0, 0])
    ax_hist.grid(True)
    ax_hist.set_axisbelow(True)
    cax = fig.add_subplot(gs[:, 1])

    # Plot the grid
    im = ax_grid.pcolormesh(
        np.linspace(0, run.delta_t, xbins),
        np.linspace(0, run.nthread - 1, run.nthread),
        grid,
        cmap="coolwarm",
    )

    # Start and end of time-step
    ax_grid.plot([0, 0], ax_grid.get_ylim(), "k--", linewidth=1)
    ax_grid.plot(
        [run.delta_t, run.delta_t], ax_grid.get_ylim(), "k--", linewidth=1
    )

    ax_grid.set_xlabel("Wall clock time [ms]")
    ax_grid.set_ylabel("Thread ID")

    # Create the colorbar
    cbar = fig.colorbar(im, cax=cax, ticks=[0, 1])
    cbar.set_ticklabels(["Inactive", "Task"])

    # Ensure the histogram aligns perfectly with the grid plot
    ax_hist.set_xlim(ax_grid.get_xlim())

    # Calculate the sum across each column for the histogram
    column_sums = grid.sum(axis=0)
    ax_hist.bar(
        np.linspace(0, run.delta_t, xbins),
        column_sums,
        width=np.diff(np.linspace(0, run.delta_t, xbins))[0],
        align="edge",
        zorder=1,
    )

    ax_hist.set_ylabel("Tasks Running")
    ax_hist.set_ylim(0, run.nthread + 0.1 * run.nthread)

    # Turn off the x-axis ticks and tick labels on the histogram
    ax_hist.tick_params(axis="x", which="both", bottom=False, top=False)
    ax_hist.set_xticklabels([])

    # Define the filename
    filename = "task_activity"
    if ci_type is not None and cj_type is not None:
        filename += f"_types{ci_type}-{cj_type}"
    if ci_subtype is not None and cj_subtype is not None:
        filename += f"_subtypes{ci_subtype}-{cj_subtype}"
    if ci_type is not None and cj_type is None:
        filename += f"_type{ci_type}"
    if ci_subtype is not None and cj_subtype is None:
        filename += f"_subtype{ci_subtype}"
    if ci_type is None and cj_type is not None:
        filename += f"_type{cj_type}"
    if ci_subtype is None and cj_subtype is not None:
        filename += f"_subtype{cj_subtype}"
    if depth is not None:
        filename += f"_depth{depth}"
    if run.mpimode:
        filename += f"_rank{run.rank}"
    if sort_threads:
        filename += "_sorted"

    fig.tight_layout()

    # Save the figure
    fig.savefig(f"{output}{run.name}_{filename}.png", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    # Define the command line arguments
    parser = argparse.ArgumentParser(
        description="Produce task analysis plots for SWIFT zoom simulations"
    )

    # Adding files argument
    parser.add_argument(
        "--files",
        nargs="+",
        help="List of files to combine on outputs",
        required=True,
    )

    # Adding labels argument
    parser.add_argument(
        "--labels", nargs="+", help="List of labels", default=[]
    )

    # Adding output directory
    parser.add_argument("--outdir", help="Output directory", default=".")

    # Adding output base name
    parser.add_argument("--outbase", help="Output base name", default="")

    # Parse the arguments
    args = parser.parse_args()
    files = args.files
    labels = args.labels
    outdir = args.outdir
    outbase = args.outbase
    output = f"{outdir}/{outbase}"
    if output[-1] != "_" and output[-1] != "/":
        output += "_"

    print("Writing outputs to:", output)

    if len(labels) == 0:
        labels = files
        print("Using filenames as labels")

    if len(files) != len(labels):
        raise ValueError("Number of files and labels must match")

    # Parse all the task files
    runs = {}
    for f, lab in zip(files, labels):
        runs[lab] = TaskParser(f, lab)

    # Below we call the functions for all common useful combinations of filters

    # Detailed stacked histograms
    make_task_hist_split(runs, output=output)
    make_task_hist_time_split(runs, output=output)
    make_task_hist_time_split(runs, order_by_count=False, output=output)

    # Counts of tasks
    make_task_hist(runs, output=output)
    make_task_hist(runs, ci_type=1, cj_type=1, output=output)
    make_task_hist(runs, ci_type=2, cj_type=2, output=output)
    make_task_hist(runs, ci_type=3, cj_type=3, output=output)
    make_task_hist(runs, ci_type=1, cj_type=2, output=output)
    make_task_hist(runs, ci_type=1, cj_type=3, output=output)
    make_task_hist(runs, ci_type=2, cj_type=3, output=output)

    # Count useful subtypes tasks to be sure there are none
    make_task_hist(runs, ci_subtype=0, output=output)  # regular
    make_task_hist(runs, ci_subtype=1, output=output)  # neighbour
    make_task_hist(runs, ci_subtype=2, output=output)  # void

    # Counts of tasks but only depth 0
    make_task_hist(runs, depth=0, output=output)
    make_task_hist(runs, ci_type=1, cj_type=1, depth=0, output=output)
    make_task_hist(runs, ci_type=2, cj_type=2, depth=0, output=output)
    make_task_hist(runs, ci_type=3, cj_type=3, depth=0, output=output)
    make_task_hist(runs, ci_type=1, cj_type=3, depth=0, output=output)
    make_task_hist(runs, ci_type=1, cj_type=2, depth=0, output=output)
    make_task_hist(runs, ci_type=2, cj_type=3, depth=0, output=output)

    # Time weighted counts of tasks
    make_task_hist_time_weighted(runs, output=output)
    make_task_hist_time_weighted(runs, ci_type=1, cj_type=1, output=output)
    make_task_hist_time_weighted(runs, ci_type=2, cj_type=2, output=output)
    make_task_hist_time_weighted(runs, ci_type=3, cj_type=3, output=output)
    make_task_hist_time_weighted(runs, ci_type=1, cj_type=3, output=output)
    make_task_hist_time_weighted(runs, ci_type=1, cj_type=2, output=output)
    make_task_hist_time_weighted(runs, ci_type=2, cj_type=3, output=output)

    # Time weighted counts of tasks but only depth 0
    make_task_hist_time_weighted(runs, depth=0, output=output)
    make_task_hist_time_weighted(
        runs, ci_type=1, cj_type=1, depth=0, output=output
    )
    make_task_hist_time_weighted(
        runs, ci_type=2, cj_type=2, depth=0, output=output
    )
    make_task_hist_time_weighted(
        runs, ci_type=3, cj_type=3, depth=0, output=output
    )
    make_task_hist_time_weighted(
        runs, ci_type=1, cj_type=3, depth=0, output=output
    )
    make_task_hist_time_weighted(
        runs, ci_type=1, cj_type=2, depth=0, output=output
    )
    make_task_hist_time_weighted(
        runs, ci_type=2, cj_type=3, depth=0, output=output
    )

    # Pair distance histograms
    make_pair_mindist_plot(runs, output=output)
    make_pair_mindist_plot(runs, ci_type=1, cj_type=1, output=output)
    make_pair_mindist_plot(runs, ci_type=2, cj_type=2, output=output)
    make_pair_mindist_plot(runs, ci_type=3, cj_type=3, output=output)
    make_pair_mindist_plot(runs, ci_type=1, cj_type=3, output=output)
    make_pair_mindist_plot(runs, ci_type=1, cj_type=2, output=output)
    make_pair_mindist_plot(runs, ci_type=2, cj_type=3, output=output)

    # Pair multipole distance histograms
    make_pair_mpoledist_plot(runs, output=output)
    make_pair_mpoledist_plot(runs, ci_type=1, cj_type=1, output=output)
    make_pair_mpoledist_plot(runs, ci_type=2, cj_type=2, output=output)
    make_pair_mpoledist_plot(runs, ci_type=3, cj_type=3, output=output)
    make_pair_mpoledist_plot(runs, ci_type=1, cj_type=3, output=output)
    make_pair_mpoledist_plot(runs, ci_type=1, cj_type=2, output=output)
    make_pair_mpoledist_plot(runs, ci_type=2, cj_type=3, output=output)

    # Make the plots which only plot a single run at a time
    for run in runs.values():
        # Distance comparison
        make_mindist_mpoledist_comp(run, output=output)
        make_mindist_mpoledist_comp(run, ci_type=1, cj_type=1, output=output)
        make_mindist_mpoledist_comp(run, ci_type=2, cj_type=2, output=output)
        make_mindist_mpoledist_comp(run, ci_type=3, cj_type=3, output=output)
        make_mindist_mpoledist_comp(run, ci_type=1, cj_type=3, output=output)
        make_mindist_mpoledist_comp(run, ci_type=1, cj_type=2, output=output)
        make_mindist_mpoledist_comp(run, ci_type=2, cj_type=3, output=output)

        # Make the task plots showing tasks per thread as a function of time
        make_task_plot(run, output=output)
        make_task_plot(run, ci_type=1, output=output)
        make_task_plot(run, ci_type=2, output=output)
        make_task_plot(run, ci_type=3, output=output)
        make_task_plot(run, ci_type=1, cj_type=1, output=output)
        make_task_plot(run, ci_type=2, cj_type=2, output=output)
        make_task_plot(run, ci_type=3, cj_type=3, output=output)
        make_task_plot(run, ci_type=1, cj_type=3, output=output)
        make_task_plot(run, ci_type=1, cj_type=2, output=output)
        make_task_plot(run, ci_type=2, cj_type=3, output=output)

        # Make task activity plots
        make_task_activity_plot(run, output=output)
        make_task_activity_plot(run, ci_type=1, output=output)
        make_task_activity_plot(run, ci_type=2, output=output)
        make_task_activity_plot(run, ci_type=3, output=output)
        make_task_activity_plot(run, ci_type=1, cj_type=1, output=output)
        make_task_activity_plot(run, ci_type=2, cj_type=2, output=output)
        make_task_activity_plot(run, ci_type=3, cj_type=3, output=output)
        make_task_activity_plot(run, ci_type=1, cj_type=3, output=output)
        make_task_activity_plot(run, ci_type=1, cj_type=2, output=output)
        make_task_activity_plot(run, ci_type=2, cj_type=3, output=output)

        # Make task activity plots but don't sort the threads
        make_task_activity_plot(run, sort_threads=False, output=output)
        make_task_activity_plot(
            run, ci_type=1, sort_threads=False, output=output
        )
        make_task_activity_plot(
            run, ci_type=2, sort_threads=False, output=output
        )
        make_task_activity_plot(
            run, ci_type=3, sort_threads=False, output=output
        )
        make_task_activity_plot(
            run, ci_type=1, cj_type=1, sort_threads=False, output=output
        )
        make_task_activity_plot(
            run, ci_type=2, cj_type=2, sort_threads=False, output=output
        )
        make_task_activity_plot(
            run, ci_type=3, cj_type=3, sort_threads=False, output=output
        )
        make_task_activity_plot(
            run, ci_type=1, cj_type=3, sort_threads=False, output=output
        )
        make_task_activity_plot(
            run, ci_type=1, cj_type=2, sort_threads=False, output=output
        )
        make_task_activity_plot(
            run, ci_type=2, cj_type=3, sort_threads=False, output=output
        )
