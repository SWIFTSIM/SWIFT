#!/usr/bin/env python3
"""
Usage:


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

from task_parser import TaskParser


def make_mask(
    run,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
):
    mask = np.ones(len(run.task_labels), dtype=bool)
    if (ci_type is not None and cj_type is None) or (
        ci_type is None and cj_type is not None
    ):
        cell_type = ci_type if ci_type is not None else cj_type
        mask = np.logical_and(
            mask,
            np.logical_or(
                run.ci_types == cell_type,
                run.cj_types == cell_type,
            ),
        )
    if ci_type is not None and cj_type is not None:
        mask = np.logical_and(
            mask,
            np.logical_or(
                np.logical_and(
                    run.ci_types == ci_type,
                    run.cj_types == cj_type,
                ),
                np.logical_and(
                    run.ci_types == cj_type,
                    run.cj_types == ci_type,
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
                run.ci_subtypes == cell_subtype,
                run.cj_subtypes == cell_subtype,
            ),
        )
    if ci_subtype is not None and cj_subtype is not None:
        mask = np.logical_and(
            mask,
            np.logical_or(
                np.logical_and(
                    run.ci_subtypes == ci_subtype,
                    run.cj_subtypes == cj_subtype,
                ),
                np.logical_and(
                    run.ci_subtypes == cj_subtype,
                    run.cj_subtypes == ci_subtype,
                ),
            ),
        )
    if depth is not None:
        mask = np.logical_and(
            mask,
            np.logical_or(run.ci_depths == depth, run.cj_depths == depth),
        )

    return mask


def make_task_hist_time_split(runs, order_by_count=True):
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
            times = np.array([time_dict[l] for l in labels])
        else:
            # Unpack the times
            times = np.array([time_dict[l] for l in labels])

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
        bars = ax.barh(
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
    fig.savefig(filename, bbox_inches="tight")

    plt.close(fig)


def make_task_hist_split(runs):
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
            labels_dict[name][i] = f"{task}:{run.tasks[i].ci_type}"
            if run.tasks[i].ci_subtype != "Regular":
                labels_dict[name][i] += f"({run.tasks[i].ci_subtype})"
            if "pair" in task:
                labels_dict[name][i] += f"->{run.tasks[i].cj_type}"
                if run.tasks[i].cj_subtype != "Regular":
                    labels_dict[name][i] += f"({run.tasks[i].cj_subtype})"
            labels_dict[name][i] += f"@{run.tasks[i].ci_depth}"

    for i, (name, run) in enumerate(runs.items()):
        labels, counts = np.unique(labels_dict[name], return_counts=True)

        # Sort the labels
        if i == 0:
            sinds = np.argsort(-counts)
        labels = labels[sinds]
        counts = counts[sinds]

        # Calculate positions for horizontal bars
        positions = np.arange(len(labels))

        # Compute the width between labels
        width = 0.8 / (len(runs) + 1)

        # Create horizontal bar plot
        bars = ax.barh(
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
    fig.savefig(filename, bbox_inches="tight")

    plt.close(fig)


def make_task_hist(
    runs,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
):
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.set_xscale("log")
    ax.grid(True)

    for i, (name, run) in enumerate(runs.items()):
        mask = make_mask(run, ci_type, cj_type, ci_subtype, cj_subtype, depth)

        labels, counts = np.unique(run.task_labels[mask], return_counts=True)

        # Sort the labels and counts by counts in descending order
        if i == 0:
            sorted_indices = np.argsort(-counts)
        labels = labels[sorted_indices]
        counts = counts[sorted_indices]

        # Calculate positions for horizontal bars
        positions = np.arange(len(labels))

        # Compute the width between labels
        width = 0.8 / (len(runs) + 1)

        # Create horizontal bar plot
        bars = ax.barh(
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
    fig.savefig(filename, bbox_inches="tight")

    plt.close(fig)


def make_task_hist_time_weighted(
    runs,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
):
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.set_xscale("log")
    ax.grid(True)

    for i, (name, run) in enumerate(runs.items()):
        mask = make_mask(run, ci_type, cj_type, ci_subtype, cj_subtype, depth)

        # Loop over tasks collecting their runtime
        labels = np.unique(run.task_labels[mask])
        counts = np.array(
            [np.sum(run.dt[mask][run.task_labels[mask] == k]) for k in labels]
        )

        # Sort the labels and counts by counts in descending order
        if i == 0:
            sorted_indices = np.argsort(-counts)
        labels = labels[sorted_indices]
        counts = counts[sorted_indices]

        # Calculate positions for horizontal bars
        positions = np.arange(len(labels))

        # Compute the width between labels
        width = 0.8 / (len(runs) + 1)

        # Create horizontal bar plot
        bars = ax.barh(
            positions + (i * width),
            counts,
            height=0.75 / len(runs),
            label=name,
            alpha=0.7,
        )

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
    fig.savefig(filename, bbox_inches="tight")

    plt.close(fig)


def make_pair_mindist_plot(
    runs,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
    nbins=30,
):
    # Make the figure
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.set_yscale("log")
    ax.grid(True)

    # Collect the distances
    dists = {}
    for i, (name, run) in enumerate(runs.items()):
        mask = make_mask(run, ci_type, cj_type, ci_subtype, cj_subtype, depth)

        # Ensure we only have pair tasks (i.e. the string "pair" is in the
        # task label)
        mask = np.logical_and(
            mask, np.array(["pair" in t for t in run.task_labels])
        )

        # Get the distances
        dists[name] = run.min_dists[mask]

    # Combine all distances
    all_dists = np.concatenate(list(dists.values()))

    # Anything to do here?
    if all_dists.size == 0:
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
    fig.savefig(filename, bbox_inches="tight")

    plt.close(fig)


def make_pair_mpoledist_plot(
    runs,
    ci_type=None,
    cj_type=None,
    ci_subtype=None,
    cj_subtype=None,
    depth=None,
    nbins=30,
):
    # Make the figure
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.set_yscale("log")
    ax.grid(True)

    # Collect the distances
    dists = {}
    for i, (name, run) in enumerate(runs.items()):
        mask = make_mask(run, ci_type, cj_type, ci_subtype, cj_subtype, depth)

        # Ensure we only have pair tasks (i.e. the string "pair" is in the
        # task label)
        mask = np.logical_and(
            mask, np.array(["pair" in t for t in run.task_labels])
        )

        # Get the distances
        dists[name] = run.mpole_dists[mask]

    # Combine all distances
    all_dists = np.concatenate(list(dists.values()))

    # Anything to do here?
    if all_dists.size == 0:
        return

    # Construct the bins
    all_dists = np.concatenate(list(dists.values()))
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
    fig.savefig(filename, bbox_inches="tight")
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

    # Parse the arguments
    args = parser.parse_args()
    files = args.files
    labels = args.labels

    if len(labels) == 0:
        labels = files
        print("Using filenames as labels")

    if len(files) != len(labels):
        raise ValueError("Number of files and labels must match")

    # Parse all the task files
    runs = {}
    for f, lab in zip(files, labels):
        runs[lab] = TaskParser(f, lab)

    make_task_hist_split(runs)
    make_task_hist_time_split(runs)
    make_task_hist_time_split(runs, order_by_count=False)

    make_task_hist(runs)
    make_task_hist(runs, ci_type=1, cj_type=1)
    make_task_hist(runs, ci_type=2, cj_type=2)
    make_task_hist(runs, ci_type=3, cj_type=3)
    make_task_hist(runs, ci_type=1, cj_type=2)
    make_task_hist(runs, ci_type=1, cj_type=3)
    make_task_hist(runs, ci_type=2, cj_type=3)

    make_task_hist(runs, depth=0)
    make_task_hist(runs, ci_type=1, cj_type=1, depth=0)
    make_task_hist(runs, ci_type=2, cj_type=2, depth=0)
    make_task_hist(runs, ci_type=3, cj_type=3, depth=0)
    make_task_hist(runs, ci_type=1, cj_type=3, depth=0)
    make_task_hist(runs, ci_type=1, cj_type=2, depth=0)
    make_task_hist(runs, ci_type=2, cj_type=3, depth=0)

    make_task_hist_time_weighted(runs)
    make_task_hist_time_weighted(runs, ci_type=1, cj_type=1)
    make_task_hist_time_weighted(runs, ci_type=2, cj_type=2)
    make_task_hist_time_weighted(runs, ci_type=3, cj_type=3)
    make_task_hist_time_weighted(runs, ci_type=1, cj_type=3)
    make_task_hist_time_weighted(runs, ci_type=1, cj_type=2)
    make_task_hist_time_weighted(runs, ci_type=2, cj_type=3)

    make_task_hist_time_weighted(runs, depth=0)
    make_task_hist_time_weighted(runs, ci_type=1, cj_type=1, depth=0)
    make_task_hist_time_weighted(runs, ci_type=2, cj_type=2, depth=0)
    make_task_hist_time_weighted(runs, ci_type=3, cj_type=3, depth=0)
    make_task_hist_time_weighted(runs, ci_type=1, cj_type=3, depth=0)
    make_task_hist_time_weighted(runs, ci_type=1, cj_type=2, depth=0)
    make_task_hist_time_weighted(runs, ci_type=2, cj_type=3, depth=0)

    make_pair_mindist_plot(runs)
    make_pair_mindist_plot(runs, ci_type=1, cj_type=1)
    make_pair_mindist_plot(runs, ci_type=2, cj_type=2)
    make_pair_mindist_plot(runs, ci_type=3, cj_type=3)
    make_pair_mindist_plot(runs, ci_type=1, cj_type=3)
    make_pair_mindist_plot(runs, ci_type=1, cj_type=2)
    make_pair_mindist_plot(runs, ci_type=2, cj_type=3)

    make_pair_mpoledist_plot(runs)
    make_pair_mpoledist_plot(runs, ci_type=1, cj_type=1)
    make_pair_mpoledist_plot(runs, ci_type=2, cj_type=2)
    make_pair_mpoledist_plot(runs, ci_type=3, cj_type=3)
    make_pair_mpoledist_plot(runs, ci_type=1, cj_type=3)
    make_pair_mpoledist_plot(runs, ci_type=1, cj_type=2)
    make_pair_mpoledist_plot(runs, ci_type=2, cj_type=3)
