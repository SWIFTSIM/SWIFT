#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# This file contains data that is hardcoded into swift that needs to be
# reproduced exactly in order for analysis outputs to make any sense. The data
# is used in other scripts in this directory, this script is intended for
# imports only.
# Additionally, we add some setup work:
# - check if a file "task_labels_task_types.txt" exists. If it does, verify that
#   the hardcoded data in this file corresponds to the output written by SWIFT
#   in that file.
# - set up a list of useful task type/subtype combinations. Note that this list
#   needs to be verified manually for completeness.
# - assing colours to all task types, subtypes, and combinations thereof.
# -----------------------------------------------------------------------------

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
    "drift_sink",
    "drift_bpart",
    "drift_gpart",
    "drift_gpart_out",
    "hydro_end_force",
    "kick1",
    "kick2",
    "timestep",
    "timestep_limiter",
    "timestep_sync",
    "collect",
    "send",
    "recv",
    "pack",
    "unpack",
    "grav_long_range",
    "grav_mm",
    "grav_down_in",
    "grav_down",
    "grav_end_force",
    "cooling",
    "cooling_in",
    "cooling_out",
    "star_formation",
    "star_formation_in",
    "star_formation_out",
    "star_formation_sink",
    "csds",
    "stars_in",
    "stars_out",
    "stars_ghost_in",
    "stars_density_ghost",
    "stars_ghost_out",
    "stars_prep_ghost1",
    "hydro_prep_ghost1",
    "stars_prep_ghost2",
    "stars_sort",
    "stars_resort",
    "bh_in",
    "bh_out",
    "bh_density_ghost",
    "bh_swallow_ghost1",
    "bh_swallow_ghost2",
    "bh_swallow_ghost3",
    "fof_self",
    "fof_pair",
    "fof_attach_self",
    "fof_attach_pair",
    "neutrino_weight",
    "sink_in",
    "sink_density_ghost",
    "sink_ghost1",
    "sink_ghost2",
    "sink_out",
    "rt_in",
    "rt_out",
    "sink_formation",
    "rt_ghost1",
    "rt_ghost2",
    "rt_transport_out",
    "rt_tchem",
    "rt_advance_cell_time",
    "rt_sort",
    "rt_collect_times",
    #  "count",
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
    "part_swallow",
    "bpart_merger",
    "gpart",
    "spart_density",
    "part_prep1",
    "spart_prep2",
    "stars_density",
    "stars_prep1",
    "stars_prep2",
    "stars_feedback",
    "sf_counts",
    "bpart_rho",
    "bpart_feedback",
    "bh_density",
    "bh_swallow",
    "do_gas_swallow",
    "do_bh_swallow",
    "bh_feedback",
    "sink_do_sink_swallow",
    "sink_swallow",
    "sink_do_gas_swallow",
    "rt_gradient",
    "rt_transport",
    #  "count",
]

#  Task/subtypes of interest.
FULLTYPES = [
    "self/limiter",
    "self/force",
    "self/gradient",
    "self/density",
    "self/grav",
    "sub_self/limiter",
    "sub_self/force",
    "sub_self/gradient",
    "sub_self/density",
    "pair/limiter",
    "pair/force",
    "pair/gradient",
    "pair/density",
    "pair/grav",
    "sub_pair/limiter",
    "sub_pair/force",
    "sub_pair/gradient",
    "sub_pair/density",
    "recv/xv",
    "send/xv",
    "recv/rho",
    "send/rho",
    "recv/tend_part",
    "send/tend_part",
    "recv/tend_gpart",
    "send/tend_gpart",
    "recv/tend_spart",
    "send/tend_spart",
    "recv/tend_bpart",
    "send/tend_bpart",
    "recv/gpart",
    "send/gpart",
    "recv/spart",
    "send/spart",
    "send/sf_counts",
    "recv/sf_counts",
    "recv/bpart",
    "send/bpart",
    "recv/limiter",
    "send/limiter",
    "pack/limiter",
    "unpack/limiter",
    "self/stars_density",
    "pair/stars_density",
    "sub_self/stars_density",
    "sub_pair/stars_density",
    "self/stars_prep1",
    "pair/stars_prep1",
    "sub_self/stars_prep1",
    "sub_pair/stars_prep1",
    "self/stars_prep2",
    "pair/stars_prep2",
    "sub_self/stars_prep2",
    "sub_pair/stars_prep2",
    "self/stars_feedback",
    "pair/stars_feedback",
    "sub_self/stars_feedback",
    "sub_pair/stars_feedback",
    "self/bh_density",
    "pair/bh_density",
    "sub_self/bh_density",
    "sub_pair/bh_density",
    "self/bh_swallow",
    "pair/bh_swallow",
    "sub_self/bh_swallow",
    "sub_pair/bh_swallow",
    "self/do_swallow",
    "pair/do_swallow",
    "sub_self/do_swallow",
    "sub_pair/do_swallow",
    "self/bh_feedback",
    "pair/bh_feedback",
    "sub_self/bh_feedback",
    "sub_pair/bh_feedback",
    "self/rt_gradient",
    "pair/rt_gradient",
    "sub_self/rt_gradient",
    "sub_pair/rt_gradient",
    "self/rt_transport",
    "pair/rt_transport",
    "sub_self/rt_transport",
    "sub_pair/rt_transport",
    "self/sink_density",
    "pair/sink_density",
    "sub_self/sink_density",
    "sub_pair/sink_density",
    "self/sink_swallow",
    "pair/sink_swallow",
    "sub_self/sink_swallow",
    "sub_pair/sink_swallow",
    "self/sink_do_swallow",
    "pair/sink_do_swallow",
    "sub_self/sink_do_swallow",
    "sub_pair/sink_do_swallow",
    "self/sink_do_gas_swallow",
    "pair/sink_do_gas_swallow",
    "sub_self/sink_do_gas_swallow",
    "sub_pair/sink_do_gas_swallow",
]

#  A number of colours for the various types. Recycled when there are
#  more task types than colours...
colours = [
    "cyan",
    "lightgray",
    "darkblue",
    "yellow",
    "tan",
    "dodgerblue",
    "sienna",
    "aquamarine",
    "bisque",
    "blue",
    "green",
    "lightgreen",
    "brown",
    "purple",
    "moccasin",
    "olivedrab",
    "chartreuse",
    "olive",
    "darkgreen",
    "green",
    "mediumseagreen",
    "mediumaquamarine",
    "darkslategrey",
    "mediumturquoise",
    "black",
    "cadetblue",
    "skyblue",
    "red",
    "slategray",
    "gold",
    "slateblue",
    "blueviolet",
    "mediumorchid",
    "firebrick",
    "magenta",
    "hotpink",
    "pink",
    "orange",
    "lightgreen",
]
maxcolours = len(colours)

#  Set colours of task/subtype.
TASKCOLOURS = {}
ncolours = 0
for task in TASKTYPES:
    TASKCOLOURS[task] = colours[ncolours]
    ncolours = (ncolours + 1) % maxcolours

SUBCOLOURS = {}
for task in FULLTYPES:
    SUBCOLOURS[task] = colours[ncolours]
    ncolours = (ncolours + 1) % maxcolours

for task in SUBTYPES:
    SUBCOLOURS[task] = colours[ncolours]
    ncolours = (ncolours + 1) % maxcolours


# check if label files are found that have (possibly different) labels
# output by SWIFT itself
import os

if os.path.exists("task_labels_task_types.txt"):
    print("SWIFT task label file 'task_labels_task_types.txt' found, reading it.")
    with open("task_labels_task_types.txt", "r") as file:
        NEW_TASKTYPES = []
        for line in file.readlines()[1:]:
            NEW_TASKTYPES.append(line.split()[1])
    if not len(NEW_TASKTYPES) == len(TASKTYPES):
        print(
            "Hardcoded labels do not match SWIFT labels,"
            " consider updating the hardcoded labels (different list size)."
        )
    else:
        for ilabel in range(len(TASKTYPES)):
            if not NEW_TASKTYPES[ilabel] == TASKTYPES[ilabel]:
                print(
                    "Hardcoded labels do not match SWIFT labels,"
                    " consider updating the hardcoded labels "
                    " (different labels: TASKTYPES[{0}]: {1} vs {2}).".format(
                        ilabel, TASKTYPES[ilabel], NEW_TASKTYPES[ilabel]
                    )
                )
    TASKTYPES = NEW_TASKTYPES
if os.path.exists("task_labels_task_subtypes.txt"):
    print("SWIFT task label file 'task_labels_task_subtypes.txt' found, reading it.")
    with open("task_labels_task_subtypes.txt", "r") as file:
        NEW_SUBTYPES = []
        for line in file.readlines()[1:]:
            NEW_SUBTYPES.append(line.split()[1])
    if not len(NEW_SUBTYPES) == len(SUBTYPES):
        print(
            "Hardcoded labels do not match SWIFT labels,"
            " consider updating the hardcoded labels (different list size)."
        )
    else:
        for ilabel in range(len(SUBTYPES)):
            if not NEW_SUBTYPES[ilabel] == SUBTYPES[ilabel]:
                print(
                    "Hardcoded labels do not match SWIFT labels,"
                    " consider updating the hardcoded labels"
                    " (different labels: SUBTYPES[{0}]: {1} vs {2}).".format(
                        ilabel, SUBTYPES[ilabel], NEW_SUBTYPES[ilabel]
                    )
                )
    SUBTYPES = NEW_SUBTYPES

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
