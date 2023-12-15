#!/usr/bin/env python3

# ------------------------------------------------------------
# This file contains data that is hardcoded into swift
# that needs to be reproduced exactly in order for analysis
# outputs to make any sense.
# The data is used in other scripts in this directory, this
# script is intended for imports only.
# ------------------------------------------------------------

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
