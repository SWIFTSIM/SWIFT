import numpy as np
import argparse

# import hardcoded data
from swift_hardcoded_data import TASKTYPES, SUBTYPES

parser = argparse.ArgumentParser(description="Plot task graphs")
parser.add_argument("input", help="Thread data file (-y output)")
args = parser.parse_args()

step = args.input

minticmaxtoc = np.loadtxt(
    "thread_info-step{0}.dat".format(step), max_rows=1, dtype="l", usecols=[4, 5]
)
tasks = np.loadtxt(
    "thread_info-step{0}.dat".format(step),
    skiprows=1,
    dtype={
        "names": (
            "rid",
            "type",
            "subtype",
            "is_pair",
            "tic",
            "toc",
            "nhydro_i",
            "nhydro_j",
            "ngrav_i",
            "ngrav_j",
            "sid",
        ),
        "formats": ("i", "i", "i", "i", "l", "l", "i", "i", "i", "i", "i"),
    },
)
tasklines = np.loadtxt(
    "task_line-step{0}.dat".format(step),
    dtype={"names": ("tid", "line"), "formats": ("U14", "i")},
)

mintic = minticmaxtoc[0]
maxtoc = minticmaxtoc[1]
tottoc = maxtoc - mintic

taskgraph = {}
for taskline in tasklines:
    if taskline["line"] < 0:
        # implicit task
        taskgraph[taskline["tid"]] = {
            "cost": 0,
            "parents": [],
            "cumulcost": 0,
            "name": "implicit_task",
            "rid": -1,
            "tic": -1,
            "toc": -1,
        }
    else:
        task = tasks[taskline["line"]]
        taskgraph[taskline["tid"]] = {
            "cost": task["toc"] - task["tic"],
            "parents": [],
            "cumulcost": 0,
            "name": "{0}/{1}".format(
                TASKTYPES[task["type"]], SUBTYPES[task["subtype"]]
            ),
            "rid": task["rid"],
            "tic": (task["tic"] - mintic) / tottoc,
            "toc": (task["toc"] - mintic) / tottoc,
        }

graph = np.loadtxt(
    "task_graph-step{0}.dat".format(step),
    dtype={"names": ("pid", "cid"), "formats": ("U14", "U14")},
)
maxcumulcost = 0
for connection in graph:
    cid = connection["cid"]
    pid = connection["pid"]
    taskgraph[cid]["parents"].append(pid)
    taskgraph[cid]["cumulcost"] = max(
        taskgraph[cid]["cumulcost"], taskgraph[pid]["cost"]
    )
    if taskgraph[cid]["cumulcost"] > maxcumulcost:
        maxcumulcost = taskgraph[cid]["cumulcost"]
        maxcumultid = cid

path = [maxcumultid]
parents = taskgraph[maxcumultid]["parents"]
while len(parents) > 0:
    maxcumulpid = parents[0]
    maxcumulcost = taskgraph[maxcumulpid]["cumulcost"]
    for p in parents:
        pcost = taskgraph[p]["cumulcost"]
        if pcost > maxcumulcost:
            maxcumulcost = pcost
            maxcumulpid = p
    path.append(p)
    parents = taskgraph[p]["parents"]

path = np.array(path)
for step in path[::-1]:
    task = taskgraph[step]
    print(task["rid"], task["name"], task["tic"], task["toc"])
