import numpy as np
import argparse

# import hardcoded data
from swift_hardcoded_data import TASKTYPES, SUBTYPES

parser = argparse.ArgumentParser(description="Plot task graphs")
parser.add_argument("input", help="Thread data file (-y output)")
args = parser.parse_args()

step = args.input

minticmaxtoc = np.loadtxt(
    "thread_info-step{0}.dat".format(step), max_rows=1, dtype="l", usecols=[4, 5, 10]
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
            "id",
        ),
        "formats": ("i", "i", "i", "i", "l", "l", "i", "i", "i", "i", "i", "i"),
    },
)

mintic = minticmaxtoc[0]
maxtoc = minticmaxtoc[1]
cpufreq = minticmaxtoc[2] * 0.001
tottoc = maxtoc - mintic

taskgraph = {}
for task in tasks:
    taskgraph[task["id"]] = {
        "cost": task["toc"] - task["tic"],
        "parents": [],
        "children": [],
        "cost_to_task": 0,
        "cost_from_task": 0,
        "first_ancestor": 0,
        "last_offspring": 0,
        "name": "{0}/{1}".format(TASKTYPES[task["type"]], SUBTYPES[task["subtype"]]),
        "rid": task["rid"],
        "tic": (task["tic"] - mintic),  # / tottoc,
        "toc": (task["toc"] - mintic),  # / tottoc,
    }

graph = np.loadtxt(
    "task_graph-step{0}.dat".format(step),
    dtype={"names": ("pid", "cid"), "formats": ("l", "l")},
)
for connection in graph:
    cid = connection["cid"]
    pid = connection["pid"]
    taskgraph[cid]["parents"].append(pid)
    taskgraph[pid]["children"].append(cid)


def update_children(task, taskgraph, cost):
    last_offspring = taskgraph[task]["toc"]
    if len(taskgraph[task]["children"]) > 0:
        for child in taskgraph[task]["children"]:
            taskgraph[child]["cost_to_task"] = max(
                taskgraph[child]["cost_to_task"], cost
            )
            update_children(child, taskgraph, cost + taskgraph[child]["cost"])
            last_offspring = max(last_offspring, taskgraph[child]["last_offspring"])
    taskgraph[task]["last_offspring"] = last_offspring


def update_parents(task, taskgraph, cost):
    first_ancestor = taskgraph[task]["tic"]
    if len(taskgraph[task]["parents"]) > 0:
        for parent in taskgraph[task]["parents"]:
            taskgraph[parent]["cost_from_task"] = max(
                taskgraph[parent]["cost_from_task"], cost
            )
            update_parents(parent, taskgraph, cost + taskgraph[parent]["cost"])
            first_ancestor = min(first_ancestor, taskgraph[parent]["first_ancestor"])
    taskgraph[task]["first_ancestor"] = first_ancestor


first_tasks = []
last_tasks = []
for task in taskgraph:
    if len(taskgraph[task]["parents"]) == 0:
        # follow the trail
        update_children(task, taskgraph, taskgraph[task]["cost"])
        first_tasks.append(task)
    if len(taskgraph[task]["children"]) == 0:
        update_parents(task, taskgraph, taskgraph[task]["cost"])
        last_tasks.append(task)

for task in taskgraph:
    taskgraph[task]["cost_to_task"] += taskgraph[task]["cost"]
    taskgraph[task]["cost_from_task"] += taskgraph[task]["cost"]

maxfirst = 0
maxi = 0
print("First tasks:")
maxtime = 0
for first in first_tasks:
    maxtime = max(
        maxtime, taskgraph[first]["last_offspring"] - taskgraph[first]["first_ancestor"]
    )
    print(first, taskgraph[first])
    maxthis = taskgraph[first]["cost_from_task"]
    if maxthis > maxfirst:
        maxfirst = maxthis
        maxi = first
print(maxtime / cpufreq)

next = maxi
print(next, taskgraph[next]["name"])
while len(taskgraph[next]["children"]) > 0:
    maxi = 0
    maxcost = 0
    for child in taskgraph[next]["children"]:
        maxthis = taskgraph[child]["cost_from_task"]
        if maxthis > maxcost:
            maxcost = maxthis
            maxi = child
    print(maxi, taskgraph[maxi]["name"])
    next = maxi

maxlast = 0
maxi = 0
print("Last tasks:")
maxtime = 0
for last in last_tasks:
    maxtime = max(
        maxtime, taskgraph[last]["last_offspring"] - taskgraph[last]["first_ancestor"]
    )
    print(last, taskgraph[last])
    maxthis = taskgraph[last]["cost_to_task"]
    if maxthis > maxlast:
        maxlast = maxthis
        maxi = last
print(maxtime / cpufreq)

next = maxi
print(next, taskgraph[next]["name"])
while len(taskgraph[next]["parents"]) > 0:
    maxi = 0
    maxcost = 0
    for parent in taskgraph[next]["parents"]:
        maxthis = taskgraph[parent]["cost_to_task"]
        if maxthis > maxcost:
            maxcost = maxthis
            maxi = parent
    print(maxi, taskgraph[maxi]["name"])
    next = maxi

print(maxfirst / cpufreq, maxlast / cpufreq)
