#!/usr/bin/python
"""
Usage:
  ./plot_task_level.py task_level.txt

Description:
  Plot the number of tasks for each depth level and each type of task.
"""


import pandas as pd
import matplotlib.pyplot as plt
import sys

# get filename
filename = sys.argv[-1]

# Column names
names = ["type", "subtype", "depth", "count"]

# read file
data = pd.read_csv(filename, sep=" ", comment="#", names=names)

# generate color map
cmap = plt.get_cmap("hsv")
N = data["depth"].max() + 5

# plot data
for i in range(data["depth"].max()):
    ind = data["depth"] == i
    label = "depth = %i" % i
    c = cmap(i / N)
    plt.plot(
        data["type"][ind] + "_" + data["subtype"][ind],
        data["count"][ind],
        ".",
        label=label,
        color=c,
    )

# modify figure parameters and show it
plt.gca().set_yscale("log")
plt.xticks(rotation=45)
plt.ylabel("Number of Tasks")
plt.gcf().subplots_adjust(bottom=0.15)
plt.legend()
plt.show()
