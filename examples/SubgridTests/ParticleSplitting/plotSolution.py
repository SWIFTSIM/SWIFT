"""
Plots the result from the splitting binary tree
in the simple example.
"""

from swiftsimio import load
import matplotlib.pyplot as plt
import numpy as np
import sys

try:
    plt.style.use("../../../tools/stylesheets/mnras.mplstyle")
except:
    pass


def add_arrow(
    line, position=None,
):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    line.axes.annotate(
        "",
        xytext=(xdata[0], ydata[0]),
        xy=(xdata.mean(), ydata.mean()),
        arrowprops=dict(arrowstyle="->", color=color),
        size=10,
    )


data = load("particleSplitting_0001.hdf5")

have_split = data.gas.split_counts > 0
counts = data.gas.split_counts[have_split]
split_trees = data.gas.split_trees[have_split]
formatted_split_trees = np.array(
    [f"{tree:b}".zfill(count) for tree, count in zip(split_trees.v, counts)],
    dtype=object,
)

special_coordinates = data.gas.coordinates[have_split].value.T

fig, ax = plt.subplots(figsize=(4, 4))

for particle, tree in enumerate(formatted_split_trees):
    ax.text(
        special_coordinates[0][particle],
        special_coordinates[1][particle],
        tree,
        ha="center",
        va="top",
        zorder=20,
    )

    for generation, item in enumerate(tree):
        if item == "0":
            continue
        else:
            parent = list(tree)
            parent[generation] = "0"
            parent = "".join(parent)
            which_parent = formatted_split_trees == parent

            (line,) = ax.plot(
                [
                    special_coordinates[0][which_parent],
                    special_coordinates[0][particle],
                ],
                [
                    special_coordinates[1][which_parent],
                    special_coordinates[1][particle],
                ],
                color=f"C{generation + 1}",
            )

            add_arrow(line)

            break

ax.scatter(special_coordinates[0], special_coordinates[1], zorder=10)

from matplotlib.lines import Line2D

custom_lines = [
    Line2D([0], [0], color=f"C{generation + 1}") for generation in range(len(tree))
]
custom_labels = [
    f"Generation {generation}" for generation in reversed(range(len(tree)))
]

ax.legend(custom_lines, custom_labels)

fig.subplots_adjust(0, 0, 1, 1)
ax.axis("off")

plt.savefig("particle_splitting.png")
