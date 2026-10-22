"""Which time bins are active at which point of the integer time-line.

A bin is active at time t when t is a multiple of the length of its step.
The chart is drawn with the functions of src/timeline.h.
"""

import os

if os.path.exists("active_bins.png"):
    # do not generate the plot again
    exit()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from timeline_helpers import (
    PALETTE,
    get_integer_timestep,
    get_integer_time_end,
    get_max_active_bin,
)

n_bins = 8
# Every time a step can end at is a multiple of the shortest step (4 ticks).
times = np.arange(4, 520, 4)

active = np.zeros((n_bins, len(times)), dtype=bool)
for j, t in enumerate(times):
    max_active = get_max_active_bin(int(t))
    for b in range(1, n_bins + 1):
        # Same statement, two ways: end of the step of the bin equals t,
        # and the bin is not above the highest active bin.
        by_end = get_integer_time_end(int(t), b) == t
        assert by_end == (b <= max_active)
        active[b - 1, j] = by_end

fig, ax = plt.subplots(
    2, 1, figsize=(10, 5.6), sharex=True, gridspec_kw={"height_ratios": [4, 1]},
    constrained_layout=True,
)
for b in range(1, n_bins + 1):
    for j, t in enumerate(times):
        if active[b - 1, j]:
            ax[0].add_patch(
                plt.Rectangle(
                    (t - 2, b - 0.4), 4, 0.8, color=PALETTE["blue"], lw=0
                )
            )
        else:
            ax[0].add_patch(
                plt.Rectangle((t - 2, b - 0.4), 4, 0.8, color="0.92", lw=0)
            )
ax[0].set_xlim(0, times[-1] + 2)
ax[0].set_ylim(0.4, n_bins + 1.0)
ax[0].set_yticks(range(1, n_bins + 1))
ax[0].set_yticklabels(
    ["bin %d (%d ticks)" % (b, get_integer_timestep(b)) for b in range(1, n_bins + 1)]
)
ax[0].set_title("Bins that end a step (blue) at each integer time t")

max_bins = [get_max_active_bin(int(t)) for t in times]
ax[1].step(times, max_bins, where="mid", color=PALETTE["vermillion"], lw=1.4)
ax[1].set_ylabel("highest\nactive bin")
ax[1].set_xlabel("Integer time t [ticks]")
ax[1].set_yticks([1, 2, 4, 6, 8])
ax[1].grid(True, color="0.85")

fig.savefig("active_bins.png", dpi=140)
