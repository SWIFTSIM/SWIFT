"""Kick-drift-kick timeline of one particle over two steps, to scale.

The particle takes a step of bin 3 (16 ticks) followed by a step of bin 4
(32 ticks). The growth from bin 3 to bin 4 is only possible because the end of
the first step is a multiple of 32 ticks; this is checked with the rule of
make_integer_timestep() (see timeline_helpers.py).
"""

import os

if os.path.exists("kdk_timeline.png"):
    # do not generate the plot again
    exit()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from timeline_helpers import (
    PALETTE,
    apply_growth_rules,
    get_integer_time_begin,
    get_integer_time_end,
    get_integer_timestep,
)

bin_a, bin_b = 3, 4
t0 = 16  # start of the first step
dti_a = get_integer_timestep(bin_a)
t1 = t0 + dti_a
dti_b = get_integer_timestep(bin_b)
t2 = t1 + dti_b

# The first step is a valid step of bin 3, and its end allows the growth.
assert get_integer_time_begin(t0 + 1, bin_a) == t0
assert get_integer_time_end(t1, bin_a) == t1
assert apply_growth_rules(64, bin_a, t1) == dti_b

kick1_col, kick2_col = PALETTE["blue"], PALETTE["orange"]
fig, ax = plt.subplots(figsize=(10, 4.6), constrained_layout=True)

rows = {"steps": 3.0, "kick": 2.0, "drift": 1.0, "v": 0.0}


def bar(a, b, y, colour, label=None, hatch=None):
    ax.barh(y, b - a, left=a, height=0.5, color=colour, edgecolor="white", lw=1,
            label=label, hatch=hatch)


for (beg, dti, b), name in (((t0, dti_a, bin_a), "A"), ((t1, dti_b, bin_b), "B")):
    half = beg + dti // 2
    end = beg + dti
    ax.barh(rows["steps"], dti, left=beg, height=0.5, color="0.9",
            edgecolor="0.3", lw=1.2)
    ax.text(beg + dti / 2, rows["steps"],
            "step %s: bin %d, %d ticks" % (name, b, dti), ha="center", va="center")
    bar(beg, half, rows["kick"], kick1_col,
        "kick1 (first half)" if name == "A" else None)
    bar(half, end, rows["kick"], kick2_col,
        "kick2 (second half)" if name == "A" else None)
    ax.annotate("", (end, rows["drift"]), (beg, rows["drift"]),
                arrowprops=dict(arrowstyle="->", lw=2, color="0.25"))
    ax.text(beg + dti / 2, rows["drift"] + 0.12, "drift: $x = x + v_{full}\\,\\Delta t$",
            ha="center", va="bottom", fontsize=9)
    # velocity time-stamps: after kick1 at the half-step, after kick2 at the end
    ax.plot([half], [rows["v"]], "o", color=kick1_col, ms=8)
    ax.plot([end], [rows["v"]], "o", mfc="white", color=kick2_col, ms=8, mew=2)

ax.plot([], [], "o", color=kick1_col, ms=8, label="$v_{full}$ after kick1")
ax.plot([], [], "o", mfc="white", color=kick2_col, ms=8, mew=2,
        label="$v_{full}$ after kick2")

ax.set_yticks(list(rows.values()))
ax.set_yticklabels(["time-steps", "kicks", "drift", "velocity\ntime-stamp"])
ax.set_xticks([t0, t0 + dti_a // 2, t1, t1 + dti_b // 2, t2])
ax.set_xlim(t0 - 2, t2 + 2)
ax.set_ylim(-1.5, 4.3)
ax.set_xlabel("Integer time [ticks]")
ax.grid(True, axis="x", color="0.85")
ax.axvline(t1, color="0.3", ls=":", lw=1)
ax.text(t1 + 0.6, 4.2,
        "at this tick, in this order:\nforce, kick2 (A), timestep, kick1 (B)",
        va="top", fontsize=9)
ax.legend(loc="lower center", ncol=4, fontsize=9, framealpha=0.95)
fig.savefig("kdk_timeline.png", dpi=140)
