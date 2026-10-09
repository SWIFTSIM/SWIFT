"""A particle interrupted by the time-step limiter, to scale.

A particle of bin 5 (step from 0 to 64 ticks, first half-kick already applied)
is woken up at t = 20 by an active neighbour of bin 1. The new bin is
-wakeup + 2 = 3. The interval arithmetic below is the one of
timestep_limit_part() in src/timestep_limiter.h (inactive-particle branch).
"""

import os

if os.path.exists("limiter_interrupt.png"):
    # do not generate the plot again
    exit()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from timeline_helpers import (
    PALETTE,
    get_integer_time_begin,
    get_integer_time_end,
    get_integer_timestep,
    get_max_active_bin,
)

ti_current = 20
old_bin = 5
waker_bin = 1
wakeup = -waker_bin  # what runner_iact_nonsym_limiter stores
new_bin = -wakeup + 2

ti_beg_old = get_integer_time_begin(ti_current, old_bin)
ti_end_old = get_integer_time_end(ti_current, old_bin)
dti_old = ti_end_old - ti_beg_old
dti_new = get_integer_timestep(new_bin)
k = 0
while ti_beg_old + k * dti_new <= ti_current:
    k += 1
ti_beg_new = ti_beg_old + (k - 1) * dti_new
ti_end_half_old = ti_beg_old + dti_old // 2
ti_end_half_new = ti_beg_new + dti_new // 2
apply_missing_kick1 = new_bin > get_max_active_bin(ti_current)
ti_end_new = ti_beg_new + dti_new if apply_missing_kick1 else ti_current + dti_new

# Values that appear in the caption of the documentation.
assert (ti_beg_old, ti_end_old, ti_beg_new, ti_end_half_new) == (0, 64, 16, 24)
assert apply_missing_kick1 and ti_end_new == 32

blue, orange, red = PALETTE["blue"], PALETTE["orange"], PALETTE["vermillion"]
fig, ax = plt.subplots(figsize=(10, 4.8), constrained_layout=True)


def bar(a, b, y, colour, hatch=None, alpha=1.0, label=None):
    ax.barh(y, b - a, left=a, height=0.5, color=colour, edgecolor="white", lw=1,
            hatch=hatch, alpha=alpha, label=label)


# Row 3: what was planned
ax.barh(3, dti_old, left=ti_beg_old, height=0.5, color="0.9", edgecolor="0.3")
ax.text(32, 3, "old step: bin %d, %d ticks" % (old_bin, dti_old), ha="center")
bar(ti_beg_old, ti_end_half_old, 2, blue, label="kick1 (applied)")
bar(ti_end_half_old, ti_end_old, 2, orange, alpha=0.35, hatch="//",
    label="kick2 (planned, never done)")

# Row 1: what happens at ti_current
bar(ti_beg_old, ti_end_half_old, 1, red, hatch="xx", alpha=0.6,
    label="kick1 undone (negative dt)")
bar(ti_beg_old, ti_beg_new, 0.4, blue, label="kick to the new start")
bar(ti_beg_new, ti_end_half_new, -0.2, blue, hatch="..",
    label="missing kick1 of the new step")

# Row -1: the new step
ax.barh(-1, dti_new, left=ti_beg_new, height=0.5, color="0.9", edgecolor="0.3")
ax.text(ti_beg_new + dti_new / 2, -1, "new step: bin %d, %d ticks" % (new_bin, dti_new),
        ha="center")

ax.axvline(ti_current, color="0.2", ls="--", lw=1.2)
ax.text(ti_current + 0.5, 4.4, "wake-up by a neighbour\nof bin %d at t = %d" % (waker_bin, ti_current),
        va="top", fontsize=9)

ax.set_yticks([3, 2, 1, 0.4, -0.2, -1])
ax.set_yticklabels(["before", "kicks before", "step 1: undo", "step 2: redo",
                    "step 3: complete", "after"])
ax.set_xticks([0, 16, 20, 24, 32, 48, 64])
ax.set_xlim(-2, 66)
ax.set_ylim(-1.6, 4.5)
ax.set_xlabel("Integer time [ticks]")
ax.grid(True, axis="x", color="0.88")
ax.legend(loc="lower right", fontsize=8.5, framealpha=0.95)
fig.savefig("limiter_interrupt.png", dpi=140)
