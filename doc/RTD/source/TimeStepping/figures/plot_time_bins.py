"""Time bins and the length of the time-steps they stand for.

Uses the definitions of src/timeline.h (see timeline_helpers.py). The example
run length and the dt_min/dt_max limits are those of the Sedov blast example.
"""

import os

if os.path.exists("time_bins.png"):
    # do not generate the plot again
    exit()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from timeline_helpers import (
    NUM_TIME_BINS,
    MAX_NR_TIMESTEPS,
    PALETTE,
    get_integer_timestep,
    get_timestep,
)

# Sedov blast example: TimeIntegration in examples/HydroTests/SedovBlast_3D
t_begin, t_end, dt_min, dt_max = 0.0, 5e-2, 1e-9, 1e-2
time_base = (t_end - t_begin) / MAX_NR_TIMESTEPS

bins = list(range(1, NUM_TIME_BINS + 1))
ticks = [get_integer_timestep(b) for b in bins]
phys = [get_timestep(b, time_base) for b in bins]


def largest_step_below(limit):
    """The loop that engine_config() uses to report the limits (in float32)."""
    dt = np.float32(t_end - t_begin)
    while dt > np.float32(limit):
        dt = np.float32(dt / np.float32(2.0))
    return float(dt)


# The two values printed at the start of a run ("on time-line") and their bins.
step_low, step_high = largest_step_below(dt_min), largest_step_below(dt_max)
bin_low = min(bins, key=lambda b: abs(phys[b - 1] / step_low - 1.0))
bin_high = min(bins, key=lambda b: abs(phys[b - 1] / step_high - 1.0))

fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)

# Left: integer ticks
ax[0].semilogy(bins, ticks, "o-", color=PALETTE["blue"], ms=4, lw=1.2)
ax[0].set_xlabel("Time bin")
ax[0].set_ylabel("Length of the time-step [integer ticks]")
ax[0].set_title("Integer time-line")
ax[0].annotate(
    "bin 1: 4 ticks",
    (1, ticks[0]),
    (3, 1e6),
    arrowprops=dict(arrowstyle="->", color="0.3"),
    color="0.2",
)
ax[0].annotate(
    "bin 56: $2^{57}$ ticks\n(the whole time-line)",
    (56, ticks[-1]),
    (12, 1e14),
    arrowprops=dict(arrowstyle="->", color="0.3"),
    color="0.2",
)
ax[0].grid(True, which="major", color="0.85")

# Right: physical length for the example run
ax[1].semilogy(bins, phys, "o-", color=PALETTE["orange"], ms=4, lw=1.2)
ax[1].axhspan(dt_min, dt_max, color=PALETTE["sky"], alpha=0.25, lw=0)
ax[1].axhline(dt_max, color=PALETTE["blue"], ls="--", lw=1)
ax[1].axhline(dt_min, color=PALETTE["blue"], ls="--", lw=1)
ax[1].text(2, dt_max * 1.6, "dt_max = %g" % dt_max, color=PALETTE["blue"])
ax[1].text(2, dt_min * 1.6, "dt_min = %g" % dt_min, color=PALETTE["blue"])
ax[1].set_xlabel("Time bin")
ax[1].set_ylabel("Length of the time-step [internal time units]")
ax[1].plot([bin_low, bin_high], [phys[bin_low - 1], phys[bin_high - 1]], "o",
           color=PALETTE["vermillion"], ms=8, mfc="none", mew=2,
           label="limits reported by SWIFT (bins %d and %d)" % (bin_low, bin_high))
ax[1].legend(loc="lower right", fontsize=9)
ax[1].set_title("Example run of length %g" % (t_end - t_begin), fontsize=11)
ax[1].grid(True, which="major", color="0.85")

fig.savefig("time_bins.png", dpi=140)
