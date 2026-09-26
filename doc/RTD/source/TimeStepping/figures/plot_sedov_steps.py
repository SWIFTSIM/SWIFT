"""Step table of a real run: the Sedov blast example.

Reads the timesteps.txt written by SWIFT (columns as in the header of the
file) and shows how the length of the global step, the number of updated
particles and the range of active bins evolve.
"""

import os

if os.path.exists("sedov_steps.png") and os.path.exists("sedov_step_cost.png"):
    # do not generate the plots again
    exit()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from timeline_helpers import PALETTE

data = np.loadtxt("sedov_3d_timesteps.txt", comments="#")
step, time, dt = data[:, 0], data[:, 1], data[:, 4]
min_bin, max_bin = data[:, 5], data[:, 6]
updates, wall_ms, props = data[:, 7], data[:, 12], data[:, 13].astype(int)

n_part = int(updates.max())  # the first step updates every particle
t_end = time[-1]
snapshot = (props & 16) > 0  # engine_step_prop_snapshot
rebuild = ((props & 1) > 0) & ~snapshot  # engine_step_prop_rebuild
all_active = updates == n_part
body = step > 0  # step 0 is the initialisation

# Cost of the run against a run in which every particle uses the shortest step.
dt_shortest = dt[body].min()
n_steps_single = t_end / dt_shortest
updates_single = n_steps_single * n_part
updates_real = updates[body].sum()

fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True, constrained_layout=True)

ax[0].step(time[body], dt[body], where="post", color=PALETTE["blue"], lw=1.4)
ax[0].set_yscale("log")
ax[0].set_ylabel("Length of the\nglobal step")
ax[0].grid(True, which="major", color="0.88")

ax[1].scatter(time[body & ~all_active], updates[body & ~all_active] / n_part,
              s=8, color=PALETTE["blue"], label="some particles active")
ax[1].scatter(time[all_active & body], updates[all_active & body] / n_part,
              s=28, color=PALETTE["vermillion"], marker="D",
              label="all particles active (synchronised step)")
ax[1].set_yscale("log")
ax[1].set_ylabel("Fraction of the\nparticles updated")
ax[1].legend(loc="lower right", fontsize=9)
ax[1].grid(True, which="major", color="0.88")

ax[2].fill_between(time[body], min_bin[body], max_bin[body], step="post",
                   color=PALETTE["sky"], alpha=0.6, lw=0)
ax[2].step(time[body], min_bin[body], where="post", color=PALETTE["blue"], lw=1,
           label="lowest active bin")
ax[2].step(time[body], max_bin[body], where="post", color=PALETTE["vermillion"],
           lw=1, label="highest active bin")
ax[2].set_ylabel("Active bins")
ax[2].set_xlabel("Time [internal units]")
ax[2].legend(loc="upper left", fontsize=9)
ax[2].grid(True, which="major", color="0.88")

fig.suptitle(
    "Sedov blast, %d particles, %d steps. Particle updates: %.1f N instead of %d N"
    " with one global step" % (n_part, body.sum(), updates_real / n_part,
                                round(n_steps_single)),
    fontsize=11,
)
fig.savefig("sedov_steps.png", dpi=140)

# Cost of a step
fig, ax = plt.subplots(figsize=(7.5, 4.6), constrained_layout=True)
regular = body & ~rebuild & ~snapshot
ax.scatter(updates[regular], wall_ms[regular], s=10,
           color=PALETTE["blue"], label="regular step")
ax.scatter(updates[body & snapshot], wall_ms[body & snapshot], s=30,
           color=PALETTE["purple"], marker="D", label="step that writes a snapshot")
ax.scatter(updates[body & rebuild], wall_ms[body & rebuild], s=22,
           color=PALETTE["vermillion"], marker="s", label="step with a tree rebuild")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Particles updated in the step")
ax.set_ylabel("Wall-clock time of the step [ms]")
ax.grid(True, which="major", color="0.88")
ax.legend(fontsize=9)
fig.savefig("sedov_step_cost.png", dpi=140)
print("shortest step %.4e, single-step run would need %d steps" % (dt_shortest, round(n_steps_single)))
print("updates: real %.3e, single-step %.3e, ratio %.1f" % (updates_real, updates_single, updates_single / updates_real))
