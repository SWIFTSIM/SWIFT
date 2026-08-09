"""Figure 4: the sorted-window prune is conservative for ANY sorted axis.

For a unit projection axis n (one of the 13 `runner_shift` directions), the
projected separation |((x_star - x_gas)).n| never exceeds the true 3D
separation |x_star - x_gas| (Cauchy-Schwarz on a unit vector). So a true
neighbour (r <= R) always falls inside the [di_star - R, di_star + R] window
on EVERY axis; the window can admit extra, farther candidates that the exact
r^2 < R^2 test then rejects, but it can never wrongly exclude a true one.
This is what justifies falling back to any sorted axis of cj
(runner_radiation_feedback.c:744-757) instead of only the carried sid.
"""

from __future__ import annotations

import numpy as np
import matplotlib.lines as mlines
import matplotlib.pyplot as plt

from style import (
    BLACK,
    CANDIDATE_COLOR,
    GRAY,
    STAR_COLOR,
    new_figure,
    save_figure,
    star_marker,
)

AXIS_COLOR_X = "#0072B2"
AXIS_COLOR_DIAG = "#CC79A7"


def make_figure(output_path: str) -> None:
    """Render and save fig4 (PNG + PDF) to ``output_path`` (no extension)."""
    fig, ax = new_figure((6.6, 5.8))

    star = np.array([0.6, 3.2])
    gas = np.array([4.4, 0.9])
    sep = gas - star
    r = np.linalg.norm(sep)

    n_x = np.array([1.0, 0.0])
    n_diag = np.array([1.0, 1.0]) / np.sqrt(2)

    proj_x = abs((gas - star) @ n_x)
    proj_diag = abs((gas - star) @ n_diag)

    ax.plot([star[0], gas[0]], [star[1], gas[1]], color=BLACK, lw=1.8, zorder=5)
    star_marker(ax, star, size=260)
    ax.scatter(
        *gas,
        marker="o",
        s=100,
        color=CANDIDATE_COLOR,
        edgecolor=BLACK,
        lw=0.9,
        zorder=6,
    )
    mid = (star + gas) / 2
    ax.text(
        mid[0] - 0.75, mid[1] + 0.25, f"r = {r:.2f}", fontsize=9.5, fontweight="bold"
    )

    # x-axis projection: drop verticals onto a horizontal strip below.
    y_x = -0.55
    for p in (star, gas):
        ax.plot([p[0], p[0]], [p[1], y_x], color=AXIS_COLOR_X, lw=0.8, ls=":", zorder=2)
    ax.annotate(
        "",
        xy=(gas[0], y_x),
        xytext=(star[0], y_x),
        arrowprops=dict(arrowstyle="-", color=AXIS_COLOR_X, lw=3.0),
        zorder=4,
    )
    ax.text(
        (star[0] + gas[0]) / 2,
        y_x - 0.3,
        f"x-axis projection = {proj_x:.2f}",
        ha="center",
        va="top",
        fontsize=8.8,
        color=AXIS_COLOR_X,
        zorder=4,
    )

    # Diagonal-sid axis: the line through the star along n_diag, drawn only
    # as far as the projection itself needs (short reference segment, not a
    # line spanning the whole figure). The star already lies on its own
    # axis by construction, so only the projected segment is shown.
    t_gas = (gas - star) @ n_diag
    axis_line = (
        star[None, :] + np.array([-0.35, t_gas + 0.35])[:, None] * n_diag[None, :]
    )
    ax.plot(
        axis_line[:, 0],
        axis_line[:, 1],
        color=AXIS_COLOR_DIAG,
        lw=0.9,
        ls="--",
        zorder=1,
        alpha=0.6,
    )

    foot_gas = star + t_gas * n_diag
    ax.annotate(
        "",
        xy=foot_gas,
        xytext=star,
        arrowprops=dict(arrowstyle="-", color=AXIS_COLOR_DIAG, lw=3.0),
    )
    ax.text(
        foot_gas[0] + 0.18,
        foot_gas[1] + 0.1,
        f"diagonal-sid\nprojection = {proj_diag:.2f}",
        fontsize=8.8,
        color=AXIS_COLOR_DIAG,
        ha="left",
        va="bottom",
    )

    ax.text(
        1.6,
        -1.55,
        f"Both projections ({proj_x:.2f}, {proj_diag:.2f}) <= true r ({r:.2f}): the slab never\n"
        "excludes a real neighbour, on ANY sorted axis. The exact r^2 < R^2\n"
        "test is what actually decides membership.",
        ha="center",
        va="top",
        fontsize=8.8,
        bbox=dict(boxstyle="round", fc="white", ec=GRAY, lw=0.8),
        zorder=3,
    )

    legend_handles = [
        mlines.Line2D([], [], color=BLACK, lw=1.8, label="true separation r"),
        mlines.Line2D([], [], color=AXIS_COLOR_X, lw=3.0, label="x-axis projection"),
        mlines.Line2D(
            [], [], color=AXIS_COLOR_DIAG, lw=3.0, label="diagonal-sid projection"
        ),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8.4, frameon=False)

    ax.set_xlim(-1.8, 5.4)
    ax.set_ylim(-3.4, 6.1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(
        "Fig. 4 -- Any-axis conservativeness of the sorted-window prune",
        fontsize=11,
    )

    save_figure(fig, output_path)


if __name__ == "__main__":
    make_figure("fig4_any_axis_conservative")
