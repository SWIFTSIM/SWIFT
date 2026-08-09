"""Figure 3 (KEY): sorted-array scan of a non-adjacent receiver cell C.

Demonstrates that the sorted-window prune operates on a GLOBAL 1D axis of
ABSOLUTE projections (runner_sort.c: entries[k].d = px . runner_shift[sid],
px = parts[k].x, an absolute position). The star is projected onto that same
axis (runner_radiation_feedback.c:998-1004); adjacency between A and C plays
no role -- only whether the window [di_star - R - dx_max, di_star + R +
dx_max] intersects C's own d-range.

Panel (a) is the micro-example from the README: star at x=1.2, C spanning
[3, 4], R=2.1, dx_max=0.05 -> window [-0.95, 3.35]. Panel (b) shows the
periodic-wrap variant, where the star's image is shifted onto C's side.
"""

from __future__ import annotations

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

from style import (
    BLACK,
    CANDIDATE_COLOR,
    GRAY,
    STAR_COLOR,
    WINDOW_COLOR,
    draw_box,
    new_figure,
    save_figure,
    star_marker,
)


def _axis_panel(
    ax,
    *,
    star_x: float,
    c_lo: float,
    c_hi: float,
    r: float,
    dx_max: float,
    title: str,
    star_label: str,
    wrapped_ghost: float | None = None,
):
    """Draw one cell row + global-axis-projection panel."""
    a_lo, a_hi = star_x - 1.0, star_x - 0.05
    b_lo, b_hi = a_hi, c_lo

    row_y = 1.0
    draw_box(ax, (a_lo, row_y), a_hi - a_lo, 0.8, lw=1.4)
    draw_box(ax, (b_lo, row_y), b_hi - b_lo, 0.8, lw=1.4)
    draw_box(ax, (c_lo, row_y), c_hi - c_lo, 0.8, lw=1.9, edgecolor=CANDIDATE_COLOR)
    ax.text(
        (a_lo + a_hi) / 2,
        row_y + 0.4,
        "A",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        (b_lo + b_hi) / 2,
        row_y + 0.4,
        "B",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        (c_lo + c_hi) / 2,
        row_y + 0.4,
        "C",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        color=CANDIDATE_COLOR,
    )

    star_marker(ax, (star_x, row_y + 0.4), size=180)

    # Global axis below the row.
    axis_y = 0.0
    axis_lo, axis_hi = a_lo - 0.3, c_hi + 0.3
    ax.plot([axis_lo, axis_hi], [axis_y, axis_y], color=BLACK, lw=1.2, zorder=2)
    ax.annotate(
        "",
        xy=(axis_hi, axis_y),
        xytext=(axis_hi - 0.01, axis_y),
        arrowprops=dict(arrowstyle="-|>", color=BLACK, lw=1.2),
    )
    ax.text(axis_hi + 0.05, axis_y, "d (global 1D axis)", va="center", fontsize=8)

    # C's gas particles as absolute-d ticks.
    n_ticks = 11
    ticks = np.linspace(c_lo, c_hi, n_ticks)
    di_min = star_x - r - dx_max
    di_max = star_x + r + dx_max
    for d in ticks:
        scanned = d <= di_max
        color = CANDIDATE_COLOR if scanned else GRAY
        ax.plot([d, d], [axis_y - 0.08, axis_y + 0.08], color=color, lw=1.6, zorder=4)
        ax.plot([d], [axis_y], marker="o", ms=3.2, color=color, zorder=5)

    # Star's projection di_star.
    ax.plot([star_x], [axis_y], marker="v", ms=9, color=STAR_COLOR, zorder=6)
    ax.text(
        star_x,
        axis_y - 0.22,
        star_label,
        ha="center",
        va="top",
        fontsize=8,
        color=STAR_COLOR,
    )

    # Window band.
    ax.axvspan(di_min, di_max, color=WINDOW_COLOR, alpha=0.22, zorder=1)
    ax.plot(
        [di_min, di_min],
        [axis_y - 0.5, axis_y + 1.9],
        color=WINDOW_COLOR,
        lw=1.0,
        ls=":",
        zorder=1,
    )
    ax.plot(
        [di_max, di_max],
        [axis_y - 0.5, axis_y + 1.9],
        color=WINDOW_COLOR,
        lw=1.0,
        ls=":",
        zorder=1,
    )
    ax.text(
        di_max,
        axis_y - 0.34,
        f"di_max={di_max:g}",
        ha="center",
        va="top",
        fontsize=7.5,
        color=WINDOW_COLOR,
    )
    ax.text(
        di_min,
        axis_y - 0.34,
        f"di_min={di_min:g}",
        ha="center",
        va="top",
        fontsize=7.5,
        color=WINDOW_COLOR,
    )

    if wrapped_ghost is not None:
        star_marker(ax, (wrapped_ghost, row_y + 0.4), size=140)
        ax.plot(
            [wrapped_ghost],
            [axis_y],
            marker="v",
            ms=9,
            color=STAR_COLOR,
            alpha=0.5,
            zorder=6,
        )
        ax.annotate(
            "",
            xy=(wrapped_ghost, row_y + 0.15),
            xytext=(star_x, row_y + 0.15),
            arrowprops=dict(
                arrowstyle="-|>",
                color=STAR_COLOR,
                lw=1.1,
                ls="--",
                connectionstyle="arc3,rad=-0.3",
            ),
        )
        ax.text(
            (star_x + wrapped_ghost) / 2,
            row_y - 0.55,
            "pix = si->x - shift",
            ha="center",
            fontsize=7.5,
            color=STAR_COLOR,
        )

    ax.set_title(title, fontsize=9.5)
    ax.set_xlim(axis_lo - 0.1, axis_hi + 1.0)
    ax.set_ylim(-0.75, row_y + 1.05)
    ax.set_aspect("auto")
    ax.axis("off")


def make_figure(output_path: str) -> None:
    """Render and save fig3 (PNG + PDF) to ``output_path`` (no extension)."""
    fig = plt.figure(figsize=(7.0, 5.4), dpi=150)
    fig.patch.set_facecolor("white")
    ax1 = fig.add_axes((0.07, 0.565, 0.86, 0.30))
    ax2 = fig.add_axes((0.07, 0.115, 0.86, 0.30))

    _axis_panel(
        ax1,
        star_x=1.2,
        c_lo=3.0,
        c_hi=4.0,
        r=2.1,
        dx_max=0.05,
        title="(a) non-wrapped: shift = 0, di_star = si->x directly",
        star_label="di_star = 1.2",
    )

    # Wrapped variant: box of size L=4.6, star near the right edge, C near
    # the left edge -- shift = -L brings the star's image next to C.
    box_len = 4.6
    star_x_wrapped = box_len - 0.3  # 4.3, near the right edge.
    ghost = star_x_wrapped - box_len  # -0.3, image on C's side.
    _axis_panel(
        ax2,
        star_x=ghost,
        c_lo=-2.3,
        c_hi=-1.3,
        r=1.6,
        dx_max=0.05,
        title="(b) periodic wrap: shift = -L folds the star's real position\n"
        "(dashed star, right) onto the image used for di_star (solid, left)",
        star_label="di_star (wrapped)",
        wrapped_ghost=None,
    )
    # Draw the real (unwrapped) star position too, off to the right, with an
    # arrow showing the -shift fold used to compute pix = si->x - shift.
    real_x = ghost + box_len
    star_marker(ax2, (real_x, 1.4), size=140)
    ax2.plot([real_x], [0.0], marker="v", ms=8, color=STAR_COLOR, alpha=0.4, zorder=6)
    ax2.annotate(
        "",
        xy=(ghost, 0.0),
        xytext=(real_x, 0.0),
        arrowprops=dict(
            arrowstyle="-|>",
            color=STAR_COLOR,
            lw=1.1,
            ls="--",
            connectionstyle="arc3,rad=0.25",
        ),
    )
    ax2.text(
        (ghost + real_x) / 2,
        1.9,
        "real position si->x (unwrapped)",
        ha="center",
        fontsize=7.5,
        color=STAR_COLOR,
        alpha=0.7,
    )
    ax2.set_xlim(-3.0, real_x + 0.6)

    legend_handles = [
        mlines.Line2D(
            [],
            [],
            color=CANDIDATE_COLOR,
            marker="|",
            linestyle="",
            markersize=12,
            markeredgewidth=2,
            label="C tick: d <= di_max -- scanned",
        ),
        mlines.Line2D(
            [],
            [],
            color=GRAY,
            marker="|",
            linestyle="",
            markersize=12,
            markeredgewidth=2,
            label="C tick: d > di_max -- loop already stopped",
        ),
        mlines.Line2D(
            [],
            [],
            color=STAR_COLOR,
            marker="v",
            linestyle="",
            markersize=8,
            label="di_star (star's projection)",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=1,
        fontsize=8.2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.0),
    )

    fig.text(
        0.5,
        0.985,
        "Fig. 3 -- Sorted scan against a non-adjacent receiver C\n"
        "adjacency is irrelevant: only window/d-range overlap on the global axis decides what is scanned",
        ha="center",
        va="top",
        fontsize=11,
    )
    save_figure(fig, output_path)


if __name__ == "__main__":
    make_figure("fig3_sorted_scan_nonadjacent")
