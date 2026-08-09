"""Figure 2: Case 1, |A|B|C| with both B and C wired and independently tested.

Star sits in leaf A; the search radius overlaps into B and C. B and C are
each reached through their OWN link in radiation_level->stars.radiation_in
and tested independently (no B-to-C shadowing or ordering). Candidates found
in either cell merge into one global distance-sorted buffer.
"""

from __future__ import annotations

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from style import (
    BLACK,
    CANDIDATE_COLOR,
    REDDISH_PURPLE,
    WINDOW_COLOR,
    WIRED_COLOR,
    draw_arrow,
    draw_box,
    new_figure,
    save_figure,
    star_marker,
)

B_COLOR = WIRED_COLOR
C_COLOR = REDDISH_PURPLE


def make_figure(output_path: str) -> None:
    """Render and save fig2 (PNG + PDF) to ``output_path`` (no extension)."""
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(6.4, 6.6), dpi=150, gridspec_kw={"height_ratios": [2.1, 1.0]}
    )
    for ax in (ax_top, ax_bot):
        ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    # --- Top panel: cell row + search radius -------------------------------
    w, h = 1.2, 1.2
    a0, b0, c0 = 0.0, w, 2 * w
    for x0, name in ((a0, "A"), (b0, "B"), (c0, "C")):
        draw_box(ax_top, (x0, 0), w, h, lw=1.6, zorder=4)
        ax_top.text(
            x0 + w / 2,
            -0.16,
            name,
            ha="center",
            va="top",
            fontsize=10,
            fontweight="bold",
        )

    star_xy = (a0 + 0.75 * w, 0.5 * h)
    radius = 1.85

    circle = mpatches.Circle(
        star_xy, radius, fill=False, ls="--", lw=1.3, edgecolor=WINDOW_COLOR, zorder=5
    )
    ax_top.add_patch(circle)
    shade_clip = mpatches.Circle(star_xy, radius, transform=ax_top.transData)
    shade = mpatches.Rectangle(
        (b0, 0), 2 * w, h, facecolor=WINDOW_COLOR, alpha=0.22, zorder=1
    )
    shade.set_clip_path(shade_clip)
    ax_top.add_patch(shade)

    star_marker(ax_top, star_xy, size=260)
    ax_top.text(
        star_xy[0],
        star_xy[1] + 0.18,
        "star",
        ha="center",
        fontsize=8.5,
        color="#D55E00",
    )

    draw_arrow(
        ax_top,
        star_xy,
        (b0 + w / 2, h + 0.55),
        color=B_COLOR,
        lw=1.2,
        connectionstyle="arc3,rad=-0.15",
        zorder=6,
    )
    draw_arrow(
        ax_top,
        star_xy,
        (c0 + w / 2, h + 0.55),
        color=C_COLOR,
        lw=1.2,
        connectionstyle="arc3,rad=-0.25",
        zorder=6,
    )
    ax_top.text(
        b0 + w / 2,
        h + 0.62,
        "link -> B",
        ha="center",
        va="bottom",
        fontsize=8,
        color=B_COLOR,
    )
    ax_top.text(
        c0 + w / 2,
        h + 0.62,
        "link -> C",
        ha="center",
        va="bottom",
        fontsize=8,
        color=C_COLOR,
    )

    ax_top.set_xlim(-0.3, 3 * w + 0.3)
    ax_top.set_ylim(-0.4, h + 1.15)
    ax_top.set_aspect("equal")
    ax_top.axis("off")
    ax_top.set_title(
        "Fig. 2 -- Case 1: both B and C wired, tested independently\n"
        "search radius overlaps B and C; neither test depends on the other's outcome",
        fontsize=10.2,
    )

    # --- Bottom panel: merge into one distance-sorted buffer ----------------
    rng = np.random.default_rng(7)
    n_b, n_c = 5, 4
    r2_b = np.sort(rng.uniform(0.05, 0.9, n_b))
    r2_c = np.sort(rng.uniform(0.3, 1.0, n_c))
    entries = sorted(
        [(r, "B") for r in r2_b] + [(r, "C") for r in r2_c], key=lambda t: t[0]
    )
    max_ngbs_demo = 6

    x_positions = np.arange(len(entries))
    for i, (r2, origin) in enumerate(entries):
        color = B_COLOR if origin == "B" else C_COLOR
        kept = i < max_ngbs_demo
        draw_box(
            ax_bot,
            (i, 0),
            0.85,
            0.85,
            edgecolor=color,
            facecolor=color,
            alpha=0.85 if kept else 0.18,
            lw=1.2,
        )
        ax_bot.text(
            i + 0.425,
            0.425,
            origin,
            ha="center",
            va="center",
            fontsize=9,
            color="white" if kept else color,
            fontweight="bold",
        )

    ax_bot.axvline(max_ngbs_demo - 0.075, color=BLACK, ls=":", lw=1.2)
    ax_bot.text(
        max_ngbs_demo - 0.075, 1.05, "nearest-max_ngbs cut", ha="center", fontsize=8
    )
    ax_bot.annotate(
        "",
        xy=(len(entries) - 0.3, 0.425),
        xytext=(-0.35, 0.425),
        arrowprops=dict(arrowstyle="-|>", color=BLACK, lw=1.0),
    )
    ax_bot.text(
        len(entries) / 2,
        -0.32,
        "increasing r2 -- one global buffer, sorted across the B/C union",
        ha="center",
        fontsize=8.5,
    )

    legend_handles = [
        mlines.Line2D(
            [],
            [],
            color=B_COLOR,
            marker="s",
            linestyle="",
            markersize=9,
            label="candidate from B",
        ),
        mlines.Line2D(
            [],
            [],
            color=C_COLOR,
            marker="s",
            linestyle="",
            markersize=9,
            label="candidate from C",
        ),
    ]
    ax_bot.legend(handles=legend_handles, loc="upper right", fontsize=8, frameon=False)

    ax_bot.set_xlim(-0.6, len(entries) + 0.3)
    ax_bot.set_ylim(-0.55, 1.3)
    ax_bot.set_aspect("equal")
    ax_bot.axis("off")

    fig.tight_layout()
    save_figure(fig, output_path)


if __name__ == "__main__":
    make_figure("fig2_case1_independent_merge")
