"""Figure 5: Case 2, C two regions away when radiation_level sits at leaf level.

When radiation_level coarsening has not (yet) happened, each leaf is its own
radiation_level region, so the wired stencil is only the leaf's up-to-26
immediate neighbours (runner_radiation_feedback.c:266-274). A cell two
regions away is structurally absent from radiation_in -- not filtered by the
bounds check, simply never an edge in the graph -- even if the search radius
geometrically reaches into it. Three guards bound how bad this gets, and a
rebuild coarsens the stencil so C eventually enters it.
"""

from __future__ import annotations

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from style import (
    BLACK,
    GRAY,
    UNWIRED_COLOR,
    WINDOW_COLOR,
    WIRED_COLOR,
    draw_arrow,
    draw_box,
    new_figure,
    save_figure,
    star_marker,
)


def make_figure(output_path: str) -> None:
    """Render and save fig5 (PNG + PDF) to ``output_path`` (no extension)."""
    fig = plt.figure(figsize=(7.0, 7.6), dpi=150)
    fig.patch.set_facecolor("white")

    ax_main = fig.add_axes((0.06, 0.50, 0.90, 0.40))
    ax_before = fig.add_axes((0.10, 0.08, 0.35, 0.30))
    ax_after = fig.add_axes((0.55, 0.08, 0.35, 0.30))

    # --- Main panel: A|B|C, each its own radiation_level (leaf level) ------
    w = 1.3
    a0, b0, c0 = 0.0, w, 2 * w
    names = [("A", a0, WIRED_COLOR), ("B", b0, WIRED_COLOR), ("C", c0, UNWIRED_COLOR)]
    for name, x0, color in names:
        draw_box(ax_main, (x0, 0), w, w, edgecolor=color, lw=1.8, zorder=3)
        ax_main.text(
            x0 + w / 2,
            -0.22,
            name,
            ha="center",
            va="top",
            fontsize=11,
            fontweight="bold",
            color=color,
        )

    star_xy = (a0 + 0.7 * w, 0.55 * w)
    star_marker(ax_main, star_xy, size=240)

    radius = 2.05
    circle = mpatches.Circle(
        star_xy, radius, fill=False, ls="--", lw=1.3, edgecolor=WINDOW_COLOR, zorder=4
    )
    ax_main.add_patch(circle)

    draw_arrow(
        ax_main,
        star_xy,
        (b0 + w / 2, w + 0.5),
        color=WIRED_COLOR,
        lw=1.3,
        connectionstyle="arc3,rad=-0.15",
        zorder=5,
    )
    ax_main.text(
        b0 + w / 2,
        w + 0.58,
        "wired\n(1 region away)",
        ha="center",
        va="bottom",
        fontsize=7.8,
        color=WIRED_COLOR,
    )

    ax_main.text(
        c0 + w / 2,
        w + 0.58,
        "NOT wired\n(2 regions away --\nabsent from radiation_in,\nnot bounds-check-skipped)",
        ha="center",
        va="bottom",
        fontsize=7.8,
        color=UNWIRED_COLOR,
    )
    ax_main.plot(
        [c0 + w / 2],
        [w + 0.15],
        marker="x",
        ms=11,
        mew=2.2,
        color=UNWIRED_COLOR,
        zorder=5,
    )

    ax_main.set_xlim(-0.3, 3 * w + 0.3)
    ax_main.set_ylim(-0.55, w + 1.35)
    ax_main.set_aspect("equal")
    ax_main.axis("off")
    ax_main.set_title(
        "Fig. 5 -- Case 2: radiation_level at leaf level, C two regions away\n"
        "radius sliver-covers C, but C is structurally invisible to the star's traversal",
        fontsize=10.3,
    )

    guard_text = (
        "Three guards bound this gap:\n"
        "1. Split gate -- a region stops splitting once "
        "factor . kernel_gamma . h_hii_max >= 0.5 . dmin\n"
        "2. Rebuild criterion -- whole-tree rebuild fires once reach + drift > dmin\n"
        "3. In-pass clamp -- dynamic_search_radius capped at interaction_limit + dmin this pass;\n"
        "   the outer shell is deferred to the next rebuild, never lost"
    )
    fig.text(
        0.5,
        0.485,
        guard_text,
        ha="center",
        va="top",
        fontsize=8.6,
        bbox=dict(boxstyle="round", fc="white", ec=GRAY, lw=0.8),
    )

    # --- Before/after inset --------------------------------------------------
    for ax, title, wired_c in (
        (ax_before, "before rebuild", False),
        (ax_after, "after rebuild", True),
    ):
        ww = 0.55
        if not wired_c:
            for i, (name, color) in enumerate(
                zip("ABC", [WIRED_COLOR, WIRED_COLOR, UNWIRED_COLOR])
            ):
                draw_box(ax, (i * ww, 0), ww, ww, edgecolor=color, lw=1.6)
                ax.text(
                    i * ww + ww / 2,
                    ww + 0.08,
                    name,
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color=color,
                )
        else:
            draw_box(ax, (0, 0), 2 * ww, ww, edgecolor=WIRED_COLOR, lw=2.0)
            ax.plot([ww, ww], [0, ww], color=GRAY, lw=0.8, ls=":")
            draw_box(ax, (2 * ww, 0), ww, ww, edgecolor=WIRED_COLOR, lw=1.6)
            ax.text(
                ww,
                ww + 0.08,
                "A+B (coarsened)",
                ha="center",
                va="bottom",
                fontsize=7.6,
                color=WIRED_COLOR,
            )
            ax.text(
                2 * ww + ww / 2,
                ww + 0.08,
                "C",
                ha="center",
                va="bottom",
                fontsize=8,
                color=WIRED_COLOR,
            )
            draw_arrow(
                ax,
                (ww, ww / 2),
                (2 * ww + ww / 2, ww / 2),
                color=WIRED_COLOR,
                lw=1.1,
                connectionstyle="arc3,rad=-0.3",
                zorder=4,
            )
        ax.set_xlim(-0.15, 3 * ww + 0.15)
        ax.set_ylim(-0.15, ww + 0.35)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(title, fontsize=9)

    fig.text(
        0.5,
        0.055,
        "Rebuild coarsens radiation_level: A and B merge into one region "
        "whose 26-neighbour stencil now reaches C.",
        ha="center",
        va="top",
        fontsize=8.4,
        style="italic",
    )

    save_figure(fig, output_path)


if __name__ == "__main__":
    make_figure("fig5_case2_stencil_gap")
