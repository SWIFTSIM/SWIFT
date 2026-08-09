"""Figure 1: the region-wiring traversal structure.

Draws a star's leaf cell A inside its ``radiation_level`` region, the self
link descending that region's whole subtree, and the (2D analogue of the)
26 pair links, each descending an independently-wired neighbour region down
to leaves under the per-subtree box prune. See README.md fig1 section for
file:line anchors.
"""

from __future__ import annotations

import matplotlib.lines as mlines
import matplotlib.pyplot as plt

from style import (
    BLACK,
    CANDIDATE_COLOR,
    GRAY,
    STAR_COLOR,
    WIRED_COLOR,
    draw_arrow,
    draw_box,
    new_figure,
    save_figure,
    star_marker,
)

REGION = 2.0  # Side length of one radiation_level region.
LEAF = REGION / 2.0  # 2x2 leaves per region (2D stand-in for an octree).
GAP = 0.55


def _draw_region_with_leaves(ax, x0: float, y0: float, *, highlight: bool):
    """Draw one radiation_level region subdivided into its 2x2 leaves."""
    edge = WIRED_COLOR if highlight else GRAY
    lw = 2.2 if highlight else 1.1
    draw_box(ax, (x0, y0), REGION, REGION, edgecolor=edge, lw=lw, zorder=3)
    for i in range(2):
        for j in range(2):
            draw_box(
                ax,
                (x0 + i * LEAF, y0 + j * LEAF),
                LEAF,
                LEAF,
                edgecolor=GRAY,
                lw=0.7,
                ls=":",
                zorder=2,
            )


def make_figure(output_path: str) -> None:
    """Render and save fig1 (PNG + PDF) to ``output_path`` (no extension)."""
    fig, ax = new_figure((6.6, 7.2))

    centers = [
        (gi, gj, gi * (REGION + GAP), gj * (REGION + GAP))
        for gi in range(3)
        for gj in range(3)
    ]
    mid = (1, 1)
    example = (2, 1)

    for gi, gj, x0, y0 in centers:
        _draw_region_with_leaves(ax, x0, y0, highlight=(gi, gj) == mid)

    cx0, cy0 = REGION + GAP, REGION + GAP  # Center region lower-left.
    ax.text(
        cx0 + REGION / 2,
        cy0 + REGION + 0.08,
        "radiation_level region\n(contains star's leaf A)",
        ha="center",
        va="bottom",
        fontsize=7.8,
        color=WIRED_COLOR,
    )
    ex_x0, ex_y0 = example[0] * (REGION + GAP), example[1] * (REGION + GAP)
    ax.text(
        ex_x0 + REGION / 2,
        ex_y0 + REGION + 0.08,
        "one of up to 26 neighbour regions",
        ha="center",
        va="bottom",
        fontsize=7.8,
        color=CANDIDATE_COLOR,
    )

    star_xy = (cx0 + 0.32 * LEAF, cy0 + 1.65 * LEAF)
    star_marker(ax, star_xy, size=240)
    ax.text(
        star_xy[0] + 0.10,
        star_xy[1],
        "star\n(leaf A)",
        fontsize=7.5,
        color=STAR_COLOR,
        va="center",
    )

    # Self link: descends the WHOLE subtree of the star's own region.
    region_center = (cx0 + REGION / 2, cy0 + REGION / 2)
    for i in range(2):
        for j in range(2):
            leaf_c = (cx0 + (i + 0.5) * LEAF, cy0 + (j + 0.5) * LEAF)
            draw_arrow(ax, region_center, leaf_c, color=WIRED_COLOR, lw=1.0, zorder=4)

    # Pair links: one arrow from the center region to EACH of the other 8
    # region boxes (2D analogue of up to 26 in 3D).
    for gi, gj, x0, y0 in centers:
        if (gi, gj) == mid:
            continue
        ncenter = (x0 + REGION / 2, y0 + REGION / 2)
        draw_arrow(
            ax,
            region_center,
            ncenter,
            color=CANDIDATE_COLOR,
            lw=1.0,
            connectionstyle="arc3,rad=0.05",
            zorder=3,
        )

    # For the example neighbour: show its own independent recursion to leaves.
    ex_center = (ex_x0 + REGION / 2, ex_y0 + REGION / 2)
    for i in range(2):
        for j in range(2):
            leaf_c = (ex_x0 + (i + 0.5) * LEAF, ex_y0 + (j + 0.5) * LEAF)
            draw_arrow(ax, ex_center, leaf_c, color=CANDIDATE_COLOR, lw=0.8, zorder=4)

    legend_handles = [
        mlines.Line2D(
            [],
            [],
            color=WIRED_COLOR,
            marker=">",
            linestyle="-",
            lw=1.4,
            label="self link -- descends own subtree",
        ),
        mlines.Line2D(
            [],
            [],
            color=CANDIDATE_COLOR,
            marker=">",
            linestyle="-",
            lw=1.4,
            label="pair link -- visits neighbour, recurses independently",
        ),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.02),
        fontsize=7.6,
        frameon=False,
        handlelength=2.2,
    )

    ax.set_xlim(-0.4, 3 * (REGION + GAP) - GAP + 0.4)
    ax.set_ylim(-0.5, 3 * (REGION + GAP) - GAP + 0.55)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(
        "Fig. 1 -- Traversal structure\n"
        "no neighbour-of-neighbour hops: visits follow region wiring only",
        fontsize=10.5,
        pad=10,
    )

    save_figure(fig, output_path)


if __name__ == "__main__":
    make_figure("fig1_traversal_structure")
