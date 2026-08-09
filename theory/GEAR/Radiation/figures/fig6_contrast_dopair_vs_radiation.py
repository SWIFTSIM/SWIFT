"""Figure 6: classic two-sided DOPAIR vs the one-sided radiation gather.

Left: the classic hydro sorted pair kernel (runner_doiact_functions_hydro.h)
walks BOTH cells' sort arrays in lockstep along their shared sid, using
rshift = shift . runner_shift[sid] to align the two cells' frames -- this
requires ci and cj to be sorted along the pair's own sid, which only their
own pair tasks populate, so it is adjacency-bound by construction.

Right: the radiation gather (runner_dopair_stars_hii_ionization_feedback,
:967-1083) only ever reads cj's sort array; the star contributes a single
projected position, not a second sorted list. It prefers the carried sid but
falls back to any axis cj is sorted along, and to an O(N) naive scan as a
last resort.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from style import (
    BLACK,
    CANDIDATE_COLOR,
    GRAY,
    STAR_COLOR,
    WIRED_COLOR,
    draw_arrow,
    draw_box,
    save_figure,
    star_marker,
)


def _sort_ticks(
    ax, x0: float, x1: float, y: float, color: str, n: int = 7, seed: int = 0
):
    rng = np.random.default_rng(seed)
    xs = np.sort(rng.uniform(x0 + 0.05, x1 - 0.05, n))
    for x in xs:
        ax.plot([x, x], [y - 0.06, y + 0.06], color=color, lw=1.4, zorder=4)
    return xs


def make_figure(output_path: str) -> None:
    """Render and save fig6 (PNG + PDF) to ``output_path`` (no extension)."""
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(8.0, 5.2), dpi=150)
    fig.patch.set_facecolor("white")

    # --- Left: classic two-sided lockstep DOPAIR -----------------------------
    w = 1.4
    draw_box(axl, (0, 1.0), w, w, edgecolor=WIRED_COLOR, lw=1.7)
    draw_box(axl, (w, 1.0), w, w, edgecolor=CANDIDATE_COLOR, lw=1.7)
    axl.text(
        w / 2,
        1.0 + w + 0.12,
        "ci",
        ha="center",
        fontsize=10,
        color=WIRED_COLOR,
        fontweight="bold",
    )
    axl.text(
        w + w / 2,
        1.0 + w + 0.12,
        "cj",
        ha="center",
        fontsize=10,
        color=CANDIDATE_COLOR,
        fontweight="bold",
    )
    axl.plot([w, w], [1.0, 1.0 + w], color=BLACK, lw=3.0, zorder=6)
    axl.text(w, 1.0 + w + 0.35, "shared face", ha="center", fontsize=8, color=BLACK)

    y_sort = 0.55
    xs_i = _sort_ticks(axl, 0, w, y_sort, WIRED_COLOR, seed=1)
    xs_j = _sort_ticks(axl, w, 2 * w, y_sort, CANDIDATE_COLOR, seed=2)
    axl.plot([0, 2 * w], [y_sort, y_sort], color=GRAY, lw=0.8, zorder=1)
    axl.text(
        -0.08, y_sort, "sort_i", ha="right", va="center", fontsize=8, color=WIRED_COLOR
    )
    axl.text(
        2 * w + 0.08,
        y_sort,
        "sort_j",
        ha="left",
        va="center",
        fontsize=8,
        color=CANDIDATE_COLOR,
    )

    # Lockstep pointers advancing inward from both ends.
    pid = xs_i[-2]
    pjd = xs_j[1]
    draw_arrow(
        axl,
        (pid, y_sort - 0.22),
        (pjd, y_sort - 0.22),
        color=BLACK,
        lw=1.2,
        connectionstyle="arc3,rad=-0.25",
        zorder=5,
    )
    axl.text(
        (pid + pjd) / 2,
        y_sort - 0.55,
        "di, dj advance inward together\n(rshift aligns the two frames)",
        ha="center",
        fontsize=7.6,
    )

    axl.text(
        w,
        -0.35,
        "two-sided: BOTH sort_i and sort_j read;\n"
        "only sorted along sids ci/cj's own pair tasks requested\n"
        "-> adjacency-bound",
        ha="center",
        va="top",
        fontsize=8.4,
    )

    axl.set_xlim(-0.9, 2 * w + 0.9)
    axl.set_ylim(-0.9, 1.0 + w + 0.65)
    axl.set_aspect("equal")
    axl.axis("off")
    axl.set_title("Classic DOPAIR (hydro)", fontsize=10.5)

    # --- Right: one-sided star-centric gather --------------------------------
    draw_box(axr, (0, 1.0), w, w, edgecolor=WIRED_COLOR, lw=1.7)
    draw_box(axr, (w + 0.9, 1.0), w, w, edgecolor=CANDIDATE_COLOR, lw=1.7)
    axr.text(
        w / 2,
        1.0 + w + 0.12,
        "ci",
        ha="center",
        fontsize=10,
        color=WIRED_COLOR,
        fontweight="bold",
    )
    axr.text(
        w + 0.9 + w / 2,
        1.0 + w + 0.12,
        "cj",
        ha="center",
        fontsize=10,
        color=CANDIDATE_COLOR,
        fontweight="bold",
    )
    axr.text(
        w + 0.45,
        1.0 + w + 0.35,
        "gap (need not be adjacent)",
        ha="center",
        fontsize=7.8,
        color=GRAY,
    )

    star_xy = (0.55 * w, 1.0 + 0.5 * w)
    star_marker(axr, star_xy, size=200)

    y_sort = 0.55
    xs_j = _sort_ticks(axr, w + 0.9, 2 * w + 0.9, y_sort, CANDIDATE_COLOR, seed=3)
    axr.plot([w + 0.9, 2 * w + 0.9], [y_sort, y_sort], color=GRAY, lw=0.8, zorder=1)
    axr.text(
        2 * w + 0.9 + 0.08,
        y_sort,
        "sort_j only",
        ha="left",
        va="center",
        fontsize=8,
        color=CANDIDATE_COLOR,
    )

    di_star = 0.9  # A point off to the left of cj's sort strip.
    axr.plot([di_star], [y_sort], marker="v", ms=9, color=STAR_COLOR, zorder=6)
    draw_arrow(
        axr,
        star_xy,
        (di_star, y_sort + 0.15),
        color=STAR_COLOR,
        lw=1.1,
        connectionstyle="arc3,rad=0.2",
        zorder=5,
    )
    axr.text(
        di_star,
        y_sort - 0.28,
        "di_star\n(one projected value)",
        ha="center",
        va="top",
        fontsize=7.6,
        color=STAR_COLOR,
    )

    axr.text(
        w + 0.45,
        -0.35,
        "one-sided: only cj's sort array is read;\n"
        "star contributes position + dx_max_sort only.\n"
        "carried sid -> any sid cj is sorted along -> naive O(N)",
        ha="center",
        va="top",
        fontsize=8.4,
    )

    axr.set_xlim(-0.9, 2 * w + 0.9 + 0.9)
    axr.set_ylim(-0.9, 1.0 + w + 0.65)
    axr.set_aspect("equal")
    axr.axis("off")
    axr.set_title("Radiation gather (this code)", fontsize=10.5)

    fig.suptitle(
        "Fig. 6 -- Two-sided lockstep vs one-sided star-centric scan",
        fontsize=11.5,
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, output_path)


if __name__ == "__main__":
    make_figure("fig6_contrast_dopair_vs_radiation")
