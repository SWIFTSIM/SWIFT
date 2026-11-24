"""Merge of the shares of one slot, and two beams that share a slot."""

import numpy as np
import isrf_style as st
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Wedge


def arrow(ax, p, q, col, lw=1.2, ls="-", ms=9):
    """Draw an arrow from p to q."""
    ax.add_patch(
        FancyArrowPatch(
            p,
            q,
            arrowstyle="-|>",
            mutation_scale=ms,
            color=col,
            lw=lw,
            ls=ls,
            shrinkA=0,
            shrinkB=0,
        )
    )


def panel_merge(ax) -> None:
    """Three shares in one slot and their energy-weighted mean origin."""
    j = np.array([0.0, 0.0])
    box = dict(fc="white", ec="none", pad=0.3)
    origins = [np.array([-4.2, 1.6]), np.array([-3.6, -0.7]), np.array([-2.9, 1.5])]
    energies = np.array([1.0, 0.5, 0.25])
    lo, hi = np.radians(150), np.radians(210)
    for t in (lo, hi):
        ax.plot(
            [0, 4.5 * np.cos(t)], [0, 4.5 * np.sin(t)], color=st.GREY, ls="--", lw=0.8
        )
    ax.add_patch(Wedge(j, 4.5, 150, 210, color=st.LIGHT, alpha=0.5, lw=0))
    ax.text(
        -4.6,
        -2.15,
        "origins whose light\nlands in one slot of $j$",
        color=st.GREY,
        fontsize=8,
    )

    for o, e in zip(origins, energies):
        ax.add_patch(Circle(o, 0.08, color=st.STAR, zorder=5))
        arrow(
            ax,
            o,
            j - 0.15 * o / np.linalg.norm(o),
            st.RECEIVER,
            lw=0.6 + 4.0 * e,
            ms=6 + 6 * e,
        )
        ax.text(*(o + [-0.1, 0.25]), f"$e_s={e:g}$", ha="center", fontsize=8)
    mean = (energies[:, None] * np.array(origins)).sum(0) / energies.sum()
    ax.add_patch(
        Circle(mean, 0.12, facecolor="white", edgecolor=st.BLACK, lw=1.2, zorder=6)
    )
    ax.text(
        *(mean + [0.1, -0.22]),
        "mean origin $\\boldsymbol{o}_k$",
        ha="center",
        va="top",
        fontsize=8,
        bbox=box,
    )
    arrow(ax, j, j - 0.45 * mean / np.linalg.norm(mean) * 3, st.BLACK, lw=1.3)
    ax.text(
        *(j - 1.55 * mean / np.linalg.norm(mean) + [0.15, 0.3]),
        r"$\boldsymbol{n}_k$",
        ha="left",
    )
    ax.add_patch(Circle(j, 0.1, color=st.RECEIVER, zorder=6))
    ax.text(0.1, -0.4, r"$j$", ha="center")
    d_mean = (energies * np.linalg.norm(origins, axis=1)).sum() / energies.sum()
    ax.text(
        -4.9,
        -2.75,
        f"$d_k={np.linalg.norm(mean):.2f}$ $\\leq$ mean distance ${d_mean:.2f}$",
        fontsize=8,
    )
    ax.set_xlim(-5.0, 1.7)
    ax.set_ylim(-3.0, 3.0)
    ax.set_title("(a) merge: one origin per slot")


def lobe(phi0: float, p: int = 4, scale: float = 2.2) -> np.ndarray:
    """Return a polar lobe r = scale cos^p of the angle to the direction phi0."""
    phi = np.linspace(phi0 - np.pi / 2, phi0 + np.pi / 2, 200)
    r = scale * np.cos(phi - phi0) ** p
    return np.c_[r * np.cos(phi), r * np.sin(phi)], phi


def panel_beams(ax) -> None:
    """Two beams at 60 degrees: the spread of beam 2 enters the slot of beam 1."""
    box = dict(fc="white", ec="none", pad=0.3)
    b1, b2 = 0.0, np.radians(60)
    lo, hi = np.radians(-38), np.radians(38)
    ax.add_patch(Wedge((0, 0), 2.9, -38, 38, color=st.LIGHT, alpha=0.5, lw=0))
    for t in (lo, hi):
        ax.plot(
            [0, 2.9 * np.cos(t)], [0, 2.9 * np.sin(t)], color=st.GREY, ls="--", lw=0.8
        )
    ax.text(2.25, -0.35, "slot of\nbeam 1", ha="center", color=st.GREY, fontsize=8)
    for phi0, col, name in [(b1, st.SKY, "beam 1"), (b2, st.PURPLE, "beam 2")]:
        pts, phi = lobe(phi0)
        ax.add_patch(
            Polygon(
                np.r_[[[0, 0]], pts],
                closed=True,
                facecolor=col,
                alpha=0.35,
                edgecolor=col,
                lw=1.0,
            )
        )
        tip = 2.45 * np.array([np.cos(phi0), np.sin(phi0)])
        arrow(ax, -tip * 0.6, -tip * 0.15, col, lw=1.6)
        ax.text(
            *(-tip * 0.75 + [0.0, 0.2 if phi0 > 0 else -0.3]),
            name,
            ha="center",
            color=col,
            bbox=box,
        )
    pts, phi = lobe(b2)
    inside = (phi > lo) & (phi < hi)
    ax.add_patch(
        Polygon(
            np.r_[[[0, 0]], pts[inside]],
            closed=True,
            facecolor="none",
            edgecolor=st.PURPLE,
            hatch="////",
            lw=0.8,
        )
    )
    ax.annotate(
        "part of beam 2\nin the slot of beam 1",
        xy=(1.05, 0.55),
        xytext=(0.2, 2.45),
        fontsize=8,
        ha="center",
        arrowprops=dict(arrowstyle="-", color=st.BLACK, lw=0.6),
        bbox=box,
    )
    ax.add_patch(Circle((0, 0), 0.08, color=st.RECEIVER, zorder=6))
    ax.text(-0.05, -0.35, r"$j$", ha="center")
    ax.text(
        -2.1, -2.2, "lobes: angular spread of each\nbeam after the hop cone", fontsize=8
    )
    ax.set_xlim(-2.2, 3.1)
    ax.set_ylim(-2.4, 3.0)
    ax.set_title("(b) two beams at $60^\\circ$")


def main() -> None:
    """Draw the figure."""
    st.setup()
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(st.TEXT_WIDTH, 2.75), gridspec_kw=dict(width_ratios=[1.25, 1])
    )
    panel_merge(a)
    panel_beams(b)
    for ax in (a, b):
        ax.set_aspect("equal")
        ax.set_axis_off()
    fig.subplots_adjust(wspace=0.05)
    st.save(fig, "isrf_merge.pdf")


if __name__ == "__main__":
    main()
