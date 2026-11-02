"""Geometry of one hop (path and receiving slot) and the reach of a sender."""

import numpy as np
import isrf_style as st
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Wedge


def arrow(ax, p, q, col, ls="-", lw=1.2, ms=9):
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


def dot(ax, p, col, label, dx=0.0, dy=0.0, **kw):
    """Draw a particle with its label."""
    ax.add_patch(Circle(p, 0.09, color=col, zorder=5, **kw))
    ax.text(p[0] + dx, p[1] + dy, label, ha="center", va="center")


def panel_path(ax) -> None:
    """Hop from sender i to receiver j: offset, path length and slot direction."""
    O = np.array([0.0, 0.0])
    n = np.array([1.0, 0.0])
    d_s = 3.0
    i = O + d_s * n
    th = np.radians(48)
    j = i + 2.0 * np.array([np.cos(th), np.sin(th)])
    d_r = np.linalg.norm(j - O)
    u = (j - O) / d_r
    perp = np.array([-u[1], u[0]])
    box = dict(fc="white", ec="none", pad=0.3)

    ax.add_patch(
        Arc(O, 2 * d_s, 2 * d_s, theta1=-10, theta2=40, color=st.GREY, ls="--", lw=0.8)
    )
    ax.plot(*np.array([O, i]).T, color=st.SENDER, lw=1.0)
    arrow(ax, i, i + 1.5 * n, st.SENDER, lw=1.0)
    ax.text(*(i + [1.6, 0.0]), r"$\boldsymbol{n}_k$", color=st.SENDER, va="center")
    ax.text(1.5, -0.3, r"$d_s$", color=st.SENDER, ha="center")
    arrow(ax, i, j, st.BLACK)
    ax.text(*(0.5 * (i + j) + [0.3, -0.15]), r"$\boldsymbol{\Delta x}_{ij}$")
    arrow(ax, j, O + 0.1 * u, st.RECEIVER, ls="--", lw=0.9)
    ax.plot(
        *np.array([d_s * u, j]).T,
        color=st.PATH,
        lw=3.2,
        solid_capstyle="butt",
        zorder=4,
    )
    mid = 0.5 * (d_s * u + j)
    ax.text(
        *(mid + 0.42 * perp),
        r"$d_r-d_s$",
        color=st.PATH,
        ha="center",
        va="center",
        bbox=box,
    )
    ax.text(
        *(0.4 * j + 0.38 * perp),
        r"$\boldsymbol{o}_j$, $|\boldsymbol{o}_j|=d_r$",
        color=st.RECEIVER,
        ha="center",
        va="center",
        rotation=np.degrees(np.arctan2(u[1], u[0])),
        bbox=box,
    )
    arrow(ax, j, j + 0.8 * u, st.RECEIVER, lw=1.0)
    ax.text(
        *(j + 0.9 * u + [-0.9, 0.3]),
        r"$-\boldsymbol{o}_j$ sets the slot",
        color=st.RECEIVER,
        ha="center",
    )
    ax.add_patch(Arc(i, 0.9, 0.9, theta1=0, theta2=48, color=st.BLACK, lw=0.8))
    ax.text(
        *(i + 0.62 * np.array([np.cos(np.radians(24)), np.sin(np.radians(24))])),
        r"$\theta$",
        ha="center",
        va="center",
    )
    ax.add_patch(Circle(O, 0.07, color=st.STAR, zorder=5))
    ax.text(0.0, -0.32, "origin", ha="center")
    dot(ax, i, st.SENDER, r"$i$", dx=0.0, dy=-0.3)
    dot(ax, j, st.RECEIVER, r"$j$", dx=0.28, dy=-0.05)
    ax.text(4.0, -0.4, r"$\mu_{ij,k}=\cos\theta$", ha="left")
    ax.set_xlim(-0.4, 6.2)
    ax.set_ylim(-0.6, 3.4)
    ax.set_title("(a) path and receiving slot of one share")


def panel_reach(ax) -> None:
    """Reach of sender i: its kernel, and the kernel of an active receiver."""
    i = np.array([0.0, 0.0])
    g_hi = 1.6
    box = dict(fc="white", ec="none", pad=0.3)
    n = np.array([1.0, 0.0])
    for t in range(-60, 61, 4):
        w = np.cos(np.radians(t)) ** 4
        ax.add_patch(
            Wedge(i, g_hi, t - 2, t + 2, color=st.SENDER, alpha=0.35 * w, lw=0)
        )
    ax.add_patch(Circle(i, g_hi, fill=False, color=st.SENDER, lw=1.0))
    ax.text(-2.15, 0.9, r"$\gamma h_i$", color=st.SENDER, bbox=box)
    arrow(ax, i, i + 0.8 * n, st.SENDER, lw=1.0)
    ax.text(0.45, 0.17, r"$\boldsymbol{n}_k$", color=st.SENDER)

    j1 = np.array([1.05, -0.55])
    j2 = np.array([2.25, 0.75])
    j3 = np.array([1.95, -1.15])
    ax.add_patch(Circle(j2, 2.45, fill=False, color=st.RECEIVER, lw=1.0, ls="--"))
    ax.text(
        *(j2 + [0.0, 2.5]),
        r"$\gamma h_{j}$",
        color=st.RECEIVER,
        ha="center",
        va="bottom",
        bbox=box,
    )
    dot(ax, i, st.SENDER, r"$i$", dx=-0.25, dy=-0.15)
    dot(ax, j1, st.RECEIVER, "", dx=0.0)
    ax.text(
        *(j1 + [-0.55, -0.75]),
        "inside $\\gamma h_i$:\nin reach",
        ha="center",
        fontsize=8,
        bbox=box,
    )
    dot(ax, j2, st.RECEIVER, "")
    ax.text(
        *(j2 + [0.2, 0.35]),
        "active $j$, $r_{ij}<\\gamma h_j$:\nin reach",
        ha="left",
        fontsize=8,
        bbox=box,
    )
    ax.add_patch(
        Circle(j3, 0.09, facecolor="white", edgecolor=st.RECEIVER, lw=1.0, zorder=5)
    )
    ax.text(
        *(j3 + [0.2, -0.35]),
        "inactive $j$, beyond $\\gamma h_i$:\nno light",
        ha="left",
        fontsize=8,
        bbox=box,
    )
    ax.text(
        -2.4,
        -2.35,
        r"weight: $V_j\,W(r_{ij},H_{ij})\,\mu^p$, $H_{ij}=\max(h_i,h_j)$",
        ha="left",
        fontsize=8,
    )
    ax.set_xlim(-2.5, 4.9)
    ax.set_ylim(-2.6, 3.65)
    ax.set_title("(b) reach of sender $i$")


def main() -> None:
    """Draw the figure."""
    st.setup()
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(st.TEXT_WIDTH, 3.0), gridspec_kw=dict(width_ratios=[1.12, 1])
    )
    panel_path(a)
    panel_reach(b)
    for ax in (a, b):
        ax.set_aspect("equal")
        ax.set_axis_off()
    fig.subplots_adjust(wspace=0.05)
    st.save(fig, "isrf_hop.pdf")


if __name__ == "__main__":
    main()
