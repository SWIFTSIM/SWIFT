"""Overview of the packet transport in two dimensions: emission, then one hop."""

import numpy as np
import isrf_style as st
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Wedge

GAMMA = 1.936492
P = 4


def wendland(q: np.ndarray) -> np.ndarray:
    """Return the Wendland C2 shape (1-q)^4 (1+4q) for q < 1."""
    q = np.clip(q, 0, 1)
    return (1 - q) ** 4 * (1 + 4 * q)


def arrow(ax, p, q, col, lw=1.0, alpha=1.0, ms=7):
    """Draw an arrow from p to q."""
    ax.add_patch(
        FancyArrowPatch(
            p,
            q,
            arrowstyle="-|>",
            mutation_scale=ms,
            color=col,
            lw=lw,
            alpha=alpha,
            shrinkA=4,
            shrinkB=4,
        )
    )


def particles() -> np.ndarray:
    """Return a jittered lattice of gas positions."""
    rng = np.random.default_rng(7)
    xs, ys = np.meshgrid(np.arange(-1.2, 6.01, 0.9), np.arange(-2.25, 2.3, 0.9))
    pos = np.c_[xs.ravel(), ys.ravel()] + rng.uniform(-0.22, 0.22, (xs.size, 2))
    return pos


def main() -> None:
    """Draw the figure."""
    st.setup()
    fig, ax = plt.subplots(figsize=(st.TEXT_WIDTH, 3.0))
    box = dict(fc="white", ec="none", pad=0.3)
    gas = particles()
    star = np.array([0.45, 0.05])
    h_star = 0.75
    r_kernel = 0.95 * GAMMA * h_star

    # Emission cells of the star: 12 equal cells in a random frame (2D schematic of HEALPix pixels).
    rot = 17.0
    for c in range(12):
        ax.add_patch(
            Wedge(
                star,
                r_kernel,
                rot + 30 * c,
                rot + 30 * (c + 1),
                facecolor=st.STAR if c % 2 else "none",
                alpha=0.18,
                edgecolor=st.STAR,
                lw=0.5,
            )
        )
    ax.add_patch(Circle(star, r_kernel, fill=False, color=st.STAR, lw=0.9))

    # One hop: sender i sends its slot down a cone of weights V W mu^p.
    i = star + np.array([1.2, 0.15])
    gas = np.r_[gas[np.linalg.norm(gas - i, axis=1) > 0.55], [i]]
    n_k = (i - star) / np.linalg.norm(i - star)
    h_i = 1.0
    dx = gas - i
    r = np.linalg.norm(dx, axis=1)
    mu = np.where(r > 0, (dx @ n_k) / np.where(r > 0, r, 1), 0)
    w = np.where(
        (r > 0) & (r < GAMMA * h_i) & (mu > 1e-3),
        wendland(r / (GAMMA * h_i)) * np.clip(mu, 0, 1) ** P,
        0,
    )
    phi = w / w.sum()
    ax.add_patch(Circle(i, GAMMA * h_i, fill=False, color=st.SENDER, lw=0.8, ls="--"))
    for g, f in zip(gas, phi):
        if f > 0.003:
            d = g - i
            d = d * (1 - 0.12 / np.linalg.norm(d))
            wd = 0.02 + 0.16 * f / phi.max()
            ax.arrow(
                *i,
                *d,
                width=wd,
                head_width=wd + 0.12,
                head_length=0.15,
                length_includes_head=True,
                color=st.SENDER,
                alpha=0.4 + 0.6 * f / phi.max(),
                lw=0,
                zorder=4,
            )
    near = np.linalg.norm(gas - star, axis=1) < r_kernel
    for g in gas[near]:
        arrow(ax, star, g, st.STAR, lw=0.9)
    j_idx = int(np.argmax(phi))
    j = gas[j_idx]

    ax.scatter(*gas.T, s=16, color=st.GREY, zorder=4)
    ax.scatter(*i, s=30, color=st.SENDER, zorder=5)
    ax.scatter(*j, s=30, color=st.RECEIVER, zorder=5)
    ax.scatter(
        *star, s=120, marker="*", color=st.STAR, edgecolor=st.BLACK, lw=0.5, zorder=6
    )

    # The receiving slot: the direction from the origin (the star) to j.
    u = (j - star) / np.linalg.norm(j - star)
    ax.plot(
        *np.array([star, j + 1.4 * u]).T, color=st.RECEIVER, lw=0.8, ls=":", zorder=3
    )
    ax.text(
        *(i + [0.0, 0.32]), r"$i$", color=st.SENDER, ha="center", bbox=box, zorder=7
    )
    ax.text(
        *(j + [0.0, 0.32]), r"$j$", color=st.RECEIVER, ha="center", bbox=box, zorder=7
    )
    ax.text(
        *(j + 1.45 * u + [0.1, 0.0]),
        "slot of $j$: direction\nfrom the origin (star) to $j$",
        color=st.RECEIVER,
        ha="left",
        va="center",
        fontsize=8,
        bbox=box,
        zorder=7,
    )
    ax.text(
        star[0] - 0.2,
        star[1] - r_kernel - 0.35,
        "star: 12 emission cells\nin its kernel",
        ha="center",
        color=st.BLACK,
        fontsize=8,
        bbox=box,
        zorder=7,
    )
    ax.text(
        i[0] + 2.3,
        i[1] - GAMMA * h_i - 0.2,
        "hop of sender $i$: weight\n$\\propto V_jW\\mu^p$ (arrow width)",
        color=st.SENDER,
        ha="center",
        va="top",
        fontsize=8,
        bbox=box,
        zorder=7,
    )
    ax.set_aspect("equal")
    ax.set_xlim(-1.6, 6.9)
    ax.set_ylim(-3.2, 2.6)
    ax.set_axis_off()
    st.save(fig, "isrf_overview.pdf")


if __name__ == "__main__":
    main()
