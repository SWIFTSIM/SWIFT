"""Life of one share: from the slot of a sender to the held packet of a receiver."""

import isrf_style as st
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


def box(ax, x, y, w, h, text, col):
    """Draw a rounded box centred on (x, y) with a text."""
    ax.add_patch(
        FancyBboxPatch(
            (x - w / 2, y - h / 2),
            w,
            h,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            fc="white",
            ec=col,
            lw=1.1,
        )
    )
    ax.text(x, y, text, ha="center", va="center", fontsize=7.5, linespacing=1.2)


def arrow(ax, p, q, col=st.BLACK, ls="-"):
    """Draw an arrow from p to q."""
    ax.add_patch(
        FancyArrowPatch(
            p,
            q,
            arrowstyle="-|>",
            mutation_scale=8,
            color=col,
            lw=0.9,
            ls=ls,
            shrinkA=0,
            shrinkB=0,
        )
    )


def main() -> None:
    """Draw the figure."""
    st.setup()
    fig, ax = plt.subplots(figsize=(st.TEXT_WIDTH, 2.55))
    ax.set_xlim(0, 12.2)
    ax.set_ylim(0, 4.3)
    ax.axis("off")

    w, h, y = 2.15, 1.05, 2.9
    xs = [1.15, 3.55, 5.95, 8.35, 10.95]
    steps = [
        ("sender $i$, slot $k$\nshare $\\phi_{ij,k}$\nof $E$ and $N$", st.SENDER),
        ("path\n$\\Delta\\ell=d_r-d_s$\n(star: $r_{\\star j}$)", st.PATH),
        ("redshift\nif $H>0$;\nLW feeds PE", st.LW),
        ("dust: keep\n$e^{-\\tau_b}$, with\n$\\kappa_b$ of $j$", st.PE),
        (
            "slot of $j$ set\nby $-\\boldsymbol{o}_j$; add\n$e$ and $e\\,\\boldsymbol{o}_j$",
            st.RECEIVER,
        ),
    ]
    for x, (t, c) in zip(xs, steps):
        box(ax, x, y, w if x < 10 else 2.35, h, t, c)
    for a, b in zip(xs[:-1], xs[1:]):
        arrow(ax, (a + w / 2 + 0.03, y), (b - (w if b < 10 else 2.35) / 2 - 0.03, y))

    # Ledgers fed on the way.
    ly = 0.95
    box(ax, xs[2], ly, 1.9, 0.62, "ledger $Y$", st.GREY)
    box(ax, xs[3], ly, 1.9, 0.62, "ledger $A$\n(and cutoffs)", st.GREY)
    arrow(ax, (xs[2], y - h / 2 - 0.03), (xs[2], ly + 0.34), st.GREY, "--")
    arrow(ax, (xs[3], y - h / 2 - 0.03), (xs[3], ly + 0.34), st.GREY, "--")

    # The merge.
    box(
        ax,
        xs[4],
        ly,
        2.35,
        1.0,
        "merge at a step\nof $j$: sum $E$, $N$;\nmean origin",
        st.BLACK,
    )
    arrow(ax, (xs[4], y - h / 2 - 0.03), (xs[4], ly + 0.53))
    ax.text(
        xs[4] - 0.08,
        (y - h / 2 + ly + 0.5) / 2,
        "when due",
        fontsize=7.5,
        ha="right",
        va="center",
        color=st.GREY,
    )
    ax.text(
        xs[0] - w / 2,
        4.05,
        "booked at the hop of $i$",
        fontsize=8,
        color=st.GREY,
        ha="left",
    )
    ax.text(
        xs[4] + 1.18,
        3.95,
        "waits: inbox, reservoir,\nnear bin or far cell",
        fontsize=8,
        color=st.GREY,
        ha="right",
    )
    st.save(fig, "isrf_share.pdf")


if __name__ == "__main__":
    main()
