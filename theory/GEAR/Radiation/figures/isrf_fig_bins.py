"""Arrival bins and far cells of one receiver slot under a finite speed of light."""

import numpy as np
import isrf_style as st
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

N_BINS, N_FAR, FAR_SCALE = 4, 2, 2


def main() -> None:
    """Draw the figure."""
    st.setup()
    fig, (top, bot) = plt.subplots(
        2,
        1,
        figsize=(st.TEXT_WIDTH, 3.6),
        sharex=True,
        gridspec_kw=dict(height_ratios=[1.25, 1], hspace=0.32),
    )
    width = 2**FAR_SCALE
    far_end0 = (N_BINS // width + 1) * width

    # Lattice of the receiver: near bins one step wide, far cells 2^n steps wide.
    for j in range(N_BINS):
        top.add_patch(
            Rectangle(
                (j, 0), 1, 0.5, facecolor=st.SKY, alpha=0.35, edgecolor=st.BLACK, lw=0.8
            )
        )
        top.text(j + 0.5, 0.25, f"bin {j}", ha="center", va="center", fontsize=8)
    for c in range(N_FAR):
        x0 = far_end0 - width + c * width
        top.add_patch(
            Rectangle(
                (x0, 0),
                width,
                0.5,
                facecolor=st.ORANGE,
                alpha=0.3,
                edgecolor=st.BLACK,
                lw=0.8,
            )
        )
        top.text(
            x0 + width / 2,
            0.25,
            f"far cell {c}: $2^n$ steps, $n={FAR_SCALE}$",
            ha="center",
            va="center",
            fontsize=8,
        )

    # A share from a sender with a step twice as long, spread evenly over that step.
    t_a, dt_s = 0.4, 2.0
    top.add_patch(
        Rectangle((t_a, 0.75), dt_s, 0.32, facecolor=st.SENDER, alpha=0.75, lw=0)
    )
    top.text(
        0.0,
        1.13,
        "share of a sender with the step $\\Delta t_s=2\\Delta t$,\nspread evenly over "
        "$[t_a,\\,t_a+\\Delta t_s]$",
        va="bottom",
        fontsize=8,
    )
    for j in range(N_BINS):
        part = max(0.0, min(j + 1, t_a + dt_s) - max(j, t_a)) / dt_s
        if part > 0:
            top.text(
                j + 0.5,
                0.58,
                f"{part:.1f}",
                ha="center",
                va="center",
                fontsize=7,
                color=st.SENDER,
            )

    # A far share, arriving from t_first to t_last.
    t_f, t_l = 5.0, 10.5
    top.add_patch(
        Rectangle((t_f, 0.75), t_l - t_f, 0.32, facecolor=st.VERMILION, alpha=0.6, lw=0)
    )
    top.text(
        0.5 * (t_f + t_l),
        0.91,
        "far share, from $t_{\\rm first}$ to $t_{\\rm last}$",
        ha="center",
        va="center",
        fontsize=8,
    )
    top.set_ylim(0, 1.6)
    top.set_yticks([])
    for side in ("top", "right", "left"):
        top.spines[side].set_visible(False)
    top.set_title(
        "(a) where the waiting light of one slot is kept, at the merge $t$ (now)"
    )

    # The released fraction of the far share at each merge never exceeds what has arrived.
    t = np.linspace(0, 12, 400)
    arrived = np.clip((t - t_f) / (t_l - t_f), 0, 1)
    bot.plot(
        t, arrived, color=st.VERMILION, label="light of the far share that has arrived"
    )
    merges = np.arange(1, 13)
    released = np.clip((merges - t_f) / (t_l - t_f), 0, 1)
    bot.step(
        np.r_[0, merges],
        np.r_[0, released],
        where="post",
        color=st.BLACK,
        lw=1.0,
        label="released at the merges (one per step)",
    )
    bot.plot(merges, released, "o", color=st.BLACK, ms=3)
    bot.set_ylabel("fraction")
    bot.set_ylim(-0.05, 1.12)
    bot.legend(loc="upper left", frameon=False)
    for side in ("top", "right"):
        bot.spines[side].set_visible(False)
    bot.set_xlim(-0.2, 12.4)
    bot.set_xticks(range(0, 13))
    bot.set_xticklabels(
        ["$t$"]
        + [f"$t+{k}\\Delta t$" if k in (1, 4, 8, 12) else "" for k in range(1, 13)]
    )
    bot.set_xlabel("time of the receiver, in steps $\\Delta t$ of the receiver")
    bot.set_title("(b) release of the far share", pad=2)
    for ax in (top, bot):
        for k in range(0, 13):
            ax.axvline(k, color=st.GREY, lw=0.3, zorder=0)
    st.save(fig, "isrf_bins.pdf")


if __name__ == "__main__":
    main()
