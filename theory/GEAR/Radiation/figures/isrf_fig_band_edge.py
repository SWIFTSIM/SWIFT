"""Sketch of the band-edge transfer of the redshift between the PE and LW bands."""

import numpy as np
import isrf_style as st
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

E_PE, E_LW, E_HI = 6.0, 11.2, 13.6
LAMBDA_E_PE, LAMBDA_E_LW = 2.154, 6.508


def spectrum(E: np.ndarray) -> np.ndarray:
    """Return a spectrum I_E whose band-edge weights are the shipped defaults.

    A flat PE band gives lambda_E(PE) - 1 = E_lo / (E_hi - E_lo) = 1.154. A
    linear LW band from 1 at 11.2 eV to i_hi at 13.6 eV gives
    lambda_E(LW) - 1 = 2 E_lo / ((1 + i_hi)(E_hi - E_lo)) = 5.508.
    """
    i_hi = 2.0 * E_LW / ((LAMBDA_E_LW - 1.0) * (E_HI - E_LW)) - 1.0
    lw = 1.0 + (i_hi - 1.0) * (E - E_LW) / (E_HI - E_LW)
    return np.where(E < E_LW, 1.0, lw)


def main() -> None:
    """Draw the figure."""
    st.setup()
    fig, ax = plt.subplots(figsize=(st.TEXT_WIDTH, 2.5))
    E = np.linspace(E_PE, E_HI, 400)
    I = spectrum(E)
    ax.fill_between(E[E <= E_LW], 0, I[E <= E_LW], color=st.PE, alpha=0.45, lw=0)
    ax.fill_between(E[E >= E_LW], 0, I[E >= E_LW], color=st.LW, alpha=0.45, lw=0)
    ax.plot(E, I, color=st.BLACK)
    ax.plot([3.6, E_PE], [1.0, 1.0], color=st.BLACK, ls=":")
    ax.plot([E_HI, 15.8], [I[-1], 1.2], color=st.BLACK, ls=":")

    ax.text(8.6, 0.16, "PE band\n$I_{\\rm PE}=\\int_{\\rm PE} I_E\\,dE$", ha="center")
    ax.text(12.4, 0.3, "LW band\n$I_{\\rm LW}$", ha="center")
    ax.text(
        4.6,
        0.45,
        "below 6 eV:\nnot transported",
        ha="center",
        color=st.GREY,
        fontsize=8,
    )
    ax.text(
        15.0, 0.45, "ionising:\nHII channel", ha="center", color=st.GREY, fontsize=8
    )

    # The strip E_lo I_E(E_lo): the energy that leaves through a lower edge.
    for E_lo, col in [(E_PE, st.PE), (E_LW, st.LW)]:
        ax.add_patch(
            Rectangle(
                (E_lo, 0),
                0.18,
                spectrum(np.array([E_lo + 0.01]))[0],
                facecolor="none",
                edgecolor=st.BLACK,
                hatch="////",
                lw=0.6,
            )
        )

    y = 1.32
    arrows = [
        (E_HI + 0.5, E_HI - 0.5, "dropped", st.GREY, "-|>"),
        (E_LW + 0.5, E_LW - 0.5, "LW $\\to$ PE", st.BLACK, "-|>"),
        (E_PE + 0.5, E_PE - 0.5, "lost", st.BLACK, "-|>"),
    ]
    for x0, x1, label, col, style in arrows:
        ax.add_patch(
            FancyArrowPatch(
                (x0, y),
                (x1, y),
                arrowstyle=style,
                mutation_scale=9,
                color=col,
                lw=1.1,
                ls="--" if label == "dropped" else "-",
            )
        )
        ax.text(0.5 * (x0 + x1), y + 0.07, label, ha="center", color=col)
    ax.text(
        8.6,
        0.76,
        "$E_{{\\rm lo},b}\\,I_E(E_{{\\rm lo},b})$:\nflux out through\nthe lower edge",
        ha="center",
        va="center",
    )
    ax.plot([7.3, E_PE + 0.09], [0.76, 0.62], color=st.BLACK, lw=0.6)
    ax.plot([9.9, E_LW + 0.09], [0.76, 0.62], color=st.BLACK, lw=0.6)
    for E_edge in (E_PE, E_LW, E_HI):
        ax.axvline(E_edge, color=st.GREY, lw=0.6, ls="-", zorder=0)
    ax.set_xlim(3.4, 16.4)
    ax.set_ylim(0, 1.66)
    ax.set_xticks([E_PE, E_LW, E_HI])
    ax.set_xticklabels(["6.0", "11.2", "13.6"])
    ax.set_yticks([])
    ax.set_xlabel("photon energy $E$ [eV]")
    ax.set_ylabel("$I_E$ (arbitrary units)")
    ax.text(
        9.9,
        1.53,
        "the redshift moves every photon to lower $E$; what crosses an edge:",
        ha="center",
        fontsize=8,
    )
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    st.save(fig, "isrf_band_edge.pdf")


if __name__ == "__main__":
    main()
