"""Shared style of the ISRF transport figures: fonts, sizes and palette."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito palette, safe for the common forms of colour blindness.
BLACK = "#000000"
ORANGE = "#E69F00"
SKY = "#56B4E9"
GREEN = "#009E73"
YELLOW = "#F0E442"
BLUE = "#0072B2"
VERMILION = "#D55E00"
PURPLE = "#CC79A7"
GREY = "#8C8C8C"
LIGHT = "#D9D9D9"

# Roles, the same in every figure.
STAR = ORANGE
SENDER = BLUE
RECEIVER = VERMILION
PE = SKY
LW = PURPLE
PATH = GREEN

# Width of the text block of the document, inches.
TEXT_WIDTH = 6.1


def setup() -> None:
    """Set the matplotlib defaults of the figures."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif"],
            "mathtext.fontset": "cm",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.7,
            "lines.linewidth": 1.2,
            "pdf.fonttype": 42,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )


def save(fig: plt.Figure, path: str) -> None:
    """Write a figure as a vector PDF and close it."""
    fig.savefig(path)
    plt.close(fig)
