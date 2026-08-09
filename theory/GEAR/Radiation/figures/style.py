"""Shared plotting style for the GEAR HII neighbour-search explainer figures.

Centralises the colorblind-safe palette, font sizes, and small drawing
helpers (boxes, arrows, axis ticks) reused across fig1..fig6 so every
figure reads as one consistent set.
"""

from __future__ import annotations

from typing import Sequence

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import FancyArrowPatch

# Okabe-Ito colorblind-safe palette.
BLACK = "#000000"
ORANGE = "#E69F00"
SKY_BLUE = "#56B4E9"
BLUISH_GREEN = "#009E73"
YELLOW = "#F0E442"
BLUE = "#0072B2"
VERMILLION = "#D55E00"
REDDISH_PURPLE = "#CC79A7"
GRAY = "#999999"

STAR_COLOR = VERMILLION
CANDIDATE_COLOR = SKY_BLUE
IONIZED_COLOR = BLUISH_GREEN
WINDOW_COLOR = ORANGE
WIRED_COLOR = BLUE
UNWIRED_COLOR = GRAY
CELL_EDGE = BLACK

FS_TITLE = 11
FS_LABEL = 9
FS_ANNOT = 8
FS_SMALL = 7.5


def new_figure(figsize: tuple[float, float]) -> tuple[plt.Figure, Axes]:
    """Create a figure/axes pair with a consistent single-column style.

    Parameters
    ----------
    figsize : tuple of float
        Figure size in inches, (width, height).

    Returns
    -------
    tuple
        The (figure, axes) pair.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    return fig, ax


def save_figure(fig: plt.Figure, path_no_ext: str) -> None:
    """Save a figure as both PNG and PDF next to each other.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    path_no_ext : str
        Output path without extension; ``.png`` and ``.pdf`` are appended.
    """
    fig.savefig(f"{path_no_ext}.png", dpi=220, bbox_inches="tight")
    fig.savefig(f"{path_no_ext}.pdf", bbox_inches="tight")
    plt.close(fig)


def draw_box(
    ax: Axes,
    xy: tuple[float, float],
    w: float,
    h: float,
    *,
    edgecolor: str = CELL_EDGE,
    facecolor: str = "none",
    lw: float = 1.3,
    ls: str = "-",
    alpha: float = 1.0,
    zorder: int = 2,
) -> mpatches.Rectangle:
    """Draw an axis-aligned cell box and add it to the axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    xy : tuple of float
        Lower-left corner (x, y).
    w, h : float
        Width and height.

    Returns
    -------
    matplotlib.patches.Rectangle
        The patch added to the axes.
    """
    rect = mpatches.Rectangle(
        xy,
        w,
        h,
        edgecolor=edgecolor,
        facecolor=facecolor,
        lw=lw,
        ls=ls,
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_patch(rect)
    return rect


def draw_arrow(
    ax: Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = BLACK,
    lw: float = 1.3,
    style: str = "-|>",
    connectionstyle: str = "arc3,rad=0.0",
    ls: str = "-",
    zorder: int = 3,
) -> FancyArrowPatch:
    """Draw a directional arrow between two points.

    Returns
    -------
    matplotlib.patches.FancyArrowPatch
        The patch added to the axes.
    """
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle=style,
        mutation_scale=11,
        color=color,
        lw=lw,
        linestyle=ls,
        connectionstyle=connectionstyle,
        zorder=zorder,
    )
    ax.add_patch(arrow)
    return arrow


def star_marker(ax: Axes, xy: tuple[float, float], size: float = 220, zorder: int = 6):
    """Plot the star (spart) marker used consistently in every figure."""
    ax.scatter(
        [xy[0]],
        [xy[1]],
        marker="*",
        s=size,
        color=STAR_COLOR,
        edgecolor=BLACK,
        linewidth=0.6,
        zorder=zorder,
    )


def legend_swatch(color: str, label: str, marker: str = "s") -> mpatches.Patch:
    """Build a small legend handle (kept separate so callers control order)."""
    if marker == "s":
        return mpatches.Patch(facecolor=color, edgecolor=BLACK, label=label)
    raise NotImplementedError
