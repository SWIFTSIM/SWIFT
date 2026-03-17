"""Plot the cell geometry and particle distribution diagnostics.

This script generates diagnostic plots of the zoom region cell structure and
particle occupancy from files produced by running SWIFT with --zoom-geometry.

To run get the files to run this analysis run the following:

    swift --zoom-geometry <extra-options> <param file>

The following plots are produced:
  1. Cell grid coloured by cell type and subtype (geometry overview).
  2. Cell grid coloured by total particle occupancy (log-scaled heatmap).
  3. Cell grid coloured by MPI rank (domain decomposition).
  4. Particle count histograms split by cell type (zoom vs background).
  5. DM vs DM_background occupancy comparison for zoom cells.

Example:
    python plot_zoom_geometry.py --metadata zoom_metadata.yml \
                                --cells zoom_cell_data.dat
"""

import argparse
import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# ---------------------------------------------------------------------------
# Column index constants matching the C output format.
# ---------------------------------------------------------------------------
COL_CELL_ID = 0
COL_TYPE = 1
COL_SUBTYPE = 2
COL_RANK = 3
COL_LOC_X = 4
COL_LOC_Y = 5
COL_LOC_Z = 6
COL_WIDTH_X = 7
COL_WIDTH_Y = 8
COL_WIDTH_Z = 9
COL_MAXDEPTH = 10
COL_GAS = 11
COL_DM = 12
COL_DM_BKG = 13
COL_SINK = 14
COL_STARS = 15
COL_BH = 16
COL_NEUTRINO = 17

# Cell type enum values (must match cell.h)
CELL_TYPE_REGULAR = 0
CELL_TYPE_ZOOM = 1
CELL_TYPE_BKG = 2

# Cell subtype enum values (must match cell.h)
CELL_SUBTYPE_REGULAR = 0
CELL_SUBTYPE_NEIGHBOUR = 1
CELL_SUBTYPE_VOID = 2

# Colours for cell types / subtypes
COLOUR_BKG = "#bbe6e4"
COLOUR_BKG_NEIGHBOUR = "#f4a261"
COLOUR_BKG_VOID = "#e9c46a"
COLOUR_ZOOM = "#084b83"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_data(meta_file, data_file):
    """Load the metadata YAML and cell data table.

    Args:
        meta_file: Path to zoom_metadata.yml.
        data_file: Path to zoom_cell_data.dat.

    Returns:
        metadata: Dictionary of global zoom properties.
        data: Numpy array with one row per cell.
    """
    with open(meta_file, "r") as f:
        metadata = yaml.safe_load(f)

    data = np.loadtxt(data_file)

    return metadata, data


# ---------------------------------------------------------------------------
# Summary statistics (printed to stdout)
# ---------------------------------------------------------------------------
def print_summary(metadata, data):
    """Print a summary of the zoom geometry and particle counts to stdout.

    Args:
        metadata: Dictionary of global zoom properties.
        data: Numpy array with one row per cell.
    """
    # Separate zoom and background cells.
    zoom_mask = data[:, COL_TYPE] == CELL_TYPE_ZOOM
    bkg_mask = data[:, COL_TYPE] == CELL_TYPE_BKG

    # Total particle counts across all cells.
    total_gas = int(data[:, COL_GAS].sum())
    total_dm = int(data[:, COL_DM].sum())
    total_dm_bkg = int(data[:, COL_DM_BKG].sum())
    total_stars = int(data[:, COL_STARS].sum())
    total_bh = int(data[:, COL_BH].sum())
    total_sink = int(data[:, COL_SINK].sum())
    total_neutrino = int(data[:, COL_NEUTRINO].sum())
    total = (
        total_gas
        + total_dm
        + total_dm_bkg
        + total_stars
        + total_bh
        + total_sink
        + total_neutrino
    )

    # Per-cell-type totals.
    zoom_total = int(data[zoom_mask, COL_GAS : COL_NEUTRINO + 1].sum())
    bkg_total = int(data[bkg_mask, COL_GAS : COL_NEUTRINO + 1].sum())

    # Print.
    sep = "-" * 60
    print(sep)
    print("ZOOM GEOMETRY DIAGNOSTIC SUMMARY")
    print(sep)
    print(f"  Box size:            {metadata['BoxSize']}")
    print(f"  Background cdim:     {metadata['BackgroundCdim']}")
    print(f"  Zoom cdim:           {metadata['ZoomCdim']}")
    print(f"  Zoom region dim:     {metadata['ZoomRegionDim']}")
    print(f"  Zoom cell depth:     {metadata['ZoomCellDepth']}")
    print(f"  Nr zoom cells:       {metadata['NrZoomCells']}")
    print(f"  Nr bkg cells:        {metadata['NrBkgCells']}")
    print(f"  Total cells:         {metadata['TotalCells']}")
    print(sep)
    print("PARTICLE COUNTS")
    print(sep)
    print(f"  {'Type':<20} {'Total':>12} {'In Zoom':>12} {'In Bkg':>12}")
    print(f"  {'-' * 20} {'-' * 12} {'-' * 12} {'-' * 12}")

    for label, col in [
        ("Gas", COL_GAS),
        ("DM", COL_DM),
        ("DM Background", COL_DM_BKG),
        ("Sinks", COL_SINK),
        ("Stars", COL_STARS),
        ("Black Holes", COL_BH),
        ("Neutrinos", COL_NEUTRINO),
    ]:
        t = int(data[:, col].sum())
        z = int(data[zoom_mask, col].sum())
        b = int(data[bkg_mask, col].sum())
        if t > 0:
            print(f"  {label:<20} {t:>12d} {z:>12d} {b:>12d}")

    print(f"  {'-' * 20} {'-' * 12} {'-' * 12} {'-' * 12}")
    print(f"  {'All'::<20} {total:>12d} {zoom_total:>12d} {bkg_total:>12d}")
    print(sep)

    # Occupancy statistics for non-empty cells.
    for label, mask in [("Zoom", zoom_mask), ("Background", bkg_mask)]:
        counts = data[mask, COL_GAS : COL_NEUTRINO + 1].sum(axis=1)
        nonempty = counts[counts > 0]
        if len(nonempty) > 0:
            print(
                f"  {label} cells: {int(mask.sum())} total, "
                f"{len(nonempty)} non-empty, "
                f"occupancy min/median/max = "
                f"{int(nonempty.min())}/{int(np.median(nonempty))}/{int(nonempty.max())}"
            )
        else:
            print(f"  {label} cells: {int(mask.sum())} total, all empty")

    # MPI rank distribution.
    ranks = np.unique(data[:, COL_RANK].astype(int))
    if len(ranks) > 1:
        print(
            f"  MPI ranks present: {len(ranks)} (ranks {int(ranks.min())}-{int(ranks.max())})"
        )
    else:
        print("  MPI ranks present: 1 (serial run)")
    print(sep)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
def _get_cell_colour_type(ctype, csubtype):
    """Return face colour for a cell based on its type and subtype.

    Args:
        ctype: Cell type integer.
        csubtype: Cell subtype integer.

    Returns:
        Colour string and z-order.
    """
    if ctype == CELL_TYPE_ZOOM:
        return COLOUR_ZOOM, 2
    elif ctype == CELL_TYPE_BKG:
        if csubtype == CELL_SUBTYPE_VOID:
            return COLOUR_BKG_VOID, 1
        elif csubtype == CELL_SUBTYPE_NEIGHBOUR:
            return COLOUR_BKG_NEIGHBOUR, 1
        else:
            return COLOUR_BKG, 0
    else:
        return None, -1


def _draw_cells(ax, data, facecolors, edgecolor="black", linewidth=0.3):
    """Draw cells as rectangles on the given axes.

    Args:
        ax: Matplotlib axes.
        data: Cell data array.
        facecolors: Array or list of colours, one per cell.
        edgecolor: Edge colour for all rectangles.
        linewidth: Line width for cell edges.
    """
    for i, row in enumerate(data):
        loc_x, loc_y = row[COL_LOC_X], row[COL_LOC_Y]
        w_x, w_y = row[COL_WIDTH_X], row[COL_WIDTH_Y]
        fc = facecolors[i]
        if fc is None:
            continue
        ax.add_patch(
            Rectangle(
                (loc_x, loc_y),
                w_x,
                w_y,
                fill=True,
                facecolor=fc,
                edgecolor=edgecolor,
                linewidth=linewidth,
                zorder=1,
            )
        )


def _setup_axes(ax, metadata):
    """Configure axes limits, labels, and aspect ratio.

    Args:
        ax: Matplotlib axes.
        metadata: Zoom metadata dictionary.
    """
    box_size = metadata["BoxSize"][0]
    ax.set_xlim(0, box_size)
    ax.set_ylim(0, box_size)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")


# ---------------------------------------------------------------------------
# Individual plot functions
# ---------------------------------------------------------------------------
def plot_cell_types(metadata, data, outdir):
    """Plot 1: Cell grid coloured by type and subtype.

    Args:
        metadata: Zoom metadata dictionary.
        data: Cell data array.
        outdir: Output directory for saved figures.

    Returns:
        The matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(7, 7))
    _setup_axes(ax, metadata)
    ax.set_title("Cell Types & Subtypes")

    # Assign colours per cell based on type/subtype.
    colours = []
    for row in data:
        c, _ = _get_cell_colour_type(int(row[COL_TYPE]), int(row[COL_SUBTYPE]))
        colours.append(c)
    _draw_cells(ax, data, colours)

    # Legend
    handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=COLOUR_BKG,
            markeredgecolor="black",
            markersize=12,
            label="Background",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=COLOUR_BKG_NEIGHBOUR,
            markeredgecolor="black",
            markersize=12,
            label="Bkg (neighbour)",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=COLOUR_BKG_VOID,
            markeredgecolor="black",
            markersize=12,
            label="Bkg (void)",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=COLOUR_ZOOM,
            markeredgecolor="black",
            markersize=12,
            label="Zoom",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.06),
        ncol=len(handles),
        frameon=False,
        fontsize=9,
    )

    fig.tight_layout()
    path = os.path.join(outdir, "zoom_cell_types.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    return fig


def plot_occupancy(metadata, data, outdir):
    """Plot 2: Cell grid coloured by total particle occupancy (log scale).

    Args:
        metadata: Zoom metadata dictionary.
        data: Cell data array.
        outdir: Output directory for saved figures.

    Returns:
        The matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(7, 7))
    _setup_axes(ax, metadata)
    ax.set_title("Total Particle Occupancy (log$_{10}$)")

    # Total count per cell.
    total_counts = data[:, COL_GAS : COL_NEUTRINO + 1].sum(axis=1)

    # Log-scale the counts (use nan for empty cells).
    with np.errstate(divide="ignore"):
        log_counts = np.log10(total_counts)
    log_counts[~np.isfinite(log_counts)] = np.nan

    # Build a colourmap.
    valid = log_counts[np.isfinite(log_counts)]
    if len(valid) > 0:
        vmin, vmax = np.nanmin(valid), np.nanmax(valid)
    else:
        vmin, vmax = 0, 1
    cmap = plt.cm.viridis
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    colours = []
    for i in range(len(data)):
        if np.isfinite(log_counts[i]):
            colours.append(cmap(norm(log_counts[i])))
        else:
            colours.append("#dddddd")  # grey for empty cells
    _draw_cells(ax, data, colours)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("log$_{10}$(particle count)")

    fig.tight_layout()
    path = os.path.join(outdir, "zoom_occupancy.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    return fig


def plot_ranks(metadata, data, outdir):
    """Plot 3: Cell grid coloured by MPI rank.

    Args:
        metadata: Zoom metadata dictionary.
        data: Cell data array.
        outdir: Output directory for saved figures.

    Returns:
        The matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(7, 7))
    _setup_axes(ax, metadata)
    ax.set_title("MPI Domain Decomposition")

    ranks = data[:, COL_RANK].astype(int)
    unique_ranks = np.unique(ranks)
    n_ranks = len(unique_ranks)

    # Use a qualitative colourmap for ranks.
    if n_ranks <= 10:
        cmap = plt.cm.tab10
    elif n_ranks <= 20:
        cmap = plt.cm.tab20
    else:
        cmap = plt.cm.turbo
    norm = mcolors.Normalize(
        vmin=unique_ranks.min() - 0.5, vmax=unique_ranks.max() + 0.5
    )

    colours = [cmap(norm(r)) for r in ranks]
    _draw_cells(ax, data, colours)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("MPI Rank")
    if n_ranks <= 20:
        cbar.set_ticks(unique_ranks)

    fig.tight_layout()
    path = os.path.join(outdir, "zoom_ranks.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    return fig


def plot_occupancy_histograms(metadata, data, outdir):
    """Plot 4: Histograms of per-cell particle counts split by cell type.

    Produces one panel per particle type that has non-zero counts.

    Args:
        metadata: Zoom metadata dictionary.
        data: Cell data array.
        outdir: Output directory for saved figures.

    Returns:
        The matplotlib Figure, or None if no particles.
    """
    particle_types = [
        ("Gas", COL_GAS),
        ("DM", COL_DM),
        ("DM Background", COL_DM_BKG),
        ("Sinks", COL_SINK),
        ("Stars", COL_STARS),
        ("Black Holes", COL_BH),
        ("Neutrinos", COL_NEUTRINO),
    ]

    # Filter to types with non-zero counts.
    active = [(label, col) for label, col in particle_types if data[:, col].sum() > 0]
    if len(active) == 0:
        print("No particles found, skipping occupancy histograms.")
        return None

    n_panels = len(active)
    ncols = min(n_panels, 3)
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False
    )

    zoom_mask = data[:, COL_TYPE] == CELL_TYPE_ZOOM
    bkg_mask = data[:, COL_TYPE] == CELL_TYPE_BKG

    for idx, (label, col) in enumerate(active):
        ax = axes[idx // ncols, idx % ncols]

        zoom_counts = data[zoom_mask, col]
        bkg_counts = data[bkg_mask, col]

        # Only histogram non-zero cells.
        zoom_nz = zoom_counts[zoom_counts > 0]
        bkg_nz = bkg_counts[bkg_counts > 0]

        # Determine common bins in log space.
        all_nz = np.concatenate([zoom_nz, bkg_nz])
        if len(all_nz) == 0:
            ax.set_title(label)
            ax.text(
                0.5,
                0.5,
                "No particles",
                transform=ax.transAxes,
                ha="center",
                va="center",
            )
            continue

        lo = np.log10(all_nz.min())
        hi = np.log10(all_nz.max())
        if lo == hi:
            lo -= 0.5
            hi += 0.5
        bins = np.logspace(lo, hi, 30)

        if len(zoom_nz) > 0:
            ax.hist(
                zoom_nz,
                bins=bins,
                alpha=0.7,
                label="Zoom",
                color=COLOUR_ZOOM,
                edgecolor="black",
                linewidth=0.5,
            )
        if len(bkg_nz) > 0:
            ax.hist(
                bkg_nz,
                bins=bins,
                alpha=0.7,
                label="Background",
                color=COLOUR_BKG,
                edgecolor="black",
                linewidth=0.5,
            )

        ax.set_xscale("log")
        ax.set_xlabel("Particles per cell")
        ax.set_ylabel("Number of cells")
        ax.set_title(label)
        ax.legend(fontsize=8)

    # Hide unused panels.
    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)

    fig.suptitle("Cell Occupancy Distributions", fontsize=14, y=1.02)
    fig.tight_layout()
    path = os.path.join(outdir, "zoom_occupancy_histograms.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    return fig


def plot_dm_split(metadata, data, outdir):
    """Plot 5: DM vs DM_background occupancy for zoom cells.

    Shows a 2D scatter of normal DM count vs background DM count per zoom
    cell, which is useful for diagnosing particle contamination in the zoom
    region.

    Args:
        metadata: Zoom metadata dictionary.
        data: Cell data array.
        outdir: Output directory for saved figures.

    Returns:
        The matplotlib Figure, or None if no DM particles.
    """
    zoom_mask = data[:, COL_TYPE] == CELL_TYPE_ZOOM
    zoom_data = data[zoom_mask]

    dm = zoom_data[:, COL_DM]
    dm_bkg = zoom_data[:, COL_DM_BKG]

    if dm.sum() == 0 and dm_bkg.sum() == 0:
        print("No DM particles in zoom cells, skipping DM split plot.")
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left panel: scatter of DM vs DM_bkg per zoom cell.
    ax = axes[0]
    # Only show cells with at least one DM particle of either type.
    mask = (dm > 0) | (dm_bkg > 0)
    if mask.sum() > 0:
        sc = ax.scatter(
            dm[mask] + 1,
            dm_bkg[mask] + 1,
            c=dm_bkg[mask] / np.maximum(dm[mask] + dm_bkg[mask], 1),
            cmap="RdYlGn_r",
            s=8,
            alpha=0.7,
            edgecolors="none",
            vmin=0,
            vmax=1,
        )
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label("Background DM fraction")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Normal DM count + 1")
    ax.set_ylabel("Background DM count + 1")
    ax.set_title("DM Contamination (zoom cells)")

    # Right panel: spatial map of background DM fraction in zoom cells.
    ax = axes[1]
    box_size = metadata["BoxSize"][0]
    ax.set_xlim(0, box_size)
    ax.set_ylim(0, box_size)
    ax.set_aspect("equal")
    ax.set_title("Bkg DM Fraction (zoom cells)")

    total_dm = dm + dm_bkg
    frac = np.where(total_dm > 0, dm_bkg / total_dm, 0.0)
    cmap = plt.cm.RdYlGn_r
    norm = mcolors.Normalize(vmin=0, vmax=max(frac.max(), 0.01))

    for i in range(len(zoom_data)):
        row = zoom_data[i]
        fc = cmap(norm(frac[i])) if total_dm[i] > 0 else "#dddddd"
        ax.add_patch(
            Rectangle(
                (row[COL_LOC_X], row[COL_LOC_Y]),
                row[COL_WIDTH_X],
                row[COL_WIDTH_Y],
                fill=True,
                facecolor=fc,
                edgecolor="black",
                linewidth=0.2,
            )
        )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Background DM fraction")

    fig.tight_layout()
    path = os.path.join(outdir, "zoom_dm_split.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    return fig


def plot_per_type_maps(metadata, data, outdir):
    """Plot 6: Spatial maps of individual particle type counts.

    One panel per particle type with non-zero counts, showing the log-scaled
    occupancy on the cell grid.

    Args:
        metadata: Zoom metadata dictionary.
        data: Cell data array.
        outdir: Output directory for saved figures.

    Returns:
        The matplotlib Figure, or None if no particles.
    """
    particle_types = [
        ("Gas", COL_GAS),
        ("DM", COL_DM),
        ("DM Background", COL_DM_BKG),
        ("Sinks", COL_SINK),
        ("Stars", COL_STARS),
        ("Black Holes", COL_BH),
        ("Neutrinos", COL_NEUTRINO),
    ]

    active = [(label, col) for label, col in particle_types if data[:, col].sum() > 0]
    if len(active) == 0:
        return None

    n_panels = len(active)
    ncols = min(n_panels, 3)
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(6 * ncols, 5.5 * nrows), squeeze=False
    )

    for idx, (label, col) in enumerate(active):
        ax = axes[idx // ncols, idx % ncols]
        box_size = metadata["BoxSize"][0]
        ax.set_xlim(0, box_size)
        ax.set_ylim(0, box_size)
        ax.set_aspect("equal")
        ax.set_title(label)

        counts = data[:, col]
        with np.errstate(divide="ignore"):
            log_counts = np.log10(counts)
        log_counts[~np.isfinite(log_counts)] = np.nan

        valid = log_counts[np.isfinite(log_counts)]
        if len(valid) > 0:
            vmin, vmax = np.nanmin(valid), np.nanmax(valid)
        else:
            vmin, vmax = 0, 1
        cmap = plt.cm.inferno
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        for i, row in enumerate(data):
            if np.isfinite(log_counts[i]):
                fc = cmap(norm(log_counts[i]))
            else:
                fc = "#dddddd"
            ax.add_patch(
                Rectangle(
                    (row[COL_LOC_X], row[COL_LOC_Y]),
                    row[COL_WIDTH_X],
                    row[COL_WIDTH_Y],
                    fill=True,
                    facecolor=fc,
                    edgecolor="black",
                    linewidth=0.2,
                    zorder=1,
                )
            )

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(f"log$_{{10}}$({label} count)")

    # Hide unused panels.
    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)

    fig.suptitle("Per-Type Particle Occupancy", fontsize=14, y=1.02)
    fig.tight_layout()
    path = os.path.join(outdir, "zoom_per_type_maps.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    return fig


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Analyse zoom geometry diagnostic files produced by "
        "SWIFT --zoom-geometry.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Output files:\n"
            "  zoom_cell_types.png            Cell grid by type/subtype\n"
            "  zoom_occupancy.png             Total particle occupancy heatmap\n"
            "  zoom_ranks.png                 MPI domain decomposition\n"
            "  zoom_occupancy_histograms.png  Per-type occupancy histograms\n"
            "  zoom_dm_split.png              DM vs DM_bkg contamination\n"
            "  zoom_per_type_maps.png         Per-type spatial occupancy maps\n"
        ),
    )
    parser.add_argument(
        "--metadata",
        default="zoom_metadata.yml",
        help="Path to the metadata YAML file (default: zoom_metadata.yml)",
    )
    parser.add_argument(
        "--cells",
        default="zoom_cell_data.dat",
        help="Path to the cell data table (default: zoom_cell_data.dat)",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory for plots (default: current directory)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show plots interactively after saving",
    )
    args = parser.parse_args()

    # Create output directory if needed.
    if args.outdir != "." and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Load data.
    metadata, data = load_data(args.metadata, args.cells)

    # Print summary statistics to stdout.
    print_summary(metadata, data)

    # Generate all plots.
    figs = []
    figs.append(plot_cell_types(metadata, data, args.outdir))
    figs.append(plot_occupancy(metadata, data, args.outdir))
    figs.append(plot_ranks(metadata, data, args.outdir))
    figs.append(plot_occupancy_histograms(metadata, data, args.outdir))
    figs.append(plot_dm_split(metadata, data, args.outdir))
    figs.append(plot_per_type_maps(metadata, data, args.outdir))

    if args.show:
        plt.show()

    # Close all figures.
    for f in figs:
        if f is not None:
            plt.close(f)


if __name__ == "__main__":
    main()
