"""Plot the cell geometry and particle distribution diagnostics.

This script generates diagnostic plots of the zoom region cell structure and
particle occupancy from files produced by running SWIFT with --zoom-geometry.

To run get the files to run this analysis run the following:

    swift --zoom-geometry <extra-options> <param file>

The following plots are produced:
  1. Cell grid coloured by cell type and subtype (geometry overview).
  2. Particle count histograms split by cell type (zoom vs background).
  3. DM vs DM_background occupancy scatter (contamination diagnostic).
  4. Padding analysis: particle extent vs zoom region with breakdown.
  5. Resolution quality: occupancy uniformity, radial profiles, and
     DM contamination as a function of distance from the zoom centre.

Example:
    python plot_zoom_geometry.py

    # Or with custom file paths and output directory:
    python plot_zoom_geometry.py --metadata my_zoom_metadata.yml \
        --cells my_zoom_cell_data.dat --outdir zoom_plots
"""

import argparse
import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# Column index constants matching the C output format.
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


def _fmt_vec(v, fmt=".4f"):
    """Format a 3-vector as a compact string."""
    return f"({v[0]:{fmt}}, {v[1]:{fmt}}, {v[2]:{fmt}})"


def _fmt_ivec(v):
    """Format an integer 3-vector as a compact string."""
    return f"({v[0]}, {v[1]}, {v[2]})"


def _compute_neighbour_shells(data, void_mask, neighbour_mask):
    """Compute the number of neighbour cells in each shell around void cells.

    A shell is defined as all neighbour cells at a given Manhattan distance
    from the void region (measured in background cell units).

    Args:
        data: Cell data array.
        void_mask: Boolean mask for void cells.
        neighbour_mask: Boolean mask for neighbour cells.

    Returns:
        List of integers giving the number of neighbour cells in each shell,
        ordered from closest to farthest from the void region.
    """
    if void_mask.sum() == 0 or neighbour_mask.sum() == 0:
        return []

    # Get void cell indices and their centres.
    void_data = data[void_mask]
    void_centres = np.column_stack(
        [
            void_data[:, COL_LOC_X] + void_data[:, COL_WIDTH_X] / 2,
            void_data[:, COL_LOC_Y] + void_data[:, COL_WIDTH_Y] / 2,
            void_data[:, COL_LOC_Z] + void_data[:, COL_WIDTH_Z] / 2,
        ]
    )

    # Get neighbour cell centres.
    neighbour_data = data[neighbour_mask]
    neighbour_centres = np.column_stack(
        [
            neighbour_data[:, COL_LOC_X] + neighbour_data[:, COL_WIDTH_X] / 2,
            neighbour_data[:, COL_LOC_Y] + neighbour_data[:, COL_WIDTH_Y] / 2,
            neighbour_data[:, COL_LOC_Z] + neighbour_data[:, COL_WIDTH_Z] / 2,
        ]
    )

    # For each neighbour cell, compute the minimum distance to any void cell
    # (in units of background cell widths).
    bkg_width = neighbour_data[0, COL_WIDTH_X]  # Assume cubic cells.
    min_distances = np.full(len(neighbour_centres), np.inf)

    for nc_idx, nc in enumerate(neighbour_centres):
        # Compute distances to all void centres.
        diffs = void_centres - nc[np.newaxis, :]
        dists = np.sqrt(np.sum(diffs**2, axis=1))
        min_distances[nc_idx] = dists.min() / bkg_width

    # Bin by distance (rounded to nearest integer = shell number).
    # Shell 1 = cells touching void, shell 2 = one cell away, etc.
    shell_indices = np.round(min_distances).astype(int)
    unique_shells = np.unique(shell_indices)
    unique_shells = unique_shells[
        unique_shells > 0
    ]  # Exclude shell 0 (shouldn't happen).

    shell_counts = []
    for shell in sorted(unique_shells):
        count = (shell_indices == shell).sum()
        shell_counts.append(count)

    return shell_counts


def _simulate_setup(box_size, part_dim, user_pad, depth, bkg_cdim):
    """Simulate zoom geometry computation for a given bkg_cdim."""
    # 1. Compute particle extent and requested padding
    ini_dim = np.max(part_dim)
    requested_dim = ini_dim * user_pad

    # 2. Simulate void snapping logic
    bkg_width = box_size / bkg_cdim

    # Void snapping
    lower = np.floor(((box_size - requested_dim) / 2.0) / bkg_width) * bkg_width
    upper = (np.floor(((box_size + requested_dim) / 2.0) / bkg_width) + 1) * bkg_width
    actual_dim = upper - lower

    actual_pad = actual_dim / ini_dim

    # 3. Compute cell counts
    # Number of zoom regions tessellating void region
    nr_zoom_regions = int(np.ceil(actual_dim / requested_dim))

    # This is a simplification based on zoom_get_geometry
    # zoom_cdim = nr_parents * 2^depth
    # nr_parents = actual_dim / bkg_width
    cdim = int(np.floor((actual_dim + (0.1 * bkg_width)) / bkg_width) * (2**depth))

    nr_zoom_cells = cdim**3
    nr_bkg_cells = bkg_cdim**3

    return actual_pad, nr_bkg_cells + nr_zoom_cells


def _find_optimal_bkg_cdim(metadata):
    """Find the optimal bkg_cdim to minimize padding waste and total cells."""
    if "ParticleDim" not in metadata:
        return None

    box_size = metadata["BoxSize"][0]
    part_dim = np.array(metadata["ParticleDim"])
    user_pad = metadata["UserRegionPadFactor"]
    depth = metadata["ZoomCellDepth"]
    current_bkg = metadata["BackgroundCdim"][0]

    # Sweep range
    possible_cdims = np.arange(max(8, current_bkg // 2), current_bkg * 2 + 1)

    results = []
    for cdim in possible_cdims:
        actual_pad, total_cells = _simulate_setup(
            box_size, part_dim, user_pad, depth, cdim
        )

        # Cost function: normalize padding ratio (want ~1.0) and cell count
        # Padding waste ratio
        pad_waste = actual_pad / user_pad

        # Simple weighted cost function
        cost = pad_waste + (total_cells / 1e5)

        results.append((cdim, pad_waste, total_cells, cost))

    # Return best
    return min(results, key=lambda x: x[3])


def print_summary(metadata, data):
    """Print a summary of the zoom geometry and particle counts to stdout.

    Args:
        metadata: Dictionary of global zoom properties.
        data: Numpy array with one row per cell.
    """
    # Separate zoom and background cells.
    zoom_mask = data[:, COL_TYPE] == CELL_TYPE_ZOOM
    bkg_mask = data[:, COL_TYPE] == CELL_TYPE_BKG

    # Sub-type masks for background cells.
    void_mask = (data[:, COL_TYPE] == CELL_TYPE_BKG) & (
        data[:, COL_SUBTYPE] == CELL_SUBTYPE_VOID
    )
    neighbour_mask = (data[:, COL_TYPE] == CELL_TYPE_BKG) & (
        data[:, COL_SUBTYPE] == CELL_SUBTYPE_NEIGHBOUR
    )
    regular_bkg_mask = (data[:, COL_TYPE] == CELL_TYPE_BKG) & (
        data[:, COL_SUBTYPE] == CELL_SUBTYPE_REGULAR
    )

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

    # Formatting constants - compute dynamically based on content.
    # We'll measure sample lines to determine the needed width.
    LBL = 36  # label column width

    # Build sample lines to measure maximum width needed.
    sample_lines = []

    # Cell grid samples
    sample_lines.append(f"  {'Box size':<{LBL}} {_fmt_vec(metadata['BoxSize'])}")
    sample_lines.append(
        f"  {'Background cdim':<{LBL}} {_fmt_ivec(metadata['BackgroundCdim'])}"
    )
    sample_lines.append(
        f"  {'Background cell width':<{LBL}} {_fmt_vec(metadata['BackgroundCellWidth'])}"
    )

    # Zoom region samples
    zoom_dim = np.array(metadata["ZoomRegionDim"])
    sample_lines.append(f"  {'Zoom region dim':<{LBL}} {_fmt_vec(zoom_dim)}")

    # Cell counts (with potential large numbers)
    nr_zoom = metadata["NrZoomCells"]
    nr_bkg = metadata["NrBkgCells"]
    nr_void = metadata.get("NrVoidCells", int(void_mask.sum()))
    nr_neighbour = metadata.get("NrNeighbourCells", int(neighbour_mask.sum()))
    nr_regular_bkg = int(regular_bkg_mask.sum())
    sample_lines.append(f"  {'Nr background cells':<{LBL}} {nr_bkg}")

    # Particle counts table (different format)
    sample_lines.append(f"  {'Type':<18} {'Total':>12} {'In Zoom':>12} {'In Bkg':>12}")
    sample_lines.append(
        f"  {'DM Background':<18} {total_dm_bkg:>12,d} {total_dm_bkg:>12,d} {total_dm_bkg:>12,d}"
    )
    sample_lines.append(
        f"  {'All':<18} {total:>12,d} {zoom_total:>12,d} {bkg_total:>12,d}"
    )

    # Occupancy statistics (typically the longest)
    zoom_counts = data[zoom_mask, COL_GAS : COL_NEUTRINO + 1].sum(axis=1)
    zoom_nonempty = zoom_counts[zoom_counts > 0]
    if len(zoom_nonempty) > 0:
        sample_lines.append(
            f"  {'Zoom':<16} {int(zoom_mask.sum()):>6} cells "
            f"({len(zoom_nonempty):>6} non-empty)  "
            f"min {int(zoom_nonempty.min()):>8,d}  "
            f"med {int(np.median(zoom_nonempty)):>8,d}  "
            f"max {int(zoom_nonempty.max()):>8,d}"
        )
    bkg_counts = data[bkg_mask, COL_GAS : COL_NEUTRINO + 1].sum(axis=1)
    bkg_nonempty = bkg_counts[bkg_counts > 0]
    if len(bkg_nonempty) > 0:
        sample_lines.append(
            f"  {'Background':<16} {int(bkg_mask.sum()):>6} cells "
            f"({len(bkg_nonempty):>6} non-empty)  "
            f"min {int(bkg_nonempty.min()):>8,d}  "
            f"med {int(np.median(bkg_nonempty)):>8,d}  "
            f"max {int(bkg_nonempty.max()):>8,d}"
        )

    # Resolution quality samples
    if len(zoom_nonempty) > 1:
        mean_occ = np.mean(zoom_nonempty)
        std_occ = np.std(zoom_nonempty)
        cv = std_occ / mean_occ if mean_occ > 0 else 0.0
        sample_lines.append(
            f"  {'Coeff. of variation':<{LBL}} {cv:.3f}  (0 = perfectly uniform)"
        )
        ratio = zoom_nonempty.max() / zoom_nonempty.min()
        sample_lines.append(
            f"  {'** WARNING **':<{LBL}} Max/min occupancy ratio = {ratio:,.0f}"
        )
        sample_lines.append(
            f"  {'':<{LBL}} Large spread may indicate poor zoom region centering."
        )

    # Compute the maximum width needed
    max_width = max(len(line) for line in sample_lines)
    # Add a small buffer to be safe
    W = max_width + 2

    sep = "=" * W
    thin = "-" * W

    # ---- Header ----
    print()
    print(sep)
    print("  ZOOM GEOMETRY DIAGNOSTIC SUMMARY".center(W))
    print(sep)

    # ---- Cell grid ----
    print()
    print("  Cell Grid".ljust(W))
    print(f"  {thin[2:]}")
    print(f"  {'Box size':<{LBL}} {_fmt_vec(metadata['BoxSize'])}")
    print(f"  {'Background cdim':<{LBL}} {_fmt_ivec(metadata['BackgroundCdim'])}")
    print(
        f"  {'Background cell width':<{LBL}} "
        f"{_fmt_vec(metadata['BackgroundCellWidth'])}"
    )
    print(f"  {'Zoom cdim':<{LBL}} {_fmt_ivec(metadata['ZoomCdim'])}")
    print(f"  {'Zoom cell width':<{LBL}} {_fmt_vec(metadata['ZoomCellWidth'])}")
    print(f"  {'Zoom cell depth':<{LBL}} {metadata['ZoomCellDepth']}")
    print()

    # Compute neighbor shells around void cells.
    neighbour_shells = _compute_neighbour_shells(data, void_mask, neighbour_mask)

    print(f"  {'Nr zoom cells':<{LBL}} {nr_zoom}")
    print(f"  {'Nr background cells':<{LBL}} {nr_bkg}")
    print(f"    {'Void':<{LBL - 4}} {nr_void}")
    print(f"    {'Neighbour':<{LBL - 4}} {nr_neighbour}")
    if len(neighbour_shells) > 0:
        shells_str = ", ".join([f"{n}" for n in neighbour_shells])
        print(
            f"      {'Neighbour shells':<{LBL - 8}} {len(neighbour_shells)} ({shells_str} cells/shell)"
        )
    print(f"    {'Regular':<{LBL - 4}} {nr_regular_bkg}")
    print(f"  {'Total cells':<{LBL}} {metadata['TotalCells']}")

    # ---- Zoom region geometry & padding ----
    print()
    print("  Zoom Region Geometry".ljust(W))
    print(f"  {thin[2:]}")

    zoom_dim = np.array(metadata["ZoomRegionDim"])
    zoom_lower = np.array(metadata["ZoomRegionLowerBounds"])
    zoom_upper = np.array(metadata["ZoomRegionUpperBounds"])
    zoom_com = np.array(metadata["ZoomRegionCoM"])
    applied_shift = np.array(metadata["AppliedZoomShift"])
    actual_pad = metadata["RegionPadFactor"]

    print(f"  {'Zoom region dim':<{LBL}} {_fmt_vec(zoom_dim)}")
    print(f"  {'Zoom region lower':<{LBL}} {_fmt_vec(zoom_lower)}")
    print(f"  {'Zoom region upper':<{LBL}} {_fmt_vec(zoom_upper)}")
    print(f"  {'Zoom region CoM':<{LBL}} {_fmt_vec(zoom_com)}")
    print(f"  {'Applied shift':<{LBL}} {_fmt_vec(applied_shift)}")

    # Padding analysis (if metadata available).
    has_padding_info = "ParticleDim" in metadata
    if has_padding_info:
        part_dim = np.array(metadata["ParticleDim"])
        user_pad = metadata["UserRegionPadFactor"]
        max_part_dim = np.max(part_dim)

        # Padding breakdown: particle extent -> user-padded -> actual (cell-aligned)
        user_padded_dim = max_part_dim * user_pad
        actual_dim = np.max(zoom_dim)

        # Padding shell thickness on each side.
        padding_each_side = (zoom_dim - part_dim) / 2.0
        padding_frac = padding_each_side / (zoom_dim / 2.0) * 100.0

        print()
        print("  Padding Analysis".ljust(W))
        print(f"  {thin[2:]}")
        print(f"  {'Particle extent':<{LBL}} {_fmt_vec(part_dim)}")
        print(f"  {'Max particle extent':<{LBL}} {max_part_dim:.4f}")
        print(f"  {'Requested pad factor':<{LBL}} {user_pad:.4f}")
        print(f"  {'Padded extent (requested)':<{LBL}} {user_padded_dim:.4f}")
        print(f"  {'Final region extent':<{LBL}} {actual_dim:.4f}")
        print(f"  {'Actual pad factor':<{LBL}} {actual_pad:.4f}")
        print(f"  {'Pad factor ratio':<{LBL}} {actual_pad / user_pad:.2f}x requested")
        print(f"  {'Padding shell (per side)':<{LBL}} {_fmt_vec(padding_each_side)}")
        print(f"  {'Padding fraction':<{LBL}} {_fmt_vec(padding_frac, fmt='.1f')}%")

        # Volume efficiency: what fraction of zoom volume is particles vs padding.
        part_vol = np.prod(part_dim)
        zoom_vol = np.prod(zoom_dim)
        vol_eff = part_vol / zoom_vol * 100.0
        print(f"  {'Volume efficiency':<{LBL}} {vol_eff:.1f}%")

    # ---- Particle counts ----
    print()
    print("  Particle Counts".ljust(W))
    print(f"  {thin[2:]}")
    hdr = f"  {'Type':<18} {'Total':>12} {'In Zoom':>12} {'In Bkg':>12}"
    print(hdr)
    print(f"  {'-' * 18} {'-' * 12} {'-' * 12} {'-' * 12}")

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
            print(f"  {label:<18} {t:>12,d} {z:>12,d} {b:>12,d}")

    print(f"  {'-' * 18} {'-' * 12} {'-' * 12} {'-' * 12}")
    print(f"  {'All':<18} {total:>12,d} {zoom_total:>12,d} {bkg_total:>12,d}")

    # ---- DM contamination summary ----
    zoom_dm_total = int(data[zoom_mask, COL_DM].sum())
    zoom_dm_bkg_total = int(data[zoom_mask, COL_DM_BKG].sum())
    if zoom_dm_total + zoom_dm_bkg_total > 0:
        contam_frac = zoom_dm_bkg_total / (zoom_dm_total + zoom_dm_bkg_total)
        n_contaminated = int(((data[zoom_mask, COL_DM_BKG]) > 0).sum())
        print()
        print("  DM Contamination (zoom cells)".ljust(W))
        print(f"  {thin[2:]}")
        print(
            f"  {'Bkg DM in zoom region':<{LBL}} "
            f"{zoom_dm_bkg_total:,d} / {zoom_dm_total + zoom_dm_bkg_total:,d}"
        )
        print(f"  {'Contamination fraction':<{LBL}} {contam_frac:.6f}")
        print(
            f"  {'Contaminated cells':<{LBL}} {n_contaminated} / {int(zoom_mask.sum())}"
        )

    # ---- Occupancy statistics ----
    print()
    print("  Occupancy Statistics".ljust(W))
    print(f"  {thin[2:]}")
    for label, mask in [
        ("Zoom", zoom_mask),
        ("Background", bkg_mask),
        ("  Void", void_mask),
        ("  Neighbour", neighbour_mask),
    ]:
        counts = data[mask, COL_GAS : COL_NEUTRINO + 1].sum(axis=1)
        nonempty = counts[counts > 0]
        n_total = int(mask.sum())
        if n_total == 0:
            continue
        if len(nonempty) > 0:
            print(
                f"  {label:<16} {n_total:>6} cells "
                f"({len(nonempty):>6} non-empty)  "
                f"min {int(nonempty.min()):>8,d}  "
                f"med {int(np.median(nonempty)):>8,d}  "
                f"max {int(nonempty.max()):>8,d}"
            )
        else:
            print(f"  {label:<16} {n_total:>6} cells (all empty)")

    # ---- Resolution quality for zoom cells ----
    zoom_counts = data[zoom_mask, COL_GAS : COL_NEUTRINO + 1].sum(axis=1)
    zoom_nonempty = zoom_counts[zoom_counts > 0]
    if len(zoom_nonempty) > 1:
        mean_occ = np.mean(zoom_nonempty)
        std_occ = np.std(zoom_nonempty)
        cv = std_occ / mean_occ if mean_occ > 0 else 0.0
        iqr = np.percentile(zoom_nonempty, 75) - np.percentile(zoom_nonempty, 25)
        print()
        print("  Resolution Quality (zoom cells)".ljust(W))
        print(f"  {thin[2:]}")
        print(f"  {'Mean occupancy':<{LBL}} {mean_occ:,.1f}")
        print(f"  {'Std deviation':<{LBL}} {std_occ:,.1f}")
        print(f"  {'Coeff. of variation':<{LBL}} {cv:.3f}  (0 = perfectly uniform)")
        print(f"  {'IQR (25th-75th pctl)':<{LBL}} {iqr:,.1f}")
        print(f"  {'10th percentile':<{LBL}} {np.percentile(zoom_nonempty, 10):,.1f}")
        print(f"  {'90th percentile':<{LBL}} {np.percentile(zoom_nonempty, 90):,.1f}")

        # Flag extreme outliers.
        ratio = zoom_nonempty.max() / zoom_nonempty.min()
        if ratio > 100:
            print(f"  {'** WARNING **':<{LBL}} Max/min occupancy ratio = {ratio:,.0f}")
            print(
                f"  {'':<{LBL}} Large spread may indicate poor zoom region centering."
            )

    # ---- MPI ----
    ranks = np.unique(data[:, COL_RANK].astype(int))
    print()
    if len(ranks) > 1:
        zoom_ranks = np.unique(data[zoom_mask, COL_RANK].astype(int))
        bkg_ranks = np.unique(data[bkg_mask, COL_RANK].astype(int))
        print(
            f"  {'MPI ranks':<{LBL}} {len(ranks)} "
            f"(ranks {int(ranks.min())}-{int(ranks.max())})"
        )
        print(f"  {'Ranks with zoom cells':<{LBL}} {len(zoom_ranks)}")
        print(f"  {'Ranks with bkg cells':<{LBL}} {len(bkg_ranks)}")
    else:
        print(f"  {'MPI ranks':<{LBL}} 1 (serial run)")

    # ---- Setup Optimization ----
    best_bkg = _find_optimal_bkg_cdim(metadata)
    if best_bkg is not None:
        best_cdim, best_pad, best_total, _ = best_bkg
        current_bkg = metadata["BackgroundCdim"][0]
        if best_cdim != current_bkg:
            print()
            print(sep)
            print("  *** OPTIMAL ZOOM SETUP RECOMMENDATION ***".center(W))
            print(sep)
            print()
            print(f"  Your current setup:")
            print(f"    bkg_top_level_cells:     {current_bkg}")
            print(f"    Padding ratio:           {metadata['RegionPadFactor']:.4f}")
            print()
            print(f"  Recommended setup:")
            print(f"    bkg_top_level_cells:     {best_cdim}")
            print(f"    Expected padding ratio:  {best_pad:.4f}")
            print()
            print(f"  This reduces padding waste while minimizing total cell count.")
            print(sep)

    print()
    print(sep)


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

    # Draw cells with proper z-ordering (background first, then zoom on top).
    # Get colors and z-orders for each cell.
    for i, row in enumerate(data):
        c, zorder = _get_cell_colour_type(int(row[COL_TYPE]), int(row[COL_SUBTYPE]))
        if c is None:
            continue
        ax.add_patch(
            Rectangle(
                (row[COL_LOC_X], row[COL_LOC_Y]),
                row[COL_WIDTH_X],
                row[COL_WIDTH_Y],
                fill=True,
                facecolor=c,
                edgecolor="black",
                linewidth=0.3,
                zorder=zorder,
            )
        )

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

        # Add gridlines
        ax.grid(True, which="major", axis="both", alpha=0.3, zorder=0)
        ax.set_axisbelow(True)

        ax.set_xscale("log")
        ax.set_xlabel("Particles per cell")
        ax.set_ylabel("Number of cells")
        ax.set_title(label)
        ax.legend(fontsize=8, framealpha=0.9)

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
    """Plot: DM vs DM_background occupancy for zoom cells.

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

    fig, ax = plt.subplots(figsize=(7, 6))

    # Scatter of DM vs DM_bkg per zoom cell.
    # Only show cells with at least one DM particle of either type.
    mask = (dm > 0) | (dm_bkg > 0)
    if mask.sum() > 0:
        sc = ax.scatter(
            dm[mask] + 1,
            dm_bkg[mask] + 1,
            c=dm_bkg[mask] / np.maximum(dm[mask] + dm_bkg[mask], 1),
            cmap="RdYlGn_r",
            s=12,
            alpha=0.7,
            edgecolors="none",
            vmin=0,
            vmax=1,
        )
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Background DM fraction")

    # Add gridlines
    ax.grid(True, which="major", axis="both", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Normal DM count + 1")
    ax.set_ylabel("Background DM count + 1")
    ax.set_title("DM Contamination (zoom cells)")

    fig.tight_layout()
    path = os.path.join(outdir, "zoom_dm_split.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    return fig


def plot_padding_analysis(metadata, data, outdir):
    """Plot 7: Padding analysis showing particle extent vs zoom region.

    Left panel: 2D cell map with concentric rectangles showing the particle
    extent, user-requested padded extent, and final (cell-aligned) zoom region,
    zoomed in to the zoom region with a small margin.

    Right panel: 1D bar chart decomposing the zoom region extent into particle
    extent, user padding, and cell-alignment overhead.

    Args:
        metadata: Zoom metadata dictionary.
        data: Cell data array.
        outdir: Output directory for saved figures.

    Returns:
        The matplotlib Figure, or None if padding metadata is unavailable.
    """
    if "ParticleDim" not in metadata:
        print(
            "Padding metadata not available (run with newer SWIFT), "
            "skipping padding analysis plot."
        )
        return None

    part_dim = np.array(metadata["ParticleDim"])
    zoom_dim = np.array(metadata["ZoomRegionDim"])
    zoom_lower = np.array(metadata["ZoomRegionLowerBounds"])
    zoom_upper = np.array(metadata["ZoomRegionUpperBounds"])
    zoom_com = np.array(metadata["ZoomRegionCoM"])
    user_pad = metadata["UserRegionPadFactor"]
    box_size = np.array(metadata["BoxSize"])

    max_part_dim = np.max(part_dim)
    user_padded_dim = max_part_dim * user_pad

    # Compute centres and extents for concentric rectangles.
    zoom_centre = (zoom_lower + zoom_upper) / 2.0

    # Particle extent rectangle (centred on CoM).
    part_lower = zoom_com - part_dim / 2.0
    part_upper = zoom_com + part_dim / 2.0

    # User-padded rectangle (centred on CoM, isotropic from max extent).
    user_pad_lower = zoom_com - user_padded_dim / 2.0
    user_pad_upper = zoom_com + user_padded_dim / 2.0

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # ---- Left panel: spatial overview ----
    ax = axes[0]

    # Draw zoom cells (faded).
    zoom_mask = data[:, COL_TYPE] == CELL_TYPE_ZOOM
    zoom_data = data[zoom_mask]
    for row in zoom_data:
        ax.add_patch(
            Rectangle(
                (row[COL_LOC_X], row[COL_LOC_Y]),
                row[COL_WIDTH_X],
                row[COL_WIDTH_Y],
                fill=True,
                facecolor=COLOUR_ZOOM,
                edgecolor="black",
                linewidth=0.15,
                alpha=0.2,
                zorder=0,
            )
        )

    # Draw void cells.
    void_mask = (data[:, COL_TYPE] == CELL_TYPE_BKG) & (
        data[:, COL_SUBTYPE] == CELL_SUBTYPE_VOID
    )
    for row in data[void_mask]:
        ax.add_patch(
            Rectangle(
                (row[COL_LOC_X], row[COL_LOC_Y]),
                row[COL_WIDTH_X],
                row[COL_WIDTH_Y],
                fill=True,
                facecolor=COLOUR_BKG_VOID,
                edgecolor="black",
                linewidth=0.3,
                alpha=0.3,
                zorder=0,
            )
        )

    # Particle extent.
    ax.add_patch(
        Rectangle(
            (part_lower[0], part_lower[1]),
            part_dim[0],
            part_dim[1],
            fill=False,
            edgecolor="#e63946",
            linewidth=2.0,
            linestyle="-",
            label="Particle extent",
            zorder=3,
        )
    )

    # User-requested padded extent.
    ax.add_patch(
        Rectangle(
            (user_pad_lower[0], user_pad_lower[1]),
            user_padded_dim,
            user_padded_dim,
            fill=False,
            edgecolor="#457b9d",
            linewidth=2.0,
            linestyle="--",
            label=f"User padding ({user_pad:.2f}x)",
            zorder=3,
        )
    )

    # Final zoom region.
    ax.add_patch(
        Rectangle(
            (zoom_lower[0], zoom_lower[1]),
            zoom_dim[0],
            zoom_dim[1],
            fill=False,
            edgecolor="#2a9d8f",
            linewidth=2.0,
            linestyle="-.",
            label=f"Final region ({metadata['RegionPadFactor']:.2f}x)",
            zorder=3,
        )
    )

    # Mark CoM.
    ax.plot(
        zoom_com[0],
        zoom_com[1],
        "x",
        color="#e63946",
        markersize=10,
        markeredgewidth=2,
        zorder=4,
        label="Zoom CoM",
    )

    # Set view to zoom region with 30% margin.
    margin = 0.3 * np.max(zoom_dim)
    ax.set_xlim(zoom_lower[0] - margin, zoom_upper[0] + margin)
    ax.set_ylim(zoom_lower[1] - margin, zoom_upper[1] + margin)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Padding Overview (x-y projection)")
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)

    # ---- Right panel: 1D breakdown bar chart ----
    ax = axes[1]

    # Break the full extent into three concentric shells per axis.
    labels = ["x", "y", "z"]
    particle_extents = part_dim
    user_pads = np.full(3, user_padded_dim) - part_dim
    alignment_pads = zoom_dim - np.full(3, user_padded_dim)
    # Clamp negative alignment pads (can happen if void snapping shrinks).
    alignment_pads = np.maximum(alignment_pads, 0.0)

    x_pos = np.arange(3)
    width = 0.5

    ax.bar(
        x_pos,
        particle_extents,
        width,
        label="Particle extent",
        color="#e63946",
        edgecolor="black",
        linewidth=0.5,
    )
    ax.bar(
        x_pos,
        user_pads,
        width,
        bottom=particle_extents,
        label="User padding",
        color="#457b9d",
        edgecolor="black",
        linewidth=0.5,
    )
    ax.bar(
        x_pos,
        alignment_pads,
        width,
        bottom=particle_extents + user_pads,
        label="Cell alignment",
        color="#2a9d8f",
        edgecolor="black",
        linewidth=0.5,
    )

    # Add gridlines
    ax.grid(True, which="major", axis="y", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Extent")
    ax.set_title("Zoom Region Extent Breakdown")
    ax.legend(fontsize=9, loc="upper left", framealpha=0.9)

    # Annotate percentages.
    for i in range(3):
        total = zoom_dim[i]
        parts = [particle_extents[i], user_pads[i], alignment_pads[i]]
        bottom = 0
        for p in parts:
            if p > 0.01 * total:
                pct = p / total * 100
                ax.text(
                    x_pos[i],
                    bottom + p / 2,
                    f"{pct:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=8,
                    fontweight="bold",
                    color="white",
                )
            bottom += p

    fig.suptitle("Zoom Region Padding Analysis", fontsize=14, y=1.02)
    fig.tight_layout()
    path = os.path.join(outdir, "zoom_padding_analysis.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    return fig


def plot_resolution_quality(metadata, data, outdir):
    """Plot 8: Resolution quality analysis of the zoom region.

    Four panels:
      1. Occupancy distribution with percentile markers and uniformity metrics.
      2. Spatial map of occupancy deviation from the median.
      3. Radial profile of occupancy from the zoom region centre.
      4. Radial profile of DM background contamination fraction.

    Args:
        metadata: Zoom metadata dictionary.
        data: Cell data array.
        outdir: Output directory for saved figures.

    Returns:
        The matplotlib Figure, or None if no zoom cells.
    """
    zoom_mask = data[:, COL_TYPE] == CELL_TYPE_ZOOM
    zoom_data = data[zoom_mask]

    if len(zoom_data) == 0:
        print("No zoom cells found, skipping resolution quality plot.")
        return None

    # Total particle count per zoom cell.
    zoom_counts = zoom_data[:, COL_GAS : COL_NEUTRINO + 1].sum(axis=1)
    zoom_nonempty_mask = zoom_counts > 0

    if zoom_nonempty_mask.sum() < 2:
        print("Too few non-empty zoom cells, skipping resolution quality plot.")
        return None

    zoom_nonempty = zoom_counts[zoom_nonempty_mask]

    # Cell centres and radial distances from zoom region centre.
    cell_centres_x = zoom_data[:, COL_LOC_X] + zoom_data[:, COL_WIDTH_X] / 2.0
    cell_centres_y = zoom_data[:, COL_LOC_Y] + zoom_data[:, COL_WIDTH_Y] / 2.0
    cell_centres_z = zoom_data[:, COL_LOC_Z] + zoom_data[:, COL_WIDTH_Z] / 2.0

    zoom_lower = np.array(metadata["ZoomRegionLowerBounds"])
    zoom_dim = np.array(metadata["ZoomRegionDim"])
    zoom_centre = zoom_lower + zoom_dim / 2.0

    # Radial distance from zoom centre (3D).
    dx = cell_centres_x - zoom_centre[0]
    dy = cell_centres_y - zoom_centre[1]
    dz = cell_centres_z - zoom_centre[2]
    radii = np.sqrt(dx**2 + dy**2 + dz**2)

    # Normalised radius (0 = centre, 1 = corner of zoom region).
    max_radius = np.sqrt(np.sum((zoom_dim / 2.0) ** 2))
    norm_radii = radii / max_radius

    fig, axes = plt.subplots(2, 2, figsize=(13, 11))

    # ---- Panel 1: Occupancy distribution ----
    ax = axes[0, 0]
    valid_counts = zoom_nonempty
    with np.errstate(divide="ignore"):
        log_counts = np.log10(valid_counts)

    lo, hi = log_counts.min(), log_counts.max()
    if lo == hi:
        lo -= 0.5
        hi += 0.5
    bins = np.logspace(lo, hi, 40)
    ax.hist(
        valid_counts,
        bins=bins,
        color=COLOUR_ZOOM,
        edgecolor="black",
        linewidth=0.4,
        alpha=0.8,
    )

    # Percentile lines.
    for pct, ls, lbl in [
        (10, ":", "10th"),
        (50, "-", "Median"),
        (90, ":", "90th"),
    ]:
        val = np.percentile(valid_counts, pct)
        ax.axvline(val, color="#e63946", linestyle=ls, linewidth=1.5, label=lbl)

    # Annotate statistics.
    mean_occ = np.mean(valid_counts)
    cv = np.std(valid_counts) / mean_occ if mean_occ > 0 else 0
    ax.axvline(mean_occ, color="#2a9d8f", linestyle="--", linewidth=1.5, label="Mean")

    # Add gridlines
    ax.grid(True, which="major", axis="both", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)

    ax.set_xscale("log")
    ax.set_xlabel("Particles per cell")
    ax.set_ylabel("Number of cells")
    ax.set_title("Occupancy Distribution (zoom cells)")

    # Extend y-axis by 10% to make room for text box
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], ylim[1] * 1.1)

    ax.text(
        0.98,
        0.95,
        f"CV = {cv:.3f}\n"
        f"Mean = {mean_occ:,.0f}\n"
        f"N(empty) = {int((~zoom_nonempty_mask).sum())}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9),
        zorder=100,
    )

    ax.legend(fontsize=7, loc="upper left", framealpha=0.9)

    # ---- Panel 2: Deviation from median histogram ----
    ax = axes[0, 1]
    median_count = np.median(zoom_nonempty)

    # Log-ratio of count to median (positive = over-occupied).
    with np.errstate(divide="ignore", invalid="ignore"):
        log_ratio = np.log10(zoom_counts / median_count)
    log_ratio[~np.isfinite(log_ratio)] = np.nan

    valid_ratio = log_ratio[np.isfinite(log_ratio)]
    if len(valid_ratio) > 1:
        # Create histogram of deviations
        ax.hist(
            valid_ratio,
            bins=30,
            color=COLOUR_ZOOM,
            edgecolor="black",
            linewidth=0.4,
            alpha=0.8,
        )

        # Mark zero (median) line
        ax.axvline(
            0, color="#e63946", linestyle="-", linewidth=2, label="Median", zorder=10
        )

        # Mark mean deviation
        mean_dev = np.mean(valid_ratio)
        ax.axvline(
            mean_dev,
            color="#2a9d8f",
            linestyle="--",
            linewidth=1.5,
            label="Mean",
            zorder=10,
        )

        # Add gridlines
        ax.grid(True, which="major", axis="both", alpha=0.3, zorder=0)
        ax.set_axisbelow(True)

        ax.set_xlabel("log$_{10}$(count / median)")
        ax.set_ylabel("Number of cells")
        ax.set_title("Occupancy Deviation Distribution (zoom cells)")
        ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

        # Add text box with statistics
        std_dev = np.std(valid_ratio)
        ax.text(
            0.02,
            0.98,
            f"Std = {std_dev:.3f}\nMean = {mean_dev:.3f}\nMedian = 0.0",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9),
            zorder=100,
        )
    else:
        ax.text(
            0.5,
            0.5,
            "Insufficient data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
        )
        ax.set_title("Occupancy Deviation Distribution (zoom cells)")

    # ---- Panel 3: Radial occupancy profile ----
    ax = axes[1, 0]

    # Bin by normalised radius.
    n_bins = 20
    r_edges = np.linspace(0, 1, n_bins + 1)
    r_centres = (r_edges[:-1] + r_edges[1:]) / 2.0
    bin_medians = np.full(n_bins, np.nan)
    bin_q10 = np.full(n_bins, np.nan)
    bin_q90 = np.full(n_bins, np.nan)

    for i in range(n_bins):
        in_bin = (norm_radii >= r_edges[i]) & (norm_radii < r_edges[i + 1])
        bin_counts = zoom_counts[in_bin]
        bin_nz = bin_counts[bin_counts > 0]
        if len(bin_nz) > 0:
            bin_medians[i] = np.median(bin_nz)
            bin_q10[i] = np.percentile(bin_nz, 10)
            bin_q90[i] = np.percentile(bin_nz, 90)

    valid = np.isfinite(bin_medians)
    ax.fill_between(
        r_centres[valid],
        bin_q10[valid],
        bin_q90[valid],
        alpha=0.25,
        color=COLOUR_ZOOM,
        label="10th-90th pctl",
    )
    ax.plot(
        r_centres[valid],
        bin_medians[valid],
        "o-",
        color=COLOUR_ZOOM,
        markersize=4,
        linewidth=1.5,
        label="Median",
    )

    # Add gridlines
    ax.grid(True, which="major", axis="both", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)

    ax.set_xlabel("Normalised radius from zoom centre")
    ax.set_ylabel("Particles per cell")
    ax.set_yscale("log")
    ax.set_xlim(0, 1)
    ax.set_title("Radial Occupancy Profile (zoom cells)")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    # ---- Panel 4: Radial contamination profile ----
    ax = axes[1, 1]

    dm_counts = zoom_data[:, COL_DM]
    dm_bkg_counts = zoom_data[:, COL_DM_BKG]
    total_dm = dm_counts + dm_bkg_counts

    has_dm = total_dm.sum() > 0

    if has_dm:
        with np.errstate(divide="ignore", invalid="ignore"):
            contam_frac = np.where(total_dm > 0, dm_bkg_counts / total_dm, 0.0)

        bin_contam_median = np.full(n_bins, np.nan)
        bin_contam_max = np.full(n_bins, np.nan)
        bin_contam_mean = np.full(n_bins, np.nan)

        for i in range(n_bins):
            in_bin = (norm_radii >= r_edges[i]) & (norm_radii < r_edges[i + 1])
            bin_dm = total_dm[in_bin]
            bin_frac = contam_frac[in_bin]
            has_particles = bin_dm > 0
            if has_particles.sum() > 0:
                bin_contam_median[i] = np.median(bin_frac[has_particles])
                bin_contam_max[i] = np.max(bin_frac[has_particles])
                bin_contam_mean[i] = np.mean(bin_frac[has_particles])

        valid = np.isfinite(bin_contam_median)
        ax.plot(
            r_centres[valid],
            bin_contam_median[valid],
            "o-",
            color=COLOUR_ZOOM,
            markersize=4,
            linewidth=1.5,
            label="Median",
        )
        ax.plot(
            r_centres[valid],
            bin_contam_max[valid],
            "^--",
            color="#e63946",
            markersize=4,
            linewidth=1.0,
            label="Max",
        )
        ax.plot(
            r_centres[valid],
            bin_contam_mean[valid],
            "s:",
            color="#2a9d8f",
            markersize=3,
            linewidth=1.0,
            label="Mean",
        )

        # Add gridlines
        ax.grid(True, which="major", axis="both", alpha=0.3, zorder=0)
        ax.set_axisbelow(True)

        ax.set_xlabel("Normalised radius from zoom centre")
        ax.set_ylabel("Background DM fraction")
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=0)
        ax.set_title("Radial DM Contamination (zoom cells)")
        ax.legend(fontsize=8, loc="upper left", framealpha=0.9)
    else:
        ax.text(
            0.5,
            0.5,
            "No DM particles",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
        )
        ax.set_title("Radial DM Contamination (zoom cells)")

    fig.tight_layout()
    path = os.path.join(outdir, "zoom_resolution_quality.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    return fig


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Analyse zoom geometry diagnostic files produced by "
        "SWIFT --zoom-geometry.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Output files:\n"
            "  zoom_cell_types.png            Cell grid by type/subtype\n"
            "  zoom_occupancy_histograms.png  Per-type occupancy histograms\n"
            "  zoom_dm_split.png              DM vs DM_bkg contamination scatter\n"
            "  zoom_padding_analysis.png      Padding extent breakdown\n"
            "  zoom_resolution_quality.png    Occupancy uniformity & radial profiles\n"
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
    figs.append(plot_occupancy_histograms(metadata, data, args.outdir))
    figs.append(plot_dm_split(metadata, data, args.outdir))
    figs.append(plot_padding_analysis(metadata, data, args.outdir))
    figs.append(plot_resolution_quality(metadata, data, args.outdir))

    if args.show:
        plt.show()

    # Close all figures.
    for f in figs:
        if f is not None:
            plt.close(f)


if __name__ == "__main__":
    main()
