"""Plot the cell geometry given a set of cell properties.

This script generates a 2D grid diagram of cells with a nested region.
The script takes parameters from a YAML file to configure the grid structure.

Example:
    python plot_zoom_geometry.py --params params.yaml
"""

import argparse
from dataclasses import dataclass

import h5py
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.lines import Line2D

# Declare the task diff grav constant
zoom_bkg_subdepth_diff_grav_default = 0
zoom_bkg_subdepth_diff_grav = zoom_bkg_subdepth_diff_grav_default


class StoreDictKeyPair(argparse.Action):
    """A class for storing key-value pairs from the command line."""

    def __call__(self, parser, namespace, values, option_string=None):
        """
        Store the key-value pair in the namespace.

        Args:
            parser (argparse.ArgumentParser):
                The parser instance.
            namespace (argparse.Namespace):
                The namespace to store the key-value pair in.
            values (str):
                The key-value pair to store.
            option_string (str):
                The option string.
        """
        key = option_string.lstrip("--")
        if not hasattr(namespace, "swaps"):
            setattr(namespace, "swaps", {})
        getattr(namespace, "swaps")[key] = values


@dataclass
class Space:
    """A dataclass to represent the space properties."""

    box_size: float
    bkg_top_level_cells: int
    zoom_cell_depth: int = 2
    bkg_cdim: np.ndarray = None
    buffer_cell_depth: int = 0
    region_pad_factor: float = 1.1
    zoom_shift: list = None
    com: list = None
    void_lower_bounds: list = None
    void_upper_bounds: list = None
    void_dim: list = None
    zoom_props: dict = None
    bkg_cells: list = None
    zoom_cells: list = None
    buffer_cells: list = None
    cdim: np.ndarray = np.zeros(3, dtype=int)
    width: np.ndarray = np.zeros(3)
    dim: np.ndarray = np.zeros(3)
    iwidth: np.ndarray = np.zeros(3)

    def __post_init__(self):
        """Post-initialization method to set default values."""
        if self.bkg_cdim is None:
            self.bkg_cdim = np.full(3, self.bkg_top_level_cells)
        if self.zoom_shift is None:
            self.zoom_shift = [0.0, 0.0, 0.0]
        if self.com is None:
            self.com = [0.0, 0.0, 0.0]
        if self.void_lower_bounds is None:
            self.void_lower_bounds = [0.0, 0.0, 0.0]
        if self.void_upper_bounds is None:
            self.void_upper_bounds = [0.0, 0.0, 0.0]
        if self.void_dim is None:
            self.void_dim = [0.0, 0.0, 0.0]
        if self.zoom_props is None:
            self.zoom_props = {
                "region_pad_factor": self.region_pad_factor,
                "neighbour_max_tree_depth": -1,
                "bkg_subdepth_diff_grav": zoom_bkg_subdepth_diff_grav_default,
                "bkg_cdim": self.bkg_cdim,
                "buffer_cell_depth": self.buffer_cell_depth,
                "zoom_cell_depth": self.zoom_cell_depth,
                "zoom_shift": self.zoom_shift,
                "com": self.com,
                "void_lower_bounds": self.void_lower_bounds,
                "void_upper_bounds": self.void_upper_bounds,
                "void_dim": self.void_dim,
                "buffer_lower_bounds": [0.0, 0.0, 0.0],
                "buffer_upper_bounds": [0.0, 0.0, 0.0],
                "buffer_dim": [0.0, 0.0, 0.0],
                "buffer_cdim": [0, 0, 0],
                "buffer_width": [0.0, 0.0, 0.0],
                "buffer_iwidth": [0.0, 0.0, 0.0],
                "region_lower_bounds": [0.0, 0.0, 0.0],
                "region_upper_bounds": [0.0, 0.0, 0.0],
                "dim": [0.0, 0.0, 0.0],
                "cdim": [0, 0, 0],
                "width": [0.0, 0.0, 0.0],
                "iwidth": [0.0, 0.0, 0.0],
                "cell_min": 0.0,
                "nr_bkg_cells": 0,
                "nr_buffer_cells": 0,
                "nr_zoom_cells": 0,
                "bkg_cell_offset": 0,
                "buffer_cell_offset": 0,
            }
        self.dim = np.array([self.box_size] * 3)
        if self.bkg_cells is None:
            self.bkg_cells = []
        if self.zoom_cells is None:
            self.zoom_cells = []
        if self.buffer_cells is None:
            self.buffer_cells = []


class Cell:
    """A class to represent a cell in the grid."""

    def __init__(self, x, y, width, ctype="background"):
        """
        Initialize the cell with the given x, y coordinates and width.

        Args:
            x (int): The x coordinate of the lower left corner of the cell.
            y (int): The y coordinate of the lower left corner of the cell.
            width (int): The width of the cell.
        """
        self.x = x
        self.y = y
        self.width = width
        self.type = ctype

    def draw_cell(self, ax):
        """
        Draw the cell on the given axis.

        Args:
            ax (matplotlib.axes.Axes): The axis to draw the cell on.
        """
        # Define the color based on the cell type
        if self.type == "background":
            color = "#bbe6e4"
            order = 0
        elif self.type == "buffer":
            color = "#42bfdd"
            order = 1
        elif self.type == "zoom":
            color = "#084b83"
            order = 2
        else:
            raise ValueError(f"Unknown cell type: {self.type}")

        ax.add_patch(
            plt.Rectangle(
                (self.x, self.y),
                self.width[0],
                self.width[1],
                fill=True,
                facecolor=color,
                edgecolor="black",
                zorder=order,
                alpha=0.9,
            )
        )


def read_params(params_file):
    """
    Read the parameters from the given YAML file.

    Args:
        params_file (str): The path to the YAML file containing the parameters.

    Returns:
        dict: A dictionary containing the parameters read from the file.
    """
    # Read the parameters from the YAML file
    with open(params_file, "r") as file:
        params = yaml.safe_load(file)

    return params


def zoom_parse_params(params):
    """
    Parse parameters for ZoomRegion properties and update Space instance.

    Args:
        params (dict): Dictionary containing the parameters.

    Returns:
        Space: An instance of Space populated with parameters.
    """
    global zoom_bkg_subdepth_diff_grav

    # Extract the initial conditions filename
    ic_filename = params["InitialConditions"]["file_name"]

    # Set the zoom cell depth
    zoom_cell_depth = params["ZoomRegion"].get("zoom_top_level_depth", 2)

    # Set the target background cdim
    bkg_cdim = params["ZoomRegion"].get("bkg_top_level_cells", 16)

    # Get the buffer cell depth
    buffer_cell_depth = params["ZoomRegion"].get("buffer_top_level_depth", 0)

    # Ensure the buffer cell depth is less than the zoom cell depth
    if buffer_cell_depth > zoom_cell_depth:
        raise ValueError("Buffer cell depth must be less than the zoom cell depth.")

    # Extract the zoom width boost factor
    region_pad_factor = params["ZoomRegion"].get("region_pad_factor", 1.1)

    # Extract the minimum difference for background cells
    zoom_bkg_subdepth_diff_grav = params["ZoomRegion"].get(
        "bkg_subdepth_diff_grav", zoom_bkg_subdepth_diff_grav_default
    )

    # Open the HDF5 file and read the box size
    with h5py.File(ic_filename, "r") as f:
        if params["InitialConditions"].get("cleanup_h_factors", False):
            h = params["Cosmology"].get("h", 0.6777)
        else:
            h = 1.0
        box_size = f["Header"].attrs["BoxSize"] / h

    # Create and return the Space instance
    return Space(
        box_size=box_size,
        bkg_top_level_cells=bkg_cdim,
        zoom_cell_depth=zoom_cell_depth,
        buffer_cell_depth=buffer_cell_depth,
        region_pad_factor=region_pad_factor,
    )


def zoom_get_cdim_at_depth(region_dim, parent_width, child_depth):
    """
    Compute the number of child cells in a single parent cell at a given depth.

    Args:
        region_dim (float): The dimension of the region.
        parent_width (float): The width of the parent cell.
        child_depth (int): The depth of the child cell within the parent.

    Returns:
        int: The number of child cells in a single parent cell.
    """
    # How many parent_widths are in the region? (ensure correct rounding)
    region_parent_cdim = np.floor((region_dim + (0.1 * parent_width)) / parent_width)

    # Calculate cdim as the number of parents times the number of children per parent
    return int(region_parent_cdim * (2 ** child_depth))


def zoom_get_region_dim_and_shift(space, params):
    """
    Compute the zoom region centre and boundaries.

    Args:
        space (Space): The space object containing the properties.
        ic_filename (str): The path to the initial conditions HDF5 file.

    Returns:
        float: The maximum side length of the zoom region.
    """
    # Are we cleaning up h factors?
    if params["InitialConditions"].get("cleanup_h_factors", False):
        h = params["Cosmology"].get("h", 0.6777)
    else:
        h = 1.0

    # Open the HDF5 file and read particle data
    with h5py.File(params["InitialConditions"]["file_name"], "r") as f:
        particles = f["PartType1/Coordinates"][:] / h
        masses = f["PartType1/Masses"][:] / h

        # Initialize boundaries and center of mass (COM) calculations
        min_bounds = np.min(particles, axis=0)
        max_bounds = np.max(particles, axis=0)
        com = np.average(particles, axis=0, weights=masses)

        # Calculate initial dimensions and midpoint
        ini_dims = max_bounds - min_bounds
        midpoint = min_bounds + (ini_dims / 2.0)
        box_mid = np.array([space.box_size / 2.0] * 3)

        # Calculate the shift needed to place the midpoint at the center
        # of the box
        zoom_shift = box_mid - midpoint

        # Adjust shift if incremental is too small
        for i in range(3):
            if abs(zoom_shift[i]) < 0.01 * space.box_size:
                zoom_shift[i] = 0.0

        # If the volume isn't periodic, we cannot shift
        for i in range(3):
            if abs(zoom_shift[i]) > 0.01 * space.box_size:
                raise ValueError(
                    "Cannot shift the zoom region to the centre of the box "
                    "when the box is not periodic."
                )

        # Update space with computed shift and COM
        space.zoom_shift = zoom_shift.tolist()
        space.zoom_props["zoom_shift"] = zoom_shift.tolist()
        space.com = (com + zoom_shift).tolist()
        space.zoom_props["com"] = (com + zoom_shift).tolist()

        # Compute the maximum side length of the zoom region
        ini_dim = max(ini_dims)

    return ini_dim


def zoom_get_void_geometry(space, region_dim):
    """
    Compute the void region geometry.

    The void region is the region covered by background cells above the zoom
    region. If the void region is sufficiently close to the zoom region size,
    then the two will be made equivalent later on. Otherwise, the void region is
    equivalent to the buffer region.

    Args:
        space (Space): The space object containing the properties.
        region_dim (float): The dimension of the zoom region before tessellating the volume.

    Returns:
        int: The number of zoom regions that tessellate the void region.
    """
    # Get the lower and upper bounds of the zoom region based on initial dimensions and padding factor
    lower_bounds = [(space.box_size / 2) - (region_dim / 2.0) for _ in range(3)]
    upper_bounds = [(space.box_size / 2) + (region_dim / 2.0) for _ in range(3)]

    # Assign these bounds to the zoom region properties temporarily
    space.zoom_props["region_lower_bounds"] = lower_bounds
    space.zoom_props["region_upper_bounds"] = upper_bounds

    # Find the background cell edges that contain these bounds
    print(lower_bounds, space.box_size, space.bkg_cdim)
    void_lower_bounds = [
        np.floor(lower_bounds[i] * space.bkg_cdim[i] / space.box_size)
        * (space.box_size / space.bkg_cdim[i])
        for i in range(3)
    ]
    void_upper_bounds = [
        (np.floor(upper_bounds[i] * space.bkg_cdim[i] / space.box_size) + 1)
        * (space.box_size / space.bkg_cdim[i])
        for i in range(3)
    ]

    # Assign the void bounds
    space.void_lower_bounds = void_lower_bounds
    space.void_upper_bounds = void_upper_bounds

    # Compute the void region dimensions
    space.void_dim = [void_upper_bounds[i] - void_lower_bounds[i] for i in range(3)]

    # Compute the number of zoom regions that tessellate the void region
    nr_zoom_regions = int(np.ceil(space.void_dim[0] / region_dim))

    return nr_zoom_regions


def zoom_get_geometry_no_buffer_cells(space):
    """
    Compute the geometry when no buffer cells are needed.

    Args:
        space (Space): The space object containing the properties.
    """
    # If we have a buffer cell depth, warn that we will ignore it
    if space.buffer_cell_depth > 0:
        print("No buffer cells are needed, ignoring buffer cell depth.")
        space.buffer_cell_depth = 0

    # Zero the buffer region properties explicitly
    space.zoom_props["buffer_lower_bounds"] = [0.0, 0.0, 0.0]
    space.zoom_props["buffer_upper_bounds"] = [0.0, 0.0, 0.0]
    space.zoom_props["buffer_dim"] = [0.0, 0.0, 0.0]
    space.zoom_props["buffer_cdim"] = [0, 0, 0]
    space.zoom_props["buffer_width"] = [0.0, 0.0, 0.0]

    # Match the zoom region bounds to the void region bounds
    space.zoom_props["region_lower_bounds"] = space.void_lower_bounds
    space.zoom_props["region_upper_bounds"] = space.void_upper_bounds

    # Compute the zoom region dimensions
    space.zoom_props["dim"] = [
        space.zoom_props["region_upper_bounds"][i]
        - space.zoom_props["region_lower_bounds"][i]
        for i in range(3)
    ]

    # Compute the number of zoom cells in the void region
    cdim = zoom_get_cdim_at_depth(
        space.zoom_props["dim"][0],
        space.box_size / space.bkg_top_level_cells,
        space.zoom_cell_depth,
    )

    # Compute the zoom cdim and cell width
    space.zoom_props["cdim"] = [cdim] * 3
    space.zoom_props["width"] = [space.zoom_props["dim"][i] / cdim for i in range(3)]
    space.zoom_props["iwidth"] = [1.0 / space.zoom_props["width"][i] for i in range(3)]


def zoom_get_geometry_with_buffer_cells(space):
    """
    Compute the geometry when buffer cells are needed.

    Args:
        space (Space): The space object containing the properties.
    """
    # Ensure we have a buffer cell depth
    if space.buffer_cell_depth == 0:
        zoom_report_cell_properties(space)
        raise ValueError(
            "Current cell structure requires buffer cells but no buffer cell depth has "
            "been given. ZoomRegion:buffer_top_level_depth must be greater than 0."
        )

    # Match the buffer region bounds to the void region bounds
    space.zoom_props["buffer_lower_bounds"] = space.void_lower_bounds[:]
    space.zoom_props["buffer_upper_bounds"] = space.void_upper_bounds[:]

    # Compute the buffer region dimensions
    space.zoom_props["buffer_dim"] = [
        space.zoom_props["buffer_upper_bounds"][i]
        - space.zoom_props["buffer_lower_bounds"][i]
        for i in range(3)
    ]

    # Compute the number of buffer cells in the void region
    buffer_cdim = zoom_get_cdim_at_depth(
        space.zoom_props["buffer_dim"][0],
        space.box_size / space.bkg_top_level_cells,
        space.buffer_cell_depth,
    )

    # Compute the buffer cdim and cell width
    space.zoom_props["buffer_cdim"] = [buffer_cdim] * 3
    space.zoom_props["buffer_width"] = [
        space.zoom_props["buffer_dim"][i] / buffer_cdim for i in range(3)
    ]
    space.zoom_props["buffer_iwidth"] = [
        1.0 / space.zoom_props["buffer_width"][i] for i in range(3)
    ]

    # Find the buffer cell edges that contain the zoom region bounds
    region_lower_bounds = []
    region_upper_bounds = []
    for i in range(3):
        lower = int(
            np.floor(
                (
                    space.zoom_props["region_lower_bounds"][i]
                    - space.zoom_props["buffer_lower_bounds"][i]
                )
                * space.zoom_props["buffer_iwidth"][i]
            )
        )
        upper = int(
            np.floor(
                (
                    space.zoom_props["region_upper_bounds"][i]
                    - space.zoom_props["buffer_lower_bounds"][i]
                )
                * space.zoom_props["buffer_iwidth"][i]
            )
        )
        region_lower_bounds.append(
            lower * space.zoom_props["buffer_width"][i]
            + space.zoom_props["buffer_lower_bounds"][i]
        )
        region_upper_bounds.append(
            (upper + 1) * space.zoom_props["buffer_width"][i]
            + space.zoom_props["buffer_lower_bounds"][i]
        )

    # Assign the new aligned zoom bounds
    space.zoom_props["region_lower_bounds"] = region_lower_bounds
    space.zoom_props["region_upper_bounds"] = region_upper_bounds

    # Compute the zoom region dimensions
    space.zoom_props["dim"] = [
        region_upper_bounds[i] - region_lower_bounds[i] for i in range(3)
    ]

    # Compute the number of zoom cells in the zoom region
    cdim = zoom_get_cdim_at_depth(
        space.zoom_props["dim"][0],
        space.zoom_props["buffer_width"][0],
        space.zoom_cell_depth - space.buffer_cell_depth,
    )

    # Compute the zoom cdim and cell width
    space.zoom_props["cdim"] = [cdim] * 3
    space.zoom_props["width"] = [space.zoom_props["dim"][i] / cdim for i in range(3)]
    space.zoom_props["iwidth"] = [1.0 / space.zoom_props["width"][i] for i in range(3)]


def zoom_report_cell_properties(space):
    """
    Report Zoom Region Properties.

    This function prints out a table containing the properties of the
    zoom region, if it is enabled. The table includes information such as
    dimensions, center, CDIM, background CDIM, buffer CDIM, region buffer
    ratio, zoom boost factor, minimum zoom cell width, background cell width,
    buffer width, and the number of wanderers.

    Args:
        space (Space): The space object containing the properties.
    """
    zoom_props = space.zoom_props

    # Cdims
    print(
        f"{'Background cdim':>28} = [{space.bkg_cdim[0]}, {space.bkg_cdim[1]}, {space.bkg_cdim[2]}]"
    )
    if zoom_props.get("with_buffer_cells"):
        print(
            f"{'Buffer cdim':>28} = [{zoom_props['buffer_cdim'][0]}, {zoom_props['buffer_cdim'][1]}, {zoom_props['buffer_cdim'][2]}]"
        )
    print(
        f"{'Zoom cdim':>28} = [{zoom_props['cdim'][0]}, {zoom_props['cdim'][1]}, {zoom_props['cdim'][2]}]"
    )

    # Dimensions
    print(
        f"{'Background Dimensions':>28} = [{space.box_size}, {space.box_size}, {space.box_size}]"
    )
    if zoom_props.get("with_buffer_cells"):
        print(
            f"{'Buffer Region Dimensions':>28} = [{zoom_props['buffer_dim'][0]}, {zoom_props['buffer_dim'][1]}, {zoom_props['buffer_dim'][2]}]"
        )
    print(
        f"{'Zoom Region Dimensions':>28} = [{zoom_props['dim'][0]}, {zoom_props['dim'][1]}, {zoom_props['dim'][2]}]"
    )

    # Cell Widths
    print(
        f"{'Background Cell Width':>28} = [{space.box_size / space.bkg_top_level_cells}, {space.box_size / space.bkg_top_level_cells}, {space.box_size / space.bkg_top_level_cells}]"
    )
    if zoom_props.get("with_buffer_cells"):
        print(
            f"{'Buffer Cell Width':>28} = [{zoom_props['buffer_width'][0]}, {zoom_props['buffer_width'][1]}, {zoom_props['buffer_width'][2]}]"
        )
    print(
        f"{'Zoom Cell Width':>28} = [{zoom_props['width'][0]}, {zoom_props['width'][1]}, {zoom_props['width'][2]}]"
    )

    # Number of Cells
    print(f"{'Number of Background Cells':>28} = {zoom_props.get('nr_bkg_cells', 0)}")
    if zoom_props.get("with_buffer_cells"):
        print(
            f"{'Number of Buffer Cells':>28} = {zoom_props.get('nr_buffer_cells', 0)}"
        )
    print(f"{'Number of Zoom Cells':>28} = {zoom_props.get('nr_zoom_cells', 0)}")

    # Bounds
    if zoom_props.get("with_buffer_cells"):
        print(
            f"{'Buffer Bounds':>28} = [{zoom_props['buffer_lower_bounds'][0]}-{zoom_props['buffer_upper_bounds'][0]}, {zoom_props['buffer_lower_bounds'][1]}-{zoom_props['buffer_upper_bounds'][1]}, {zoom_props['buffer_lower_bounds'][2]}-{zoom_props['buffer_upper_bounds'][2]}]"
        )
    print(
        f"{'Zoom Region Bounds':>28} = [{zoom_props['region_lower_bounds'][0]}-{zoom_props['region_upper_bounds'][0]}, {zoom_props['region_lower_bounds'][1]}-{zoom_props['region_upper_bounds'][1]}, {zoom_props['region_lower_bounds'][2]}-{zoom_props['region_upper_bounds'][2]}]"
    )

    # Depths
    if zoom_props.get("with_buffer_cells"):
        print(f"{'Buffer Top Level Depth':>28} = {zoom_props['buffer_cell_depth']}")
    print(f"{'Zoom Top Level Depth':>28} = {zoom_props['zoom_cell_depth']}")

    # Assorted extra zoom properties
    print(f"{'Zoom Region Pad Factor':>28} = {zoom_props['region_pad_factor']}")
    print(
        f"{'Zoom Region Shift':>28} = [{zoom_props['zoom_shift'][0]}, {zoom_props['zoom_shift'][1]}, {zoom_props['zoom_shift'][2]}]"
    )
    print(
        f"{'Zoom Region Center':>28} = [{zoom_props['region_lower_bounds'][0] + (zoom_props['dim'][0] / 2.0)}, {zoom_props['region_lower_bounds'][1] + (zoom_props['dim'][1] / 2.0)}, {zoom_props['region_lower_bounds'][2] + (zoom_props['dim'][2] / 2.0)}]"
    )


def zoom_report_cell_properties_text(space):
    """
    Report Zoom Region Properties.

    Args:
        space (Space): The space object containing the properties.

    Returns:
        str: The formatted string containing the zoom region properties.
    """
    zoom_props = space.zoom_props
    report_lines = []

    # Cdims
    report_lines.append(
        f"{'Background cdim':>28} = [{space.bkg_cdim[0]}, {space.bkg_cdim[1]}, {space.bkg_cdim[2]}]"
    )
    if zoom_props.get("with_buffer_cells"):
        report_lines.append(
            f"{'Buffer cdim':>28} = [{zoom_props['buffer_cdim'][0]}, {zoom_props['buffer_cdim'][1]}, {zoom_props['buffer_cdim'][2]}]"
        )
    report_lines.append(
        f"{'Zoom cdim':>28} = [{zoom_props['cdim'][0]}, {zoom_props['cdim'][1]}, {zoom_props['cdim'][2]}]"
    )

    # Dimensions
    report_lines.append(
        f"{'Background Dimensions':>28} = [{space.box_size}, {space.box_size}, {space.box_size}]"
    )
    if zoom_props.get("with_buffer_cells"):
        report_lines.append(
            f"{'Buffer Region Dimensions':>28} = [{zoom_props['buffer_dim'][0]}, {zoom_props['buffer_dim'][1]}, {zoom_props['buffer_dim'][2]}]"
        )
    report_lines.append(
        f"{'Zoom Region Dimensions':>28} = [{zoom_props['dim'][0]}, {zoom_props['dim'][1]}, {zoom_props['dim'][2]}]"
    )

    # Cell Widths
    report_lines.append(
        f"{'Background Cell Width':>28} = [{space.box_size / space.bkg_top_level_cells}, {space.box_size / space.bkg_top_level_cells}, {space.box_size / space.bkg_top_level_cells}]"
    )
    if zoom_props.get("with_buffer_cells"):
        report_lines.append(
            f"{'Buffer Cell Width':>28} = [{zoom_props['buffer_width'][0]}, {zoom_props['buffer_width'][1]}, {zoom_props['buffer_width'][2]}]"
        )
    report_lines.append(
        f"{'Zoom Cell Width':>28} = [{zoom_props['width'][0]}, {zoom_props['width'][1]}, {zoom_props['width'][2]}]"
    )

    # Number of Cells
    report_lines.append(
        f"{'Number of Background Cells':>28} = {zoom_props.get('nr_bkg_cells', 0)}"
    )
    if zoom_props.get("with_buffer_cells"):
        report_lines.append(
            f"{'Number of Buffer Cells':>28} = {zoom_props.get('nr_buffer_cells', 0)}"
        )
    report_lines.append(
        f"{'Number of Zoom Cells':>28} = {zoom_props.get('nr_zoom_cells', 0)}"
    )

    # Bounds
    if zoom_props.get("with_buffer_cells"):
        report_lines.append(
            f"{'Buffer Lower Bounds':>28} = [{zoom_props['buffer_lower_bounds'][0]}, {zoom_props['buffer_lower_bounds'][1]}, {zoom_props['buffer_lower_bounds'][2]}]"
        )
        report_lines.append(
            f"{'Buffer Upper Bounds':>28} = [{zoom_props['buffer_upper_bounds'][0]}, {zoom_props['buffer_upper_bounds'][1]}, {zoom_props['buffer_upper_bounds'][2]}]"
        )
    report_lines.append(
        f"{'Zoom Region Lower Bounds':>28} = [{zoom_props['region_lower_bounds'][0]}, {zoom_props['region_lower_bounds'][1]}, {zoom_props['region_lower_bounds'][2]}]"
    )
    report_lines.append(
        f"{'Zoom Region Upper Bounds':>28} = [{zoom_props['region_upper_bounds'][0]}, {zoom_props['region_upper_bounds'][1]}, {zoom_props['region_upper_bounds'][2]}]"
    )

    # Depths
    if zoom_props.get("with_buffer_cells"):
        report_lines.append(
            f"{'Buffer Top Level Depth':>28} = {zoom_props['buffer_cell_depth']}"
        )
    report_lines.append(
        f"{'Zoom Top Level Depth':>28} = {zoom_props['zoom_cell_depth']}"
    )

    # Assorted extra zoom properties
    report_lines.append(
        f"{'Zoom Region Pad Factor':>28} = {zoom_props['region_pad_factor']}"
    )
    report_lines.append(
        f"{'Zoom Region Shift':>28} = [{zoom_props['zoom_shift'][0]}, {zoom_props['zoom_shift'][1]}, {zoom_props['zoom_shift'][2]}]"
    )
    report_lines.append(
        f"{'Zoom Region Center':>28} = [{zoom_props['region_lower_bounds'][0] + (zoom_props['dim'][0] / 2.0)}, {zoom_props['region_lower_bounds'][1] + (zoom_props['dim'][1] / 2.0)}, {zoom_props['region_lower_bounds'][2] + (zoom_props['dim'][2] / 2.0)}]"
    )

    return "\n".join(report_lines)


def zoom_region_init(space, params, verbose=False):
    """
    Initialise the zoom region geometry.

    This will compute the cell grid properties ready for cell construction
    when zoom_construct_tl_cells.

    Args:
        space (Space): The space object to update.
        params (dict): The parameters to use for the zoom region.
        verbose (bool): Whether to print verbose output.
    """
    # Compute the extent of the zoom region and calculate the shift to move
    # the zoom region to the centre of the box
    ini_dim = zoom_get_region_dim_and_shift(space, params)

    # Apply the shift to the particles
    for particle_group in ["parts", "gparts", "sparts", "bparts", "sinks"]:
        for particle in getattr(space, particle_group, []):
            particle["x"][0] += space.zoom_props["zoom_shift"][0]
            particle["x"][1] += space.zoom_props["zoom_shift"][1]
            particle["x"][2] += space.zoom_props["zoom_shift"][2]

    # Include the requested padding around the high-resolution particles
    max_dim = ini_dim * space.zoom_props["region_pad_factor"]

    # Define the background grid (we'll treat this as gospel)
    for i in range(3):
        space.cdim[i] = space.zoom_props["bkg_cdim"][i]
        space.width[i] = space.dim[i] / space.cdim[i]
        space.iwidth[i] = 1.0 / space.width[i]

    # Compute the void region bounds and number of zoom regions that tessellate it
    nr_zoom_regions = zoom_get_void_geometry(space, max_dim)

    # Check the user gave a sensible background cdim
    if nr_zoom_regions >= 64:
        zoom_report_cell_properties(space)
        raise ValueError(
            "Background cell size is too large relative to the zoom region! "
            "Increase ZoomRegion:bkg_top_level_cells (would have needed {} zoom "
            "cells in the void region).".format(nr_zoom_regions)
        )

    # Warn the user if the number of zoom regions is large but not unreasonable
    if nr_zoom_regions >= 16:
        if verbose:
            print(
                "Warning: Background cell size is large relative to the zoom region! "
                "(we'll need at least {} buffer cells which may be slow).".format(
                    nr_zoom_regions
                )
            )

    # Decide whether to use buffer cells based on the padding due to background cells
    if nr_zoom_regions <= 2:
        space.zoom_props["with_buffer_cells"] = 0
        zoom_get_geometry_no_buffer_cells(space)
    else:
        space.zoom_props["with_buffer_cells"] = 1
        zoom_get_geometry_with_buffer_cells(space)

    # Store what the true boost factor ended up being
    input_pad_factor = space.zoom_props["region_pad_factor"]
    space.zoom_props["region_pad_factor"] = space.zoom_props["dim"][0] / ini_dim

    # Ensure the zoom region is not smaller than the high-resolution particle distribution
    if space.zoom_props["dim"][0] < ini_dim:
        zoom_report_cell_properties(space)
        raise ValueError(
            "Found a zoom region smaller than the high-resolution particle "
            "distribution! Adjust the cell structure "
            "(ZoomRegion:bkg_top_level_cells, ZoomRegion:zoom_top_level_cells)"
        )

    # Warn if the size of the requested padding region has drastically changed
    if (space.zoom_props["region_pad_factor"] / input_pad_factor) >= 2:
        if verbose:
            print(
                "Warning: The pad region has to be {} times larger than requested. "
                "Either increase ZoomRegion:region_pad_factor, increase the number of "
                "background cells, or increase the depths of the zoom cells (and buffer cells if using).".format(
                    int(space.zoom_props["region_pad_factor"] / input_pad_factor)
                )
            )

    # Set the neighbour max tree depth if not explicitly provided
    if space.zoom_props["neighbour_max_tree_depth"] < 0:
        space.zoom_props["neighbour_max_tree_depth"] = space.zoom_props[
            "zoom_cell_depth"
        ]

    # Set the minimum allowed zoom cell width
    zoom_dmax = max(space.zoom_props["dim"])
    space.zoom_props["cell_min"] = 0.99 * zoom_dmax / space.zoom_props["cdim"][0]

    # Set the minimum background cell size
    dmax = max(space.dim)
    space.cell_min = 0.99 * dmax / space.cdim[0]

    # Store cell numbers and offsets
    space.zoom_props["bkg_cell_offset"] = (
        space.zoom_props["cdim"][0]
        * space.zoom_props["cdim"][1]
        * space.zoom_props["cdim"][2]
    )
    space.zoom_props["nr_zoom_cells"] = space.zoom_props["bkg_cell_offset"]
    space.zoom_props["nr_bkg_cells"] = space.cdim[0] * space.cdim[1] * space.cdim[2]
    space.zoom_props["buffer_cell_offset"] = (
        space.zoom_props["bkg_cell_offset"] + space.zoom_props["nr_bkg_cells"]
    )
    space.zoom_props["nr_buffer_cells"] = (
        space.zoom_props["buffer_cdim"][0]
        * space.zoom_props["buffer_cdim"][1]
        * space.zoom_props["buffer_cdim"][2]
    )


def construct_tl_cells(space):
    """
    Construct the top level cells.

    Args:
        space (Space): The space object containing the properties.
    """
    # Create the background cells
    space.bkg_cells = []
    for i in range(space.cdim[0]):
        for j in range(space.cdim[1]):
            cell = Cell(
                i * space.width[0], j * space.width[1], space.width, ctype="background"
            )
            space.bkg_cells.append(cell)

    # Create the zoom cells
    space.zoom_cells = []
    for i in range(space.zoom_props["cdim"][0]):
        for j in range(space.zoom_props["cdim"][1]):
            cell = Cell(
                i * space.zoom_props["width"][0]
                + space.zoom_props["region_lower_bounds"][0],
                j * space.zoom_props["width"][1]
                + space.zoom_props["region_lower_bounds"][1],
                space.zoom_props["width"],
                ctype="zoom",
            )
            space.zoom_cells.append(cell)

    # Create the buffer cells
    if space.zoom_props.get("with_buffer_cells"):
        space.buffer_cells = []
        for i in range(space.zoom_props["buffer_cdim"][0]):
            for j in range(space.zoom_props["buffer_cdim"][1]):
                cell = Cell(
                    i * space.zoom_props["buffer_width"][0]
                    + space.zoom_props["buffer_lower_bounds"][0],
                    j * space.zoom_props["buffer_width"][1]
                    + space.zoom_props["buffer_lower_bounds"][1],
                    space.zoom_props["buffer_width"],
                    ctype="buffer",
                )
                space.buffer_cells.append(cell)


def draw_grid(space, params):
    """
    Draw the grid diagram of the cells.

    Args:
        space (Space): The space object containing the properties.
    """
    # Setup the figure
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")
    ax.set_xlim(0, space.box_size)
    ax.set_ylim(0, space.box_size)

    # Draw each cell
    for cell in space.bkg_cells:
        cell.draw_cell(ax)
    for cell in space.buffer_cells:
        cell.draw_cell(ax)
    for cell in space.zoom_cells:
        cell.draw_cell(ax)

    # Define handles for the legend as square markers
    handles = [
        Line2D(
            [0],
            [0],
            marker="s",  # Use a square marker
            color="w",  # Set the color of the line to be invisible
            markerfacecolor="#bbe6e4",
            markeredgecolor="black",
            markersize=15,
            label="Background",
            alpha=0.9,
        )
    ]

    if len(space.buffer_cells) > 0:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="s",  # Use a square marker
                color="w",  # Set the color of the line to be invisible
                markerfacecolor="#42bfdd",
                markeredgecolor="black",
                markersize=15,
                label="Buffer",
                alpha=0.9,
            )
        )

    handles.append(
        Line2D(
            [0],
            [0],
            marker="s",  # Use a square marker
            color="w",  # Set the color of the line to be invisible
            markerfacecolor="#084b83",
            markeredgecolor="black",
            markersize=15,
            label="Zoom",
            alpha=0.9,
        )
    )
    # Place a legend below the x axis with the cell types
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=len(handles),
        frameon=False,
    )

    # Save the figure including relevant information in the filename
    ic_name = params["InitialConditions"]["file_name"].split("/")[-1].split(".")[0]
    fig.savefig(f"zoom_geometry_{ic_name}.png", dpi=300)
    plt.show()
    plt.close(fig)


def main():
    """
    Draw a 2D slice of the cell grid geometry.

    Example:
        python plot_zoom_geometry.py --params params.yaml
    """
    parser = argparse.ArgumentParser(
        description="Draw a 2D grid diagram of cells with a nested region."
    )
    parser.add_argument(
        dest="params", type=str, help="Path to the YAML file containing the parameters."
    )
    # Add arbitrary arguments
    parser.add_argument(
        "--",
        dest="overrides",
        action=StoreDictKeyPair,
        nargs=1,
        metavar="VALUE",
        help="Override the parameter file with the given key-value pair.",
    )

    # Unpack the arguments
    args, unknown_args = parser.parse_known_args()
    args.overrides = {}

    # Process unknown_args manually
    while unknown_args:
        if unknown_args[0].startswith("--"):
            key = unknown_args.pop(0).lstrip("--")
            if (
                unknown_args
                and not unknown_args[0].startswith("-")
                and not unknown_args[0].startswith("--")
            ):
                value = unknown_args.pop(0)
                args.overrides[key] = value
            else:
                args.overrides[key] = None
        else:
            unknown_args.pop(0)

    # Read the parameters from the YAML file
    params = read_params(args.params)

    # Override the parameters with the provided key-value pairs
    for key, value in args.overrides.items():
        if key not in params["ZoomRegion"]:
            raise ValueError(f"Invalid parameter: {key}")
        else:
            print(f"Overriding parameter: {key} = {value}")

        # Handle whether this is an integer or float
        if "." in value:
            params["ZoomRegion"][key] = float(value)
        else:
            params["ZoomRegion"][key] = int(value)

    # Parse the parameter file and populate the Space instance
    space = zoom_parse_params(params)

    # Initialize the zoom region geometry
    zoom_region_init(space, params, verbose=True)
    zoom_report_cell_properties(space)

    # Create the top level cells
    construct_tl_cells(space)

    # Draw the grid
    draw_grid(space, params)


if __name__ == "__main__":
    main()
