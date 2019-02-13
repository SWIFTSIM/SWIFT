"""
Makes an image of the Santa Barbara cluster.

Requires py-sphviewer.

Invoke as follows:

python3 makeImage.py <name of hdf5 file> \
                     <number of particle type (i.e. 0 or 1)> \
                     <colour map to use (default viridis)> \
                     <text color (default white)> \
                     <image resolution (default 2048x2048)>
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib

from sphviewer.tools import QuickView
from matplotlib.patches import Rectangle

from typing import Tuple
from collections import namedtuple


# Set up our simulation data collection to keep stuff together
SimulationData = namedtuple(
    "SimulationData",
    ["coordinates", "masses", "sph_name", "dark_matter_mass", "swift_name", "boxsize"],
)


def latex_float(f):
    """
    Taken from:
    https://stackoverflow.com/questions/13490292/format-number-using-latex-notation-in-python.
    
    Formats a float to LaTeX style.
    """

    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


def read_data_from_file(filename: str, part_type=0) -> SimulationData:
    """
    Reads the relevant data from the HDF5 file.
    """
    part_type_name = f"PartType{part_type}"

    with h5py.File(filename, "r") as file:
        coordinates, boxsize = move_box(file[f"{part_type_name}/Coordinates"][...].T)
        masses = file[f"{part_type_name}/Masses"][...]

        sph_name = file["HydroScheme"].attrs["Scheme"].decode("utf-8")
        unit_mass = (
            float(file["Units"].attrs["Unit mass in cgs (U_M)"]) / 2e33
        )  # in M_sun

        dark_matter_mass = float(file["PartType1/Masses"][0]) * unit_mass

        code_revision = file["Code"].attrs["Git Revision"].decode("utf-8")
        swift_name = f"SWIFT {code_revision}"

        data = SimulationData(
            coordinates=coordinates,
            masses=masses,
            sph_name=sph_name,
            dark_matter_mass=dark_matter_mass,
            swift_name=swift_name,
            boxsize=boxsize,
        )

    return data


def move_box(coordinates: np.ndarray) -> np.ndarray:
    """
    Takes the coordinates and moves them in the x-y plane. This moves them 20
    code units to the left/right to ensure that the zoomed-out version of the
    cluster image is nicely shown
    """

    boxsize = np.max(coordinates[0])
    coordinates[0] -= 20
    coordinates[1] -= 20
    coordinates[0] %= boxsize
    coordinates[1] %= boxsize

    return coordinates, boxsize


def generate_views(data: SimulationData, res=2048) -> Tuple[np.ndarray]:
    """
    Generates the views on the data from py-sphviewer.

    Returns the overall image for the whole box and then a zoomed region.
    """

    qv_all = QuickView(
        data.coordinates,
        data.masses,
        r="infinity",
        plot=False,
        xsize=res,
        ysize=res,
        logscale=False,
        p=0,
        np=48,
    )
    zoomed_res = (res * 6) // 10
    mask = np.logical_and(
        np.logical_and(
            data.coordinates[0] > (data.boxsize/2-4-20),
            data.coordinates[0] < (data.boxsize/2+6-20)
        ),
        np.logical_and(
            data.coordinates[1] > (data.boxsize/2-3.5-20),
            data.coordinates[1] < (data.boxsize/2+6.5-20)
        )
    )
    qv_zoomed = QuickView(
        data.coordinates.T[mask].T,
        data.masses[mask],
        r="infinity",
        plot=False,
        xsize=zoomed_res,
        ysize=zoomed_res,
        logscale=False,
        np=48,
    )

    return qv_all.get_image(), qv_zoomed.get_image()


def create_plot(data: SimulationData, res=2048, cmap="viridis", text_color="white"):
    """
    Creates a figure and axes object and returns them for you to do with what you wish.
    """

    img_all, img_zoomed = generate_views(data, res)

    fig, ax = plt.subplots(figsize=(8, 8))

    # Set up in "image" mode
    ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1)

    ax.imshow(
        np.log10(img_all + np.min(img_all[img_all != 0])),
        origin="lower",
        extent=[-1, 1, -1, 1],
        cmap=cmap,
    )

    lower_left = [(-24 / (0.5 * data.boxsize)), (-23.5 / (0.5 * data.boxsize))]
    zoom_rect = Rectangle(
        lower_left,
        10 / (0.5 * data.boxsize),
        10 / (0.5 * data.boxsize),
        linewidth=2,
        edgecolor=text_color,
        facecolor="none",
    )
    ax.add_patch(zoom_rect)

    # Remove ticks as we want "image mode"
    ax2 = fig.add_axes([0.35, 0.35, 0.6, 0.6], frame_on=True, xticks=[], yticks=[])

    ax2.imshow(
        np.log10(img_zoomed + np.min(img_zoomed[img_zoomed != 0])),
        origin="lower",
        extent=[-1, 1, -1, 1],
        cmap=cmap,
    )

    # This ugly hack sets the box around the subfigure to be white
    for child in ax2.get_children():
        if isinstance(child, matplotlib.spines.Spine):
            child.set_color(text_color)
            child.set_linewidth(2)

    # Draw lines between boxes

    # Bottom Right
    ax.plot(
        [(-14 / (0.5 * data.boxsize)), 0.9],
        [(-23.5 / (0.5 * data.boxsize)), -0.3],
        lw=2,
        color=text_color,
    )
    # Top Left
    ax.plot(
        [(-24 / (0.5 * data.boxsize)), -0.3],
        [(-13.5 / (0.5 * data.boxsize)), 0.9],
        lw=2,
        color=text_color,
    )

    ax.text(0.95, -0.95, data.swift_name, color=text_color, ha="right")
    formatted_dark_matter_mass = latex_float(data.dark_matter_mass)
    ax.text(
        -0.95,
        0.95,
        rf"M$_{{\rm DM}} = {formatted_dark_matter_mass}$ M$_\odot$",
        color=text_color,
        va="top",
    )
    ax.text(
        -0.95,
        -0.95,
        data.sph_name + "\n" + r"Santa Barbara Cluster (re-ran from Frenk+ 1999)",
        color=text_color,
    )

    return fig, ax


if __name__ == "__main__":
    import sys

    try:
        filename = sys.argv[1]
    except IndexError:
        filename = "santabarbara_0153.hdf5"

    try:
        part_type = int(sys.argv[2])
    except IndexError:
        part_type = 0

    try:
        cmap = sys.argv[3]
    except IndexError:
        cmap = "viridis"

    try:
        text_color = sys.argv[4]
    except IndexError:
        text_color = "white"

    try:
        res = int(sys.argv[5])
    except IndexError:
        res = 2048

    # Read in the data from file

    try:
        data = read_data_from_file(filename, part_type)
    except IndexError:
        # Must be a dark matter only run
        part_type = 1
        data = read_data_from_file(filename, part_type)

    # Make the plot

    fig, ax = create_plot(data, res, cmap, text_color)

    fig.savefig(
        f"SantaBarbara_{data.sph_name[:8]}_{cmap}_PartType{part_type}_res{res}.png",
        dpi=res // 8,
    )
