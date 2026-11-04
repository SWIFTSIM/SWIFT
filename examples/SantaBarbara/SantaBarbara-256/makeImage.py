"""
Makes an image of the Santa Barbara cluster.

Invoke as follows:

python3 makeImage.py <name of hdf5 file> \
                     <number of particle type (i.e. 0 or 1)> \
                     <colour map to use (default viridis)> \
                     <text color (default white)> \
                     <image resolution (default 2048x2048)>
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import swiftsimio as sw

from swiftsimio.visualisation import project_pixel_grid, generate_smoothing_lengths
from matplotlib.patches import Rectangle

from typing import Tuple
from collections import namedtuple

# Set up our simulation data collection to keep stuff together
SimulationData = namedtuple(
    "SimulationData",
    ["particles", "sph_name", "dark_matter_mass", "swift_name", "boxsize"],
)

# The cluster sits at the centre of the box. The full-box image is shifted by
# this many code units (Mpc) in x and y, using the periodicity, so that the
# cluster appears in the lower-left quadrant and leaves room for the zoom.
SHIFT = 20.0


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
    Reads the relevant data from the snapshot using swiftsimio.

    Dark matter particles carry no smoothing length, so one is generated
    from their 32 nearest neighbours.
    """
    data = sw.load(filename)

    if part_type == 0:
        particles = data.gas
    else:
        particles = data.dark_matter
        particles.smoothing_lengths = generate_smoothing_lengths(
            particles.coordinates, data.metadata.boxsize, kernel_gamma=1.8
        )

    sph_name = data.metadata.hydro_scheme["Scheme"].decode("utf-8")
    dark_matter_mass = float(data.dark_matter.masses[0].to_value("Msun"))

    code_revision = data.metadata.code["Git Revision"].decode("utf-8")
    swift_name = f"SWIFT {code_revision}"

    return SimulationData(
        particles=particles,
        sph_name=sph_name,
        dark_matter_mass=dark_matter_mass,
        swift_name=swift_name,
        boxsize=float(data.metadata.boxsize[0].to_value(particles.coordinates.units)),
    )


def make_region(data: SimulationData, x_min, x_max, y_min, y_max) -> sw.cosmo_array:
    """
    Builds a projection region (in code units) compatible with the particles.
    """
    coordinates = data.particles.coordinates
    return sw.cosmo_array(
        [x_min, x_max, y_min, y_max],
        units=coordinates.units,
        comoving=coordinates.comoving,
        scale_factor=coordinates.cosmo_factor.scale_factor,
        scale_exponent=1,
    )


def generate_views(data: SimulationData, res=2048) -> Tuple[np.ndarray]:
    """
    Generates mass projections of the data along the z axis with swiftsimio.

    Returns the overall image for the whole box and then a zoomed region.
    The images are transposed so that they can be passed directly to imshow.
    """

    # The whole box, shifted by SHIFT in x and y (wraps periodically)
    region_all = make_region(
        data, SHIFT, data.boxsize + SHIFT, SHIFT, data.boxsize + SHIFT
    )
    img_all = project_pixel_grid(
        data.particles, resolution=res, region=region_all, parallel=True
    )

    # A 10x10 code-unit region around the cluster
    zoomed_res = (res * 6) // 10
    region_zoomed = make_region(
        data,
        data.boxsize / 2 - 4,
        data.boxsize / 2 + 6,
        data.boxsize / 2 - 3.5,
        data.boxsize / 2 + 6.5,
    )
    img_zoomed = project_pixel_grid(
        data.particles, resolution=zoomed_res, region=region_zoomed, parallel=True
    )

    return (
        img_all.to_comoving_value("Msun/Mpc**2").T,
        img_zoomed.to_comoving_value("Msun/Mpc**2").T,
    )


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

    if part_type == 0 and sw.load(filename).metadata.n_gas == 0:
        # Must be a dark matter only run
        part_type = 1

    data = read_data_from_file(filename, part_type)

    # Make the plot

    fig, ax = create_plot(data, res, cmap, text_color)

    fig.savefig(
        f"SantaBarbara_{data.sph_name[:8]}_{cmap}_PartType{part_type}_res{res}.png",
        dpi=res // 8,
    )
