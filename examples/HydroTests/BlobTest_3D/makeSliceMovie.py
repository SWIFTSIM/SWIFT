"""
Makes a movie of the output of the blob test.

Josh Borrow (joshua.borrow@durham.ac.uk) 2019

LGPLv3
"""

from swiftsimio import load
from swiftsimio.visualisation import slice

from p_tqdm import p_map

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation

info_frames = 15
start_frame = 0
end_frame = 101
resolution = 1024
snapshot_name = "blob"
cmap = "Spectral_r"
text_args = dict(color="black")
# plot = "pressure"
# name = "Pressure $P$"
plot = "density"
name = "Fluid Density $\\rho$"


def get_image(n):
    """
    Gets the image for snapshot n, and also returns the associated
    SWIFT metadata object.
    """
    filename = f"{snapshot_name}_{n:04d}.hdf5"

    data = load(filename)
    boxsize = data.metadata.boxsize[0].value

    output = np.zeros((resolution, resolution * 4), dtype=float)

    x, y, z = data.gas.coordinates.value.T

    # This is an oblong box but we can only make squares!
    for box, box_edges in enumerate([[0.0, 1.1], [0.9, 2.1], [1.9, 3.1], [2.9, 4.0]]):
        mask = np.logical_and(x >= box_edges[0], x <= box_edges[1])
        masked_x = x[mask] - np.float64(box)
        masked_y = y[mask]
        masked_z = z[mask]

        try:
            hsml = data.gas.smoothing_length.value[mask]
        except:
            hsml = data.gas.smoothing_lengths.value[mask]

        if plot == "density":
            mass = data.gas.masses.value[mask]
            image = slice(
                x=masked_y,
                y=masked_x,
                z=masked_z,
                m=mass,
                h=hsml,
                z_slice=0.5,
                res=resolution,
            )
        else:
            quantity = getattr(data.gas, plot).value[mask]
            # Need to divide out the particle density for non-projected density quantities
            image = scatter(
                x=masked_y,
                y=masked_x,
                z=masked_z,
                m=quantity,
                h=hsml,
                z_slice=0.5,
                res=resolution,
            ) / scatter(
                x=masked_y,
                y=masked_x,
                z=masked_z,
                m=np.ones_like(quantity),
                h=hsml,
                z_slice=0.5,
                res=resolution,
            )

        output[:, box * resolution : (box + 1) * resolution] = image

    return output, data.metadata


def get_data_dump(metadata):
    """
    Gets a big data dump from the SWIFT metadata
    """

    try:
        viscosity = metadata.viscosity_info
    except:
        viscosity = "No info"

    try:
        diffusion = metadata.diffusion_info
    except:
        diffusion = "No info"

    output = (
        "$\\bf{Blob}$ $\\bf{Test}$\n\n"
        "$\\bf{SWIFT}$\n"
        + metadata.code_info
        + "\n\n"
        + "$\\bf{Compiler}$\n"
        + metadata.compiler_info
        + "\n\n"
        + "$\\bf{Hydrodynamics}$\n"
        + metadata.hydro_info
        + "\n\n"
        + "$\\bf{Viscosity}$\n"
        + viscosity
        + "\n\n"
        + "$\\bf{Diffusion}$\n"
        + diffusion
    )

    return output


def time_formatter(metadata):
    return f"$t = {metadata.t:2.2f}$"


# Generate the frames and unpack our variables
images, metadata = zip(*p_map(get_image, list(range(start_frame, end_frame))))

# The edges are funny because of the non-periodicity.
central_region = images[0][
    resolution // 10 : resolution - resolution // 10,
    resolution // 10 : resolution - resolution // 10,
]
norm = LogNorm(vmin=np.min(central_region), vmax=np.max(central_region), clip="black")

fig, ax = plt.subplots(figsize=(8 * 4, 8), dpi=resolution // 8)

fig.subplots_adjust(0, 0, 1, 1)
ax.axis("off")

# Set up the initial state
image = ax.imshow(np.zeros_like(images[0]), norm=norm, cmap=cmap, origin="lower")

description_text = ax.text(
    0.5,
    0.5,
    get_data_dump(metadata[0]),
    va="center",
    ha="center",
    **text_args,
    transform=ax.transAxes,
)

time_text = ax.text(
    (1 - 0.025 * 0.25),
    0.975,
    time_formatter(metadata[0]),
    **text_args,
    va="top",
    ha="right",
    transform=ax.transAxes,
)

info_text = ax.text(
    0.025 * 0.25, 0.975, name, **text_args, va="top", ha="left", transform=ax.transAxes
)


def animate(n):
    # Display just our original frames at t < 0
    if n >= 0:
        image.set_array(images[n])
        description_text.set_text("")
        time_text.set_text(time_formatter(metadata[n]))

    return (image,)


animation = FuncAnimation(
    fig, animate, range(start_frame - info_frames, end_frame), interval=40
)

animation.save(filename="blob_slice.mp4")
