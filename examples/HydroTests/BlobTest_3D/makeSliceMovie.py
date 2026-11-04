"""
Makes a movie of the output of the blob test.

Josh Borrow (joshua.borrow@durham.ac.uk) 2019

LGPLv3
"""

from swiftsimio import load, cosmo_array, cosmo_quantity
from swiftsimio.visualisation import slice_gas

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
    simulation time.
    """
    filename = f"{snapshot_name}_{n:04d}.hdf5"

    data = load(filename)
    coords = data.gas.coordinates

    output = np.zeros((resolution, resolution * 4), dtype=float)
    z_slice = cosmo_quantity(
        0.5,
        coords.units,
        comoving=coords.comoving,
        scale_factor=data.metadata.scale_factor,
        scale_exponent=1,
    )

    # This is an oblong box but we can only make squares!
    for box in range(4):
        region = cosmo_array(
            [box, box + 1.0, 0.0, 1.0],
            coords.units,
            comoving=coords.comoving,
            scale_factor=data.metadata.scale_factor,
            scale_exponent=1,
        )

        if plot == "density":
            # Density in the slice
            image = slice_gas(data, resolution, z_slice, "masses", region=region)
            image = image.to_value("g/cm**3")
        else:
            # Need to divide out the particle density for non-projected density quantities
            image = slice_gas(
                data, resolution, z_slice, plot, region=region
            ) / slice_gas(data, resolution, z_slice, None, region=region)
            image = image.to_value(getattr(data.gas, plot).units)

        # The images are indexed [x, y]; imshow wants [y, x]
        output[:, box * resolution : (box + 1) * resolution] = image.T

    return output, data.metadata.t.to_value("s")


def get_data_dump(metadata):
    """
    Gets a big data dump from the SWIFT metadata
    """

    try:
        viscosity = metadata.viscosity_info
    except KeyError:
        viscosity = "No info"

    try:
        diffusion = metadata.diffusion_info
    except KeyError:
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


def time_formatter(t):
    return f"$t = {t:2.2f}$"


# Generate the frames and unpack our variables
images, times = zip(*p_map(get_image, list(range(start_frame, end_frame))))

# Set the colour range from the central part of the first frame.
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
    get_data_dump(load(f"{snapshot_name}_{start_frame:04d}.hdf5").metadata),
    va="center",
    ha="center",
    **text_args,
    transform=ax.transAxes,
)

time_text = ax.text(
    (1 - 0.025 * 0.25),
    0.975,
    time_formatter(times[0]),
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
        time_text.set_text(time_formatter(times[n]))

    return (image,)


animation = FuncAnimation(
    fig, animate, range(start_frame - info_frames, end_frame), interval=40
)

animation.save(filename="blob_slice.mp4")
