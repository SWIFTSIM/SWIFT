"""
Makes a movie of the output of the blob test.

Josh Borrow (joshua.borrow@durham.ac.uk) 2019

LGPLv3
"""

from swiftsimio import load, cosmo_array
from swiftsimio.visualisation import project_gas

from p_tqdm import p_map
from types import SimpleNamespace

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
# plot = "pressures"
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
    coords = data.gas.coordinates

    output = np.zeros((resolution, resolution * 4), dtype=float)

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
            image = project_gas(data, resolution, "masses", region=region)
        else:
            # Need to divide out the particle density for non-projected density quantities
            image = project_gas(data, resolution, plot, region=region) / project_gas(
                data, resolution, None, region=region
            )

        # The images are indexed [x, y]; imshow wants [y, x]
        output[:, box * resolution : (box + 1) * resolution] = image.to_value().T

    return output, picklable_metadata(data.metadata)


def picklable_metadata(metadata):
    """
    The swiftsimio metadata object holds an open HDF5 file handle and cannot
    be sent back from the p_map worker processes, so keep only what we use.
    """

    def info(name):
        try:
            return getattr(metadata, name)
        except Exception:
            return "No info"

    return SimpleNamespace(
        t=metadata.t,
        code_info=metadata.code_info,
        compiler_info=metadata.compiler_info,
        hydro_info=metadata.hydro_info,
        viscosity_info=info("viscosity_info"),
        diffusion_info=info("diffusion_info"),
    )


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

animation.save(filename="blob.mp4")
