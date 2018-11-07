"""
Makes a movie of the KH 2D data.

You will need to run your movie with far higher time-resolution than usual to
get a nice movie; around 450 snapshots over 6s is required.

Edit this file near the bottom with the number of snaps you have.

Written by Josh Borrow (joshua.borrow@durham.ac.uk)
"""

import os
import h5py as h5
import numpy as np
import scipy.interpolate as si


def load_and_extract(filename):
    """
    Load the data and extract relevant info.
    """

    with h5.File(filename, "r") as f:
        x, y, _ = f["PartType0/Coordinates"][...].T
        density = f["PartType0/Density"][...]

    return x, y, density


def make_plot(filename, array, nx, ny, dx, dy):
    """
    Load the data and plop it on the grid using nearest
    neighbour searching for finding the 'correct' value of
    the density.
    """

    data_x, data_y, density = load_and_extract(filename)

    # Make the grid
    x = np.linspace(*dx, nx)
    y = np.linspace(*dy, ny)

    xv, yv = np.meshgrid(x, y)

    mesh = si.griddata((data_x, data_y), density, (xv, yv), method="nearest")

    array.set_array(mesh)

    return array,


def frame(n, *args):
    """
    Make a single frame. Requires the global variables plot and dpi.
    """

    global plot, dpi

    fn = "{}_{:04d}.hdf5".format(filename, n)

    return make_plot(fn, plot, dpi, dpi, (0, 1), (0, 1))


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")

    from tqdm import tqdm
    from matplotlib.animation import FuncAnimation
    from scipy.stats import gaussian_kde

    import matplotlib.pyplot as plt

    filename = "kelvinhelmholtz"
    dpi = 512


    # Look for the number of files in the directory.
    i = 0
    while True:
        if os.path.isfile("{}_{:04d}.hdf5".format(filename, i)):
            i += 1
        else:
            break

        if i > 10000:
            raise FileNotFoundError(
                "Could not find the snapshots in the directory")

    frames = tqdm(np.arange(0, i))

    # Creation of first frame
    fig, ax = plt.subplots(1, 1, figsize=(1, 1), frameon=False)
    ax.axis("off")  # Remove annoying black frame.

    data_x, data_y, density = load_and_extract("kelvinhelmholtz_0000.hdf5")

    x = np.linspace(0, 1, dpi)
    y = np.linspace(0, 1, dpi)
    xv, yv = np.meshgrid(x, y)

    mesh = si.griddata((data_x, data_y), density, (xv, yv), method="nearest")
    
    # Global variable for set_array
    plot = ax.imshow(mesh, extent=[0, 1, 0, 1], animated=True, interpolation="none")

    anim = FuncAnimation(fig, frame, frames, interval=40, blit=False)

    # Remove all whitespace
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)

    # Actually make the movie
    anim.save("khmovie.mp4", dpi=dpi, bitrate=4096)
