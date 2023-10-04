"""
Makes a movie of the Sedov 2D data. Adapted from

    KelvinHelmholtz_2D/makeMovie.py

You will need to run your movie with far higher time-resolution than usual to
get a nice movie; around 450 snapshots over 6s is required.

Edit this file near the bottom with the number of snaps you have.

Written by Josh Borrow (joshua.borrow@durham.ac.uk)
"""

import numpy as np

if __name__ == "__main__":
    import matplotlib

    matplotlib.use("Agg")
    params = {
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "font.size": 12,
        "legend.fontsize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "text.usetex": True,
        "figure.figsize": (3.15, 2.60),
        "figure.subplot.left": 0.17,
        "figure.subplot.right": 0.99,
        "figure.subplot.bottom": 0.08,
        "figure.subplot.top": 0.99,
        "figure.subplot.wspace": 0.0,
        "figure.subplot.hspace": 0.0,
        "lines.markersize": 6,
        "lines.linewidth": 3.0,
    }
    matplotlib.rcParams.update(params)

    from matplotlib.colors import LogNorm

    import matplotlib.pyplot as plt

    filename = "sedov"
    dpi = 1024

    # Creation of first frame
    fig, ax = plt.subplots(1, 1, frameon=False)

    with np.load("sedov_soundspeed_ratio_data.npz") as file:
        mesh = file.items()[0][1]

    # Global variable for set_array
    img = ax.imshow(
        mesh, extent=[0, 1, 0, 1], animated=True, interpolation="none", norm=LogNorm()
    )

    circle = matplotlib.patches.Circle(
        [0.5, 0.5],
        radius=0.18242863869665918,
        animated=True,
        lw=1,
        fill=False,
        ec="red",
    )

    ax.add_artist(circle)

    fig.colorbar(img, label=r"$c_{s, {\rm smoothed}}$ / $c_{s, {\rm gas}}$", pad=0.0)

    plt.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis="y", which="both", left=False, right=False, labelleft=False)

    plt.xlim(0.2, 0.8)
    plt.ylim(0.2, 0.8)

    # Actually make the movie
    plt.tight_layout()

    plt.savefig("sedov_blast_soundspeed.pdf", dpi=300)
