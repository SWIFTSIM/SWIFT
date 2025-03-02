"""
Makes a gas density projection plot. Uses the swiftsimio library.
"""

import matplotlib.pyplot as plt
import numpy as np

import swiftsimio as sw
from swiftsimio.visualisation.projection import project_gas

import unyt
from unyt import mh, cm
from tqdm import tqdm
from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation

# %%
# Constants; these could be put in the parameter file but are rarely changed.
image_resolution = 1024

# Plotting controls
cmap = "inferno"

# %%


def get_gas_mu(data: sw.SWIFTDataset) -> np.array:
    """
    Return the mean molecular weight of the gas.
    """
    # Get the cooling model
    cooling = data.metadata.subgrid_scheme["Cooling Model"]

    # Use a variable for better readability
    gas = data.gas

    # Get rho
    rho = gas.densities.to(unyt.g / unyt.cm ** 3)  # self.Rho(units='g/cm3')

    # hydrogen mass in gram
    mh.convert_to_cgs()
    mH_in_g = mh

    if cooling == b"Grackle3":  # grackle3
        nHI = gas.hi * rho / (mH_in_g)
        nHII = gas.hii * rho / (mH_in_g)
        nHeI = gas.he_i * rho / (4 * mH_in_g)
        nHeII = gas.he_ii * rho / (4 * mH_in_g)
        nHeIII = gas.he_iii * rho / (4 * mH_in_g)
        nH2I = gas.h2_i * rho / (2 * mH_in_g)
        nH2II = gas.h2_ii * rho / (2 * mH_in_g)
        nHDI = gas.hdi * rho / (3 * mH_in_g)

        nel = nHII + nHeII + 2 * nHeIII + nH2II
        mu = (
            (nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 + (nH2I + nH2II) * 2 + nHDI * 3
        ) / (nHI + nHII + nHeI + nHeII + nHeIII + nH2I + nH2II + nHDI + nel)
        return mu

    elif cooling == b"Grackle2":  # grackle2
        nHI = gas.hi * rho / (mH_in_g)
        nHII = gas.hii * rho / (mH_in_g)
        nHeI = gas.he_i * rho / (4 * mH_in_g)
        nHeII = gas.he_ii * rho / (4 * mH_in_g)
        nHeIII = gas.he_iii * rho / (4 * mH_in_g)
        nH2I = gas.h2_i * rho / (2 * mH_in_g)
        nH2II = gas.h2_ii * rho / (2 * mH_in_g)

        nel = nHII + nHeII + 2 * nHeIII + nH2II
        mu = ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 + (nH2I + nH2II) * 2) / (
            nHI + nHII + nHeI + nHeII + nHeIII + nH2I + nH2II + nel
        )
        return mu

    elif cooling == b"Grackle1":  # grackle1
        nHI = gas.hi * rho / (mH_in_g)
        nHII = gas.hii * rho / (mH_in_g)
        nHeI = gas.he_i * rho / (4 * mH_in_g)
        nHeII = gas.he_ii * rho / (4 * mH_in_g)
        nHeIII = gas.he_iii * rho / (4 * mH_in_g)
        nel = nHII + nHeII + 2 * nHeIII
        mu = ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4) / (
            nHI + nHII + nHeI + nHeII + nHeIII + nel
        )
        return mu

    else:  # Grackle0
        # print("... found info for grackle mode=0")

        from unyt.physical_constants import kboltz_cgs as k_B_cgs

        gamma = data.metadata.gas_gamma[0]

        # Get internal energy
        u = data.gas.internal_energies
        u = u.to_physical()
        u = u.to(unyt.erg / unyt.g)

        # Get hydrigen fraction
        H_frac = float(
            data.metadata.parameters["GrackleCooling:HydrogenFractionByMass"]
        )

        # Compute T/mu
        T_over_mu = (gamma - 1.0) * u.value * mH_in_g.value / k_B_cgs.value
        T_trans = 1.1e4
        mu_trans = 4.0 / (8.0 - 5.0 * (1.0 - H_frac))

        # Determine if we are ionized or not
        mu = np.ones(np.size(u))
        mask_ionized = T_over_mu > (T_trans + 1) / mu_trans
        mask_neutral = T_over_mu < (T_trans + 1) / mu_trans

        # Give the right mu
        mu[mask_ionized] = 4.0 / (8.0 - 5.0 * (1.0 - H_frac))
        mu[mask_neutral] = 4.0 / (1.0 + 3.0 * H_frac)

        return mu


def get_gas_temperatures(data: sw.SWIFTDataset) -> np.array:
    """
        Compute the temperature of the gas.
    """
    from unyt.physical_constants import kboltz_cgs as k_B

    # from unyt.physical_constants import mh

    # Convert to cgs
    mh.convert_to_cgs()

    # Get the cooling model
    cooling = data.metadata.subgrid_scheme["Cooling Model"]

    # Get internal energy and convert to physical units in cgs
    u = data.gas.internal_energies
    u = u.to_physical()
    u = u.to(unyt.erg / unyt.g)

    # Get gamm and compute mu
    gamma = data.metadata.gas_gamma[0]
    mu = get_gas_mu(data)

    # FInally compute the the temperature
    if cooling == b"Grackle3" or cooling == b"Grackle2" or cooling == b"Grackle1":
        T = mu * (gamma - 1.0) * u * mh / k_B
    else:
        a = (gamma - 1.0) * (mu * mh) / k_B * u
        T = np.where((u.value > 0), a.value, 0) * unyt.kelvin
    return T


def get_sink_and_stars_positions(filename):
    data = sw.load(filename)

    sink_pos = data.sinks.coordinates
    star_pos = data.stars.coordinates

    return sink_pos, star_pos


def make_projection(filename, image_resolution):
    """
    Compute a mass projection with swiftsimio.
    """
    data = sw.load(filename)

    # Compute projected density
    projected_mass = project_gas(
        data,
        resolution=image_resolution,
        project="masses",
        parallel=True,
        periodic=True,
    )

    boxsize = data.metadata.boxsize
    x_edges = np.linspace(0 * unyt.kpc, boxsize[0], image_resolution)
    y_edges = np.linspace(0 * unyt.kpc, boxsize[1], image_resolution)

    # Convert to 1/cm**2
    projected_mass = projected_mass.to(mh / (cm ** 2))

    return projected_mass.T, x_edges, y_edges


def setup_axes():
    """
    Creates the figure and axis object.
    """
    fig, ax = plt.subplots(1, figsize=(6, 5), dpi=300)

    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")

    return fig, ax


def make_single_image(filename, image_resolution):
    """
    Makes a single image and saves it to mass_projection_{snapshot_number}.png.

    Filename should be given _without_ hdf5 extension.
    """
    file = "{:s}.hdf5".format(filename)
    fig, ax = setup_axes()
    projected_mass, x_edges, y_edges = make_projection(file, image_resolution)

    mappable = ax.pcolormesh(
        x_edges, y_edges, projected_mass, cmap=cmap, norm=LogNorm()
    )
    fig.colorbar(mappable, label="Surface density [cm$^{-2}]$", pad=0)

    sink_pos, star_pos = get_sink_and_stars_positions(file)

    if star_pos.size != 0:
        ax.scatter(star_pos[:, 0], star_pos[:, 1], c="limegreen", zorder=1, marker="*")

    if sink_pos.size != 0:
        ax.scatter(sink_pos[:, 0], sink_pos[:, 1], c="blue", zorder=2, marker=".")

    ax.text(
        0.7,
        0.95,
        "$N_{\mathrm{star}}" + " = {}$".format(len(star_pos)),
        transform=ax.transAxes,
        fontsize=7,
        bbox=dict(facecolor="white", alpha=0.8),
    )
    ax.text(
        0.1,
        0.95,
        "$N_{\mathrm{sink}}" + " = {}$".format(len(sink_pos)),
        transform=ax.transAxes,
        fontsize=7,
        bbox=dict(facecolor="white", alpha=0.8),
    )

    fig.tight_layout()

    image = "mass_projection_{:s}.png".format(filename[-4:])
    fig.savefig(image)

    return


def make_movie(args, image_resolution):
    """
    Makes a movie and saves it to mass_projection_movie.mp4.
    """

    fig, ax = setup_axes()

    def grab_metadata(n):
        filename = "{:s}_{:04d}.hdf5".format(args["stub"], n)
        data = sw.load(filename)

        return data.metadata

    def grab_data(n):
        filename = "{:s}_{:04d}.hdf5".format(args["stub"], n)

        H, _, _ = make_projection(filename, image_resolution)

        # Need to ravel because pcolormesh's set_array takes a 1D array. Might
        # as well do it here, because 1d arrays are easier to max() than 2d.
        return H.ravel()

    def grab_sink_star_pos(n):
        filename = "{:s}_{:04d}.hdf5".format(args["stub"], n)

        sink_pos, star_pos = get_sink_and_stars_positions(filename)

        return sink_pos, star_pos

    histograms = [
        grab_data(n)
        for n in tqdm(
            range(args["initial"], args["final"] + 1),
            desc="Computing gas mass projection",
        )
    ]

    sink_star = [
        grab_sink_star_pos(n)
        for n in tqdm(
            range(args["initial"], args["final"] + 1),
            desc="Getting sink and stars positions",
        )
    ]

    sink_pos = [s[0] for s in sink_star]
    star_pos = [s[1] for s in sink_star]

    metadata = [
        grab_metadata(n)
        for n in tqdm(
            range(args["initial"], args["final"] + 1), desc="Grabbing metadata"
        )
    ]

    # Need to get a reasonable norm so that we don't overshoot.
    max_particles = max([x.max() for x in histograms])
    min_particles = max([x.min() for x in histograms])
    norm = LogNorm(vmin=min_particles, vmax=max_particles)

    # First, let's make the initial frame (we need this for our d, T values that we
    # got rid of in grab_data.
    hist, x_edges, y_edges = make_projection(
        "{:s}_{:04d}.hdf5".format(args["stub"], args["initial"]), image_resolution
    )

    mappable = ax.pcolormesh(x_edges, y_edges, hist, cmap=cmap, norm=norm)
    fig.colorbar(mappable, label="Surface density [cm$^{-2}]$", pad=0)

    fig.tight_layout()

    # Once we've rearranged the figure with tight_layout(), we can start laying
    # Down the metadata text.

    def format_metadata(metadata: sw.SWIFTMetadata):
        t = metadata.t
        t.convert_to_units(unyt.Myr)

        x = "$a$: {:2.2f}\n$z$: {:2.2f}\n$t$: {:2.2f}".format(metadata.a, metadata.z, t)

        return x

    text = ax.text(
        0.025,
        0.975,
        format_metadata(metadata[0]),
        ha="left",
        va="top",
        transform=ax.transAxes,
        color="white",
    )

    ax.text(
        0.975,
        0.975,
        metadata[0].code["Git Revision"].decode("utf-8"),
        ha="right",
        va="top",
        transform=ax.transAxes,
        color="white",
    )

    def animate(data):
        mappable.set_array(histograms[data])
        text.set_text(format_metadata(metadata[data]))

        if star_pos[data].size != 0:
            ax.scatter(
                star_pos[data][:, 0],
                star_pos[data][:, 1],
                c="limegreen",
                zorder=1,
                marker="*",
            )

        if sink_pos[data].size != 0:
            ax.scatter(
                sink_pos[data][:, 0],
                sink_pos[data][:, 1],
                c="blue",
                zorder=2,
                marker=".",
            )

        return mappable

    animation = FuncAnimation(
        fig, animate, range(len(histograms)), fargs=[], interval=1000 / 10
    )

    animation.save("mass_projection_movie.mp4")

    return


# %%
if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(
        description="""
             Plotting script for making a mass projection plot.
             Takes the filename handle, start, and (optionally) stop
             snapshots. If stop is not given, png plot is produced for
             that snapshot. If given, a movie is made.
             """
    )

    parser.add_argument(
        "-i",
        "--initial",
        help="""Initial snapshot number. Default: 0.""",
        default=0,
        required=False,
        type=int,
    )

    parser.add_argument(
        "-f",
        "--final",
        help="""Final snapshot number. Default: 0.""",
        default=0,
        required=False,
        type=int,
    )

    parser.add_argument(
        "-s",
        "--stub",
        help="""Root of the snapshots filenames (e.g. snapshot). This is
                the first part of the filename for the snapshots,
                not including the final underscore. Required.""",
        type=str,
        required=True,
    )

    args = vars(parser.parse_args())

    if args["final"] <= args["initial"]:
        # Run in single image mode.
        filename = "{:s}_{:04d}".format(args["stub"], args["initial"])

        make_single_image(filename, image_resolution=image_resolution)

    else:
        # Movie mode!
        make_movie(args, image_resolution=image_resolution)
