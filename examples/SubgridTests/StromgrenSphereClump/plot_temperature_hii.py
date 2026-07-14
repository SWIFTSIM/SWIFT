################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
"""
Plot a mass-weighted temperature slice (through the star's z-plane) and
projection (integrated along z) of a snapshot, with the star and its
HIIRegionRadii overlaid as a circle.

Usage:
    python3 plot_temperature_hii.py                    # last snapshot in snap/
    python3 plot_temperature_hii.py -s snap/snapshot -i 20
"""

import argparse
import glob

import matplotlib.pyplot as plt
import numpy as np
import swiftsimio as sw
import unyt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from swiftsimio.visualisation.projection import project_gas
from swiftsimio.visualisation.slice import slice_gas
from unyt import mh

image_resolution = 1024
cmap = "inferno"


def get_gas_mu(data: sw.SWIFTDataset) -> np.array:
    """Mean molecular weight (Grackle0 2-state mu: ionized vs. neutral,
    from the internal-energy threshold)."""
    from unyt.physical_constants import kboltz_cgs as k_B_cgs

    mh.convert_to_cgs()
    gamma = data.metadata.gas_gamma[0]

    u = data.gas.internal_energies.to_physical().to(unyt.erg / unyt.g)
    H_frac = float(data.metadata.parameters["GrackleCooling:HydrogenFractionByMass"])

    T_over_mu = (gamma - 1.0) * u.value * mh.value / k_B_cgs.value
    T_trans = 1.1e4
    mu_trans = 4.0 / (8.0 - 5.0 * (1.0 - H_frac))

    mu = np.ones(np.size(u))
    mu[T_over_mu > (T_trans + 1) / mu_trans] = 4.0 / (8.0 - 5.0 * (1.0 - H_frac))
    mu[T_over_mu < (T_trans + 1) / mu_trans] = 4.0 / (1.0 + 3.0 * H_frac)
    return mu


def get_gas_temperatures(data: sw.SWIFTDataset) -> np.array:
    from unyt.physical_constants import kboltz_cgs as k_B

    mh.convert_to_cgs()
    u = data.gas.internal_energies.to_physical().to(unyt.erg / unyt.g)
    gamma = data.metadata.gas_gamma[0]
    mu = get_gas_mu(data)

    a = (gamma - 1.0) * (mu * mh) / k_B * u
    T = np.where((u.value > 0), a.value, 0) * unyt.kelvin
    return T


def find_latest_snapshot(stub):
    files = sorted(glob.glob(f"{stub}_*.hdf5"))
    if not files:
        raise RuntimeError(f"No snapshots found matching {stub}_*.hdf5")
    return files[-1]


def get_star_position_and_hii_radius(filename, stub):
    """Star position (kpc) in `filename`, and the largest HIIRegionRadii (kpc)
    seen for that star up to and including `filename` -- HIIRegionRadii
    resets to 0 on any step that isn't an HII rebuild step, so the last
    snapshot alone usually reads 0 even after the region has grown.
    Returns (None, None) if there are no stars."""
    data = sw.load(filename)
    star_pos = np.array(data.stars.coordinates.to("kpc").value)
    if star_pos.size == 0:
        return None, None

    target_index = int(filename.split("_")[-1].split(".")[0])
    peak_r_hii = 0.0
    for f in sorted(glob.glob(f"{stub}_*.hdf5")):
        index = int(f.split("_")[-1].split(".")[0])
        if index > target_index:
            continue
        r_hii = np.max(sw.load(f).stars.hiiregion_radii.to("kpc").value)
        peak_r_hii = max(peak_r_hii, float(r_hii))

    return star_pos[0], peak_r_hii


def draw_star_and_hii_circle(ax, star_pos, r_hii):
    if star_pos is None:
        return
    ax.scatter([star_pos[0]], [star_pos[1]], c="limegreen", zorder=3, marker="*", s=60)
    if r_hii is not None and r_hii > 0:
        circle = Circle(
            (star_pos[0], star_pos[1]),
            r_hii,
            fill=False,
            edgecolor="cyan",
            linewidth=1.2,
            linestyle="--",
            zorder=3,
        )
        ax.add_patch(circle)
        ax.text(
            star_pos[0],
            star_pos[1] + r_hii * 1.05,
            f"$r_{{\\rm HII}}$ = {r_hii * 1e3:.2f} pc",
            color="cyan",
            fontsize=6,
            ha="center",
            va="bottom",
        )


def make_temperature_projection(data, image_resolution):
    data.gas.temperatures = get_gas_temperatures(data)
    data.gas.mass_weighted_temperature = data.gas.masses * data.gas.temperatures

    numerator = project_gas(
        data,
        resolution=image_resolution,
        project="mass_weighted_temperature",
        parallel=True,
        periodic=True,
    )
    denominator = project_gas(
        data,
        resolution=image_resolution,
        project="masses",
        parallel=True,
        periodic=True,
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        T_map = (numerator / denominator).to(unyt.K)
    T_map = np.where(denominator.value > 0, T_map.value, np.nan)

    boxsize = data.metadata.boxsize
    x_edges = np.linspace(0 * unyt.kpc, boxsize[0], image_resolution)
    y_edges = np.linspace(0 * unyt.kpc, boxsize[1], image_resolution)
    return T_map.T, x_edges, y_edges


def make_temperature_slice(data, image_resolution, z_slice):
    data.gas.temperatures = get_gas_temperatures(data)
    data.gas.mass_weighted_temperature = data.gas.masses * data.gas.temperatures

    numerator = slice_gas(
        data,
        resolution=image_resolution,
        z_slice=z_slice,
        project="mass_weighted_temperature",
        parallel=True,
        periodic=True,
    )
    denominator = slice_gas(
        data,
        resolution=image_resolution,
        z_slice=z_slice,
        project="masses",
        parallel=True,
        periodic=True,
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        T_map = (numerator / denominator).to(unyt.K)
    T_map = np.where(denominator.value > 0, T_map.value, np.nan)

    boxsize = data.metadata.boxsize
    x_edges = np.linspace(0 * unyt.kpc, boxsize[0], image_resolution)
    y_edges = np.linspace(0 * unyt.kpc, boxsize[1], image_resolution)
    return T_map.T, x_edges, y_edges


def plot_map(T_map, x_edges, y_edges, star_pos, r_hii, title, out_name):
    fig, ax = plt.subplots(1, figsize=(6, 5), dpi=300)
    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")
    ax.set_title(title, fontsize=9)

    mappable = ax.pcolormesh(
        x_edges, y_edges, T_map, cmap=cmap, norm=LogNorm(vmin=10, vmax=2e4)
    )
    fig.colorbar(mappable, label="Temperature [K]", pad=0)

    draw_star_and_hii_circle(ax, star_pos, r_hii)

    fig.tight_layout()
    fig.savefig(out_name)
    print(f"{out_name} saved.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--stub",
        type=str,
        default="snap/snapshot",
        help="Root of the snapshot filenames, without the index/extension "
        "(default: snap/snapshot).",
    )
    parser.add_argument(
        "-i",
        "--index",
        type=int,
        default=None,
        help="Snapshot index to plot. Default: the last snapshot found.",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        default=image_resolution,
        help=f"Image resolution (default: {image_resolution}).",
    )
    args = parser.parse_args()

    if args.index is None:
        filename = find_latest_snapshot(args.stub)
    else:
        filename = f"{args.stub}_{args.index:04d}.hdf5"

    print(f"Plotting {filename}")
    star_pos, r_hii = get_star_position_and_hii_radius(filename, args.stub)

    t_myr = float(sw.load(filename).metadata.time.to("Myr").value)
    suffix = filename.split("/")[-1].replace(".hdf5", "")

    # Projection (integrated along z).
    data_proj = sw.load(filename)
    T_proj, x_edges, y_edges = make_temperature_projection(data_proj, args.resolution)
    plot_map(
        T_proj,
        x_edges,
        y_edges,
        star_pos,
        r_hii,
        title=f"Mass-weighted temperature projection, t={t_myr:.2f} Myr",
        out_name=f"temperature_projection_{suffix}.png",
    )

    # Slice through the star's z-plane (falls back to the box center if
    # there is no star).
    data_slice = sw.load(filename)
    z_slice = (
        star_pos[2] * unyt.kpc
        if star_pos is not None
        else 0.5 * data_slice.metadata.boxsize[2]
    )
    T_slice, x_edges, y_edges = make_temperature_slice(
        data_slice, args.resolution, z_slice
    )
    plot_map(
        T_slice,
        x_edges,
        y_edges,
        star_pos,
        r_hii,
        title=f"Mass-weighted temperature slice (z={float(z_slice.to('kpc').value if hasattr(z_slice, 'to') else z_slice):.4f} kpc), "
        f"t={t_myr:.2f} Myr",
        out_name=f"temperature_slice_{suffix}.png",
    )


if __name__ == "__main__":
    main()
