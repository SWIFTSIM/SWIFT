#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
##############################################################################


# ----------------------------------------------------
# Plot slices of the hydrogen number density,
# pressure, temperature, and hydrogen mass fraction.
# ----------------------------------------------------

import sys
import swiftsimio
import gc
import unyt
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from swiftsimio.visualisation.slice import slice_gas
import stromgren_plotting_tools as spt

# Parameters users should/may tweak

# snapshot basename
snapshot_base = "output"

# parameters for imshow plots
imshow_kwargs = {"origin": "lower"}

# parameters for swiftsimio slices
slice_kwargs = {"slice": 0.5, "resolution": 1000, "parallel": True}

# -----------------------------------------------------------------------


# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True
    snapnr = -1

mpl.rcParams["text.usetex"] = True


def set_colorbar(ax, im):
    """
    Adapt the colorbar a bit for axis object <ax> and
    imshow instance <im>
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    return


def plot_result(filename):
    """
    Create and save the plot
    """
    print("working on", filename)

    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    global imshow_kwargs
    imshow_kwargs["extent"] = [
        0.0 * meta.boxsize[0].v,
        0.9 * meta.boxsize[0].v,
        0.0 * meta.boxsize[1].v,
        0.9 * meta.boxsize[1].v,
    ]
    cutoff = int(0.05 * slice_kwargs["resolution"])

    mass_map = slice_gas(data, project="masses", **slice_kwargs)
    gamma = meta.hydro_scheme["Adiabatic index"][0]

    imf = spt.get_imf(scheme, data)

    data.gas.mXHI = imf.HI * data.gas.masses
    data.gas.mXHII = imf.HII * data.gas.masses
    data.gas.mP = data.gas.pressures * data.gas.masses
    data.gas.mrho = data.gas.densities * data.gas.masses

    mu = spt.mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
    data.gas.mT = (
        spt.gas_temperature(data.gas.internal_energies, mu, gamma) * data.gas.masses
    )

    mass_weighted_HI_map = slice_gas(data, project="mXHI", **slice_kwargs)
    mass_weighted_pressure_map = slice_gas(data, project="mP", **slice_kwargs)
    mass_weighted_density_map = slice_gas(data, project="mrho", **slice_kwargs)
    mass_weighted_temperature_map = slice_gas(data, project="mT", **slice_kwargs)

    HI_map = mass_weighted_HI_map / mass_map
    HI_map = HI_map[cutoff:-cutoff, cutoff:-cutoff]

    pressure_map = mass_weighted_pressure_map / mass_map
    pressure_map = pressure_map[cutoff:-cutoff, cutoff:-cutoff]
    pressure_map = pressure_map.to("g/cm/s**2")

    density_map = mass_weighted_density_map / mass_map
    density_map = density_map[cutoff:-cutoff, cutoff:-cutoff]
    density_map = density_map.to("kg/cm**3")
    density_map = density_map / unyt.proton_mass

    temperature_map = mass_weighted_temperature_map / mass_map
    temperature_map = temperature_map[cutoff:-cutoff, cutoff:-cutoff]
    temperature_map = temperature_map.to("K")

    fig = plt.figure(figsize=(12, 12), dpi=200)
    figname = filename[:-5] + ".png"

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    try:
        im1 = ax1.imshow(
            density_map.T,
            **imshow_kwargs,
            norm=LogNorm(vmin=1e-4, vmax=1e-1),
            cmap="bone",
        )
        set_colorbar(ax1, im1)
        ax1.set_title(r"Hydrogen Number Density [cm$^{-3}$]")
    except ValueError:
        print(
            filename,
            "densities wrong? min",
            data.gas.densities.min(),
            "max",
            data.gas.densities.max(),
        )
        return

    try:
        im2 = ax2.imshow(
            HI_map.T, **imshow_kwargs, norm=LogNorm(vmin=1e-3, vmax=1.0), cmap="cividis"
        )
        set_colorbar(ax2, im2)
        ax2.set_title("HI Mass Fraction [1]")
    except ValueError:
        print(filename, "mass fraction wrong? min", imf.HI.min(), "max", imf.HI.max())
        return

    try:
        im3 = ax3.imshow(
            pressure_map.T,
            **imshow_kwargs,
            norm=LogNorm(vmin=1e-15, vmax=1e-12),
            cmap="viridis",
        )
        set_colorbar(ax3, im3)
        ax3.set_title(r"Pressure [g/cm/s$^2$]")
    except ValueError:
        print(
            filename,
            "pressures wrong? min",
            data.gas.pressures.min(),
            "max",
            data.gas.pressures.max(),
        )
        return

    try:
        im4 = ax4.imshow(
            temperature_map.T,
            **imshow_kwargs,
            norm=LogNorm(vmin=1e2, vmax=4e4),
            cmap="inferno",
        )
        set_colorbar(ax4, im4)
        ax4.set_title(r"Temperature [K]")
    except ValueError:
        print(
            filename,
            "temperatures wrong? min",
            temperature_map.min(),
            "max",
            temperature_map.max(),
        )
        return

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlabel("[kpc]")
        ax.set_ylabel("[kpc]")

    title = filename.replace("_", "\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2e}".format(meta.time.to("Myr"))
    fig.suptitle(title)

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    gc.collect()
    return


if __name__ == "__main__":
    snaplist = spt.get_snapshot_list(snapshot_base, plot_all, snapnr)
    for f in snaplist:
        plot_result(f)
