#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
#               2022 Tsang Keung Chan (chantsangkeung@gmail.com)
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
# plots 2D projection of radiation energy and fluxes
#
# Usage: give snapshot number as cmdline arg to plot
# single snapshot, otherwise this script plots
# all snapshots available in the workdir.
# ----------------------------------------------------

import sys
import os
import swiftsimio
import numpy as np
import gc
import unyt
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import SymLogNorm

# Parameters users should/may tweak

# plot all groups and all photon quantities
plot_all_data = True

# plotting propagation tests or Stromgren sphere?
do_stromgren_sphere = True

# fancy up the plots a bit?
fancy = True

# parameters for imshow plots
imshow_kwargs = {"origin": "lower", "cmap": "viridis"}

# parameters for swiftsimio projections
projection_kwargs = {"resolution": 1024, "parallel": True}

# snapshot basename
snapshot_base = "propagation_test"
if do_stromgren_sphere:
    snapshot_base = "output"

# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True

mpl.rcParams["text.usetex"] = True

# -----------------------------------------------------------------------


# get the unit for rt according to the RT scheme:
def get_units(scheme, unit_system="cgs_units"):
    if unit_system == "cgs_units":
        time_units = unyt.s
        energy_units = unyt.erg
        energy_units_str = "\\rm{erg}"
        if scheme.startswith("GEAR M1closure"):
            flux_units = 1e-10 * energy_units / unyt.cm ** 2 / unyt.s
            flux_units_str = "10^{-10} \\rm{erg} \\ \\rm{cm}^{-2} \\ \\rm{s}^{-1}"
        elif scheme.startswith("SPH M1closure"):
            flux_units = 1e10 * energy_units * unyt.cm / unyt.s
            flux_units_str = "10^{10} \\rm{erg} \\ \\rm{cm} \\ \\rm{s}^{-1}"
        else:
            print("RT scheme not identified. Exit.")
            exit()
    elif unit_system == "stromgren_units":
        time_units = unyt.Myr
        energy_units = 1e50 * unyt.erg
        energy_units_str = "10^{50} \\rm{erg}"
        if scheme.startswith("GEAR M1closure"):
            flux_units = 1e50 * unyt.erg / unyt.kpc ** 2 / unyt.Gyr
            flux_units_str = "10^{60} \\rm{erg} \\ \\rm{kpc}^{-2} \\ \\rm{Gyr}^{-1}"
        elif scheme.startswith("SPH M1closure"):
            flux_units = 1e50 * unyt.erg * unyt.kpc / unyt.Gyr
            flux_units_str = "10^{60} \\rm{erg} \\ \\rm{kpc} \\ \\rm{Gyr}^{-1}"
        else:
            print("RT scheme not identified. Exit.")
            exit()
    else:
        print("Unit system not identified. Exit.")
        exit()
    return time_units, energy_units, energy_units_str, flux_units, flux_units_str


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted 
    and return their names as list
    """

    snaplist = []

    if plot_all:
        dirlist = os.listdir()
        for f in dirlist:
            if f.startswith(snapshot_basename) and f.endswith("hdf5"):
                snaplist.append(f)

        snaplist = sorted(snaplist)

    else:
        fname = snapshot_basename + "_" + str(snapnr).zfill(4) + ".hdf5"
        if not os.path.exists(fname):
            print("Didn't find file", fname)
            quit(1)
        snaplist.append(fname)

    if len(snaplist) == 0:
        print("Didn't find any snapshots with basename '", snapshot_basename, "'")
        quit(1)

    return snaplist


def set_colorbar(ax, im):
    """
    Adapt the colorbar a bit for axis object <ax> and
    imshow instance <im>
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    return


def plot_photons(filename, energy_boundaries=None, flux_boundaries=None):
    """
    Create the actual plot.

    filename: file to work with
    energy_boundaries:  list of [E_min, E_max] for each photon group. 
                        If none, limits are set automatically.
    flux_boundaries:    list of [F_min, F_max] for each photon group. 
                        If none, limits are set automatically.
    """

    print("working on", filename)

    # Read in data first
    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    if do_stromgren_sphere:
        time_units, energy_units, energy_units_str, flux_units, flux_units_str = get_units(
            scheme, unit_system="stromgren_units"
        )
    else:
        time_units, energy_units, energy_units_str, flux_units, flux_units_str = get_units(
            scheme, unit_system="cgs_units"
        )

    ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"])
    xlabel_units_str = meta.boxsize.units.latex_representation()

    global imshow_kwargs
    imshow_kwargs["extent"] = [0, meta.boxsize[0].v, 0, meta.boxsize[1].v]

    rotation_matrix = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]])
    rotation_center = meta.boxsize * 0.5

    for g in range(ngroups):
        # workaround to access named columns data with swiftsimio visualisaiton
        # add mass weights to remove surface density dependence in images
        new_attribute_str = "mass_weighted_radiation_energy" + str(g + 1)
        en = getattr(data.gas.photon_energies, "group" + str(g + 1)) * data.gas.masses
        setattr(data.gas, new_attribute_str, en)

        if plot_all_data:
            # prepare also the fluxes
            for direction in ["X", "Y", "Z"]:
                new_attribute_str = (
                    "mass_weighted_radiation_flux" + str(g + 1) + direction
                )
                f = (
                    getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + direction)
                    * data.gas.masses
                )
                setattr(data.gas, new_attribute_str, f)

    # get mass surface density projection that we'll use to remove density dependence in image
    mass_map = swiftsimio.visualisation.projection.project_gas(
        data, project="masses", **projection_kwargs
    )

    if plot_all_data:
        fig = plt.figure(figsize=(5 * 4, 5.05 * ngroups), dpi=200)
        figname = filename[:-5] + "-radiation-projection.png"

        for g in range(ngroups):

            # get energy projection
            new_attribute_str = "mass_weighted_radiation_energy" + str(g + 1)
            mass_weighted_photon_map = swiftsimio.visualisation.projection.project_gas(
                data, project=new_attribute_str, **projection_kwargs
            )
            # mass_weighted_photon_map /= mass_map may cause RecursionErrors...
            photon_map = mass_weighted_photon_map / mass_map
            photon_map.convert_to_units(energy_units)

            ax = fig.add_subplot(ngroups, 4, g * 4 + 1)

            if energy_boundaries is not None:
                imshow_kwargs["norm"] = SymLogNorm(
                    vmin=energy_boundaries[g][0],
                    vmax=energy_boundaries[g][1],
                    linthresh=1e-3,
                    base=10,
                )
            im = ax.imshow(photon_map.T, **imshow_kwargs)
            set_colorbar(ax, im)
            ax.set_ylabel("Group {0:2d}".format(g + 1))
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            if g == 0:
                ax.set_title("Energies [$" + energy_units_str + "$]")

            # get flux X projection
            new_attribute_str = "mass_weighted_radiation_flux" + str(g + 1) + "X"
            mass_weighted_flux_map = swiftsimio.visualisation.projection.project_gas(
                data, project=new_attribute_str, **projection_kwargs
            )
            # mass_weighted_flux_map /= mass_map may cause RecursionErrors...
            photon_map = mass_weighted_flux_map / mass_map
            photon_map.convert_to_units(flux_units)

            ax = fig.add_subplot(ngroups, 4, g * 4 + 2)
            if flux_boundaries is not None:
                imshow_kwargs["norm"] = SymLogNorm(
                    vmin=flux_boundaries[g][0],
                    vmax=flux_boundaries[g][1],
                    linthresh=1e-3,
                    base=10,
                )
            im = ax.imshow(photon_map.T, **imshow_kwargs)
            set_colorbar(ax, im)
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            ax.set_ylabel("y [$" + xlabel_units_str + "$]")
            if g == 0:
                ax.set_title("Flux X [$" + flux_units_str + "$]")

            # get flux Y projection
            new_attribute_str = "mass_weighted_radiation_flux" + str(g + 1) + "Y"
            mass_weighted_flux_map = swiftsimio.visualisation.projection.project_gas(
                data, project=new_attribute_str, **projection_kwargs
            )
            # mass_weighted_flux_map /= mass_map may cause RecursionErrors...
            photon_map = mass_weighted_flux_map / mass_map
            photon_map.convert_to_units(flux_units)

            ax = fig.add_subplot(ngroups, 4, g * 4 + 3)
            im = ax.imshow(photon_map.T, **imshow_kwargs)
            set_colorbar(ax, im)
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            ax.set_ylabel("y [$" + xlabel_units_str + "$]")
            if g == 0:
                ax.set_title("Flux Y [$" + flux_units_str + "$]")

            # get flux Z projection
            new_attribute_str = "mass_weighted_radiation_flux" + str(g + 1) + "Z"
            mass_weighted_flux_map = swiftsimio.visualisation.projection.project_gas(
                data,
                project=new_attribute_str,
                **projection_kwargs,
                rotation_matrix=rotation_matrix,
                rotation_center=rotation_center,
            )
            # mass_weighted_flux_map /= mass_map may cause RecursionErrors...
            photon_map = mass_weighted_flux_map / mass_map
            photon_map.convert_to_units(flux_units)

            ax = fig.add_subplot(ngroups, 4, g * 4 + 4)
            im = ax.imshow(photon_map.T, **imshow_kwargs)
            set_colorbar(ax, im)
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            ax.set_ylabel("z [$" + xlabel_units_str + "$]")
            if g == 0:
                ax.set_title("Flux Z [$" + flux_units_str + "$]")

    else:  # plot just energies

        fig = plt.figure(figsize=(5 * ngroups, 5), dpi=200)
        figname = filename[:-5] + ".png"

        for g in range(ngroups):

            # get projection
            new_attribute_str = "mass_weighted_radiation_energy" + str(g + 1)
            mass_weighted_photon_map = swiftsimio.visualisation.projection.project_gas(
                data, project=new_attribute_str, **projection_kwargs
            )
            # mass_weighted_photon_map /= mass_map may cause RecursionErrors...
            photon_map = mass_weighted_photon_map / mass_map
            photon_map.convert_to_units(energy_units)

            ax = fig.add_subplot(1, ngroups, g + 1)

            if energy_boundaries is not None:
                imshow_kwargs["vmin"] = energy_boundaries[g][0]
                imshow_kwargs["vmax"] = energy_boundaries[g][1]
            im = ax.imshow(photon_map.T, **imshow_kwargs)
            set_colorbar(ax, im)
            ax.set_title("Group {0:2d}".format(g + 1))
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            if g == 0:
                ax.set_ylabel("Energies [$" + energy_units_str + "$]")

    # Add title
    title = filename.replace("_", "\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2e}".format(meta.time.to(time_units))
    fig.suptitle(title)

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    gc.collect()

    return


def get_minmax_vals(snaplist):
    """
    Find minimal and maximal values for energy and flux
    so you can fix axes limits over all snapshots

    snaplist: list of snapshot filenames

    returns:

    energy_boundaries: list of [E_min, E_max] for each photon group
    flux_boundaries: list of [F_min, F_max] for each photon group
    """

    emins = []
    emaxs = []
    fmins = []
    fmaxs = []

    ngroups = 0

    for filename in snaplist:

        data = swiftsimio.load(filename)
        meta = data.metadata
        scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

        if do_stromgren_sphere:
            time_units, energy_units, energy_units_str, flux_units, flux_units_str = get_units(
                scheme, unit_system="stromgren_units"
            )
        else:
            time_units, energy_units, energy_units_str, flux_units, flux_units_str = get_units(
                scheme, unit_system="cgs_units"
            )

        ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"])
        emin_group = []
        emax_group = []
        fluxmin_group = []
        fluxmax_group = []

        for g in range(ngroups):

            en = getattr(data.gas.photon_energies, "group" + str(g + 1))
            emax_group.append(en.max())
            mask = en > 0
            if mask.any():
                emin_group.append(en[mask].min())
            else:
                emin_group.append(en.max())

            dirmin = []
            dirmax = []
            for direction in ["X", "Y", "Z"]:
                f = getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + direction)
                dirmin.append(f.min())
                dirmax.append(f.max())
            fluxmin_group.append(min(dirmin))
            fluxmax_group.append(max(dirmax))

        emins.append(emin_group)
        emaxs.append(emax_group)
        fmins.append(fluxmin_group)
        fmaxs.append(fluxmax_group)

    energy_boundaries = []
    flux_boundaries = []
    for g in range(ngroups):
        # Note: snapshot 0 can have emax=0
        emin = min([emins[f][g] for f in range(len(snaplist)) if emins[f][g] > 0])
        emax = max([emaxs[f][g] for f in range(len(snaplist))])
        emin.convert_to_units(energy_units)
        emax.convert_to_units(energy_units)
        energy_boundaries.append([emin, emax])
        fmin = min([fmins[f][g] for f in range(len(snaplist))])
        fmin.convert_to_units(flux_units)
        fmax = max([fmaxs[f][g] for f in range(len(snaplist))])
        fmax.convert_to_units(flux_units)
        flux_boundaries.append([fmin, fmax])

    return energy_boundaries, flux_boundaries


if __name__ == "__main__":

    snaplist = get_snapshot_list(snapshot_base)
    if fancy:
        energy_boundaries, flux_boundaries = get_minmax_vals(snaplist)
    else:
        energy_boundaries = None
        flux_boundaries = None

    for f in snaplist:
        plot_photons(
            f, energy_boundaries=energy_boundaries, flux_boundaries=flux_boundaries
        )
