#!/usr/bin/env python3

import os
import sys
import argparse

import matplotlib as mpl
import numpy as np
import swiftsimio
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as cf

# Parameters users should/may tweak
snapshot_base = "output"  # snapshot basename
plot_physical_quantities = True

mpl.rcParams["text.usetex"] = True


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-z", "--redshift", help="Redshift domain to plot advection for", default="high"
    )

    args = parser.parse_args()
    return args


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshots that are to be plotted 
    and return their names as list
    """
    snaplist = []

    dirlist = os.listdir()
    for f in dirlist:
        if f.startswith(snapshot_basename) and f.endswith("hdf5"):
            snaplist.append(f)

    snaplist = sorted(snaplist)
    if len(snaplist) == 0:
        print(f"No snapshots with base {snapshot_basename} found!")
        sys.exit(1)
    return snaplist


def plot_param_over_time(
    snapshot_list, param="energy density", redshift_domain="high_redshift"
):
    print(f"Now plotting {param} over time")

    # Arrays to keep track of plot_param and scale factor
    plot_param = [[], []]
    scale_factor = []
    analytic_exponent = [0, 0]

    # Functions to convert between scale factor and redshift
    a2z = lambda a: 1 / a - 1
    z2a = lambda z: 1 / (z + 1)

    for file in snapshot_list:
        data = swiftsimio.load(file)
        meta = data.metadata

        # Read comoving quantities
        energy = getattr(data.gas.photon_energies, "group1")
        mass = data.gas.masses
        rho = data.gas.densities
        vol = mass / rho

        energy_density = energy / vol

        if plot_physical_quantities:
            physical_energy_density = energy_density.to_physical()
            physical_mass = mass.to_physical()
            physical_vol = vol.to_physical()
            physical_energy = physical_energy_density * physical_vol
            if param == "energy density":
                plot_param[1].append(
                    1
                    * np.sum(physical_energy_density)
                    / physical_energy_density.shape[0]
                )
                analytic_exponent[1] = -3.0
            elif param == "volume":
                plot_param[1].append(1 * np.sum(physical_vol) / physical_vol.shape[0])
                analytic_exponent[1] = 3.0
            elif param == "total energy":
                plot_param[1].append(1 * np.sum(physical_energy))
                analytic_exponent[1] = 0.0
            elif param == "mass":
                plot_param[1].append(1 * np.sum(physical_mass))
                analytic_exponent[1] = 0.0

        if param == "energy density":
            plot_param[0].append(1 * np.sum(energy_density) / energy_density.shape[0])
            analytic_exponent[0] = 0.0
        elif param == "volume":
            plot_param[0].append(1 * np.sum(vol) / vol.shape[0])
            analytic_exponent[0] = 0.0
        elif param == "total energy":
            plot_param[0].append(1 * np.sum(energy))
            analytic_exponent[0] = 0.0
        elif param == "mass":
            plot_param[0].append(1 * np.sum(mass))
            analytic_exponent[0] = 0.0

        scale_factor.append(meta.scale_factor)

    fig = plt.figure(figsize=(5.05 * (1 + plot_physical_quantities), 5.4), dpi=200)

    x = np.linspace(min(scale_factor), max(scale_factor), 1000)

    if param == "energy density":
        titles = ["Comoving energy density", "Physical energy density"]
        ylabel = "Average energy density"
        figname = f"output_energy_density_over_time-{redshift_domain}.png"
    elif param == "volume":
        titles = ["Comoving particle volume", "Physical particle volume"]
        ylabel = "Average particle volume"
        figname = f"output_volume_over_time-{redshift_domain}.png"
    elif param == "total energy":
        titles = ["Comoving total energy", "Physical total energy"]
        ylabel = "Total energy"
        figname = f"output_total_energy_over_time-{redshift_domain}.png"
    elif param == "mass":
        titles = ["Comoving total mass", "Physical total mass"]
        ylabel = "Total mass"
        figname = f"output_total_mass_over_time-{redshift_domain}.png"

    for i in range(1 + plot_physical_quantities):
        ax = fig.add_subplot(1, (1 + plot_physical_quantities), (1 + i))
        ax.scatter(scale_factor, plot_param[i], label="Simulation")

        # Analytic scale-factor relation
        analytic = x ** analytic_exponent[i]

        # Scale solution to correct offset
        analytic = analytic / analytic[0] * plot_param[i][0]
        ax.plot(
            x,
            analytic,
            c="r",
            label=f"Analytic solution $\propto a^{{{analytic_exponent[i]}}}$",
        )

        ax.legend()
        ax.set_title(titles[i])

        ax.set_xlabel("Scale factor")
        secax = ax.secondary_xaxis("top", functions=(a2z, z2a))
        secax.set_xlabel("Redshift")

        ax.yaxis.get_offset_text().set_position((-0.05, 1))

        if analytic_exponent[i] == 0.0:
            ax.set_ylim(plot_param[i][0] * 0.95, plot_param[i][0] * 1.05)
        if i == 0:
            units = plot_param[i][0].units.latex_representation()
            ax.set_ylabel(f"{ylabel} [${units}$]")

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()


if __name__ in ("__main__"):
    # Get command line args
    args = parse_args()
    domain = args.redshift.lower()
    if domain in ("low", "l", "low_redshift", "low redshift", "low-redshift"):
        redshift_domain = "low_redshift"
    elif domain in (
        "medium",
        "m",
        "medium_redshift",
        "medium redshift",
        "medium-redshift",
    ):
        redshift_domain = "medium_redshift"
    elif domain in ("high", "h", "high_redshift", "high redshift", "high-redshift"):
        redshift_domain = "high_redshift"
    else:
        print("Redshift domain not recognised!")
        sys.exit(1)

    snaplist = get_snapshot_list(snapshot_base + f"_{redshift_domain}")

    if len(snaplist) < 1:
        print("No snapshots found!")
        exit(1)

    for param in ["energy density", "volume", "total energy", "mass"]:
        plot_param_over_time(snaplist, param, redshift_domain)
