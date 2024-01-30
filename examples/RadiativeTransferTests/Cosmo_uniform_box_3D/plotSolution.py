#!/usr/bin/env python3

import os
import sys

import matplotlib as mpl
import numpy as np
import swiftsimio
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as cf

# Parameters users should/may tweak
snapshot_base = "output" # snapshot basename
plot_physical_energies = True

mpl.rcParams["text.usetex"] = True

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
    return snaplist

def plot_param_over_time(snapshot_list, param="energy density"):
    print(f"Now plotting {param} over time")

    # Arrays to keep track of plot_param and scale factor
    plot_param = [[],[]]
    scale_factor = []


    # Functions to convert between scale factor and redshift
    a2z = lambda a: 1/a - 1
    z2a = lambda z: 1/(z+1)

    for file in snapshot_list:
        data = swiftsimio.load(file)
        meta = data.metadata
        
        energy_density = data.gas.photon_energy_densities
        mass = data.gas.masses
        rho = data.gas.densities
        vol = mass / rho

        energy = energy_density * vol

        if plot_physical_energies:
            physical_energy_density = energy_density.to_physical()
            physical_vol = vol.to_physical()
            physical_energy = physical_energy_density * physical_vol
            match param:
                case "energy density":
                    plot_param[1].append(np.sum(physical_energy_density) / len(physical_energy_density))
                case "volume":
                    plot_param[1].append(np.sum(physical_vol) / len(physical_vol))
                case "total energy":
                    plot_param[1].append(np.sum(physical_energy))

        match param:
            case "energy density":
                plot_param[0].append(np.sum(energy_density) / len(energy_density))
            case "volume":
                plot_param[0].append(np.sum(vol) / len(vol))
            case "total energy":
                plot_param[0].append(np.sum(energy))
        scale_factor.append(meta.scale_factor)

    fig = plt.figure(figsize=(5.05 * (1 + plot_physical_energies), 5.4), dpi=200)

    def fit_func(x,a,b):
        return a*x+b

    x = np.linspace(min(scale_factor), max(scale_factor), 1000)
    
    match param:
        case "energy density":
            titles = ["Comoving energy density", "Physical energy density"]
            ylabel = "Average energy density"
            figname = "output_energy_density_over_time.png"
        case "volume":
            titles = ["Comoving particle volume", "Physical particle volume"]
            ylabel = "Average particle volume"
            figname = "output_volume_over_time.png"
        case "total energy":
            titles = ["Comoving total energy", "Physical total energy"]
            ylabel = "Total energy"
            figname = "output_total_energy_over_time.png"
    for i in range(1+plot_physical_energies):
        # Get exponent of scale factor
        popt, __ = cf(fit_func, np.log10(scale_factor), np.log10(plot_param[i]))

        ax = fig.add_subplot(1,(1 + plot_physical_energies), (1 + i))
        ax.scatter(scale_factor, plot_param[i], label="Simulation")
        
        if not np.isclose(popt[0], 0, atol=1e-4):
            ax.plot(x, 10**fit_func(np.log10(x), *popt), c="r", label=f"$\propto a^{{{popt[0]:.2f}}}$")
        ax.legend()
        ax.set_title(titles[i])
        # ax.set_ylim(5.2e14, 5.3e14)

        ax.set_xlabel("Scale factor")
        secax = ax.secondary_xaxis("top", functions=(a2z, z2a))
        secax.set_xlabel("Redshift")

        ax.yaxis.get_offset_text().set_position((-0.05, 1))
        
        if i == 0:
            units = plot_param[i][0].units.latex_representation()
            ax.set_ylabel(f"{ylabel} [${units}$]")

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()


if __name__ in ("__main__"):
    snaplist = get_snapshot_list(snapshot_base)

    if len(snaplist) < 1:
        print("No snapshots found!")
        exit(1)
    
    for param in ["energy density", "volume", "total energy"]:
        plot_param_over_time(snaplist, param)
