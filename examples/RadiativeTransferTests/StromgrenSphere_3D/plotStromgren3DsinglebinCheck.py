#!/usr/bin/env python3

# ----------------------------------------------------
# Stromgren 3D with grey approximation (single-frequency bin) and fixed temperature
# The test is identical to Test 1 in Iliev et al. 2006 doi:10.1111/j.1365-2966.2006.10775.x
# Analytic solution is described in Appendix C of SPHM1RT paper (https://arxiv.org/abs/2102.08404)
# Plot comparison of simulated neutral fraction with analytic solution
# ----------------------------------------------------

import sys

import matplotlib as mpl
import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt
from scipy.integrate import odeint

import stromgren_plotting_tools as spt

# Plot parameters
params = {
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "font.size": 14,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.width": 1.5,
    "ytick.major.width": 1.5,
    "axes.linewidth": 1.5,
    "text.usetex": True,
    "figure.figsize": (5, 4),
    "figure.subplot.left": 0.045,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.05,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.12,
    "lines.markersize": 1,
    "lines.linewidth": 2.0,
}
mpl.rcParams.update(params)

scatterplot_kwargs = {
    "alpha": 0.1,
    "s": 2,
    "marker": ".",
    "linewidth": 0.0,
    "facecolor": "blue",
}

# Read in cmdline arg: Are we plotting only one snapshot, or all?
# WARNING: The reference solution is comparable with snapshot_500 only
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True
    snapnr = -1
    print(
        "WARNING: plotting all snapshots, but you should compare the reference solution with the final snapshot only"
    )

snapshot_base = "output_singlebin"


def fn_neutralfraction3d(xn, rn):
    """
    ODE for Stromgren sphere 3D:
    this is the rhs of the ODE to integrate, i.e. dx/drn=fn(x,r)=x*(1-x)/(1+x)*(2/rn+x)
    """
    return xn * (1.0 - xn) / (1.0 + xn) * (2.0 / rn + xn)


# analytic solution
def Stromgren3D_neutralfraction3d(rfunc, nH, sigma, alphaB, dNinj, rini):
    """
    This function integrates the ODE for Stromgren sphere 3D
    Output the neutral fraction xn as a function of radius rfunc
    """
    xn0 = nH * alphaB * 4.0 * np.pi / sigma / dNinj * rini * rini
    rnounit = rfunc * nH * sigma
    xn = odeint(fn_neutralfraction3d, xn0, rnounit)
    return xn


def get_analytic_neutralfraction_stromgren3D(data, scheme):
    """
    This function reads the parameters from snapshot,
    and then integrates the ODE for Stromgren sphere 3D
    Output the neutral fraction xn and radius r_ana
    Work only for SPHM1RT (for now)
    """
    meta = data.metadata
    rho = data.gas.densities
    rini_value = 0.1
    r_ana = np.linspace(rini_value, 10.0, 100) * unyt.kpc
    rini = rini_value * unyt.kpc
    nH = np.mean(rho.to("g/cm**3") / unyt.proton_mass)

    if scheme.startswith("SPH M1closure"):
        sigma_cross = spt.trim_paramstr(
            meta.parameters["SPHM1RT:sigma_cross"].decode("utf-8")
        ) * unyt.unyt_array(1.0, "cm**2")
        sigma = sigma_cross[0]
        alphaB = spt.trim_paramstr(
            meta.parameters["SPHM1RT:alphaB"].decode("utf-8")
        ) * unyt.unyt_array(1.0, "cm**3/s")
    else:
        raise ValueError(
            "Error: Currently get_analytic_neutralfraction_stromgren3D can only work with SPHM1RT"
        )
    units = data.units
    unit_l_in_cgs = units.length.in_cgs()
    unit_v_in_cgs = (units.length / units.time).in_cgs()
    unit_m_in_cgs = units.mass.in_cgs()
    star_emission_rates = (
        spt.trim_paramstr(
            meta.parameters["SPHM1RT:star_emission_rates"].decode("utf-8")
        )
        * unit_m_in_cgs
        * unit_v_in_cgs ** 3
        / unit_l_in_cgs
    )
    ionizing_photon_energy_erg = (
        spt.trim_paramstr(
            meta.parameters["SPHM1RT:ionizing_photon_energy_erg"].decode("utf-8")
        )
        * unyt.erg
    )
    dNinj = star_emission_rates[1] / ionizing_photon_energy_erg[0]
    xn = Stromgren3D_neutralfraction3d(r_ana, nH, sigma, alphaB, dNinj, rini)
    return r_ana, xn


def plot_analytic_compare(filename):
    # Read in data first
    print("working on", filename)
    data = swiftsimio.load(filename)
    meta = data.metadata
    boxsize = meta.boxsize
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    xstar = data.stars.coordinates
    xpart = data.gas.coordinates
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))

    imf = spt.get_imf(scheme, data)
    xHI = imf.HI / (imf.HI + imf.HII)

    r_ana, xn = get_analytic_neutralfraction_stromgren3D(data, scheme)
    plt.scatter(r, xHI, **scatterplot_kwargs)
    plt.plot(r_ana, xn)
    plt.ylabel("Neutral Fraction")
    xlabel_units_str = meta.boxsize.units.latex_representation()
    plt.xlabel("r [$" + xlabel_units_str + "$]")
    plt.yscale("log")
    plt.xlim([0, boxsize[0] / 2.0])
    plt.tight_layout()
    figname = filename[:-5]
    figname += "-Stromgren3Dsinglebin.png"
    plt.savefig(figname, dpi=200)
    plt.close()


if __name__ == "__main__":
    snaplist = spt.get_snapshot_list(snapshot_base, plot_all, snapnr)
    for f in snaplist:
        plot_analytic_compare(f)
