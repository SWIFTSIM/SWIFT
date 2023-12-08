#!/usr/bin/env python3

"""
plotSolution.py

Plot the density, pressure, divB, velocity components and magnetic field
components for the Ryu & Jones 1A snapshot input and create a figure with the
given output name.

Usage:
  python3 plotSolution.py SNAPSHOT OUTPUT

Also plots the "exact" solution that is hard-coded in exact_solution.py.
"""

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import argparse
from exact_solution import exact_RyuJones_1A

# parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
args = argparser.parse_args()

# read the variables of interest from the snapshot file
gamma = None
boxsize = None
t = None
x = None
rho = None
v = None
P = None
u = None
B = None
divB = None
with h5py.File(args.input, "r") as handle:
    gamma = handle["HydroScheme"].attrs["Adiabatic index"][0]
    boxsize = handle["Header"].attrs["BoxSize"][0]
    t = handle["Header"].attrs["Time"][0]
    pos = handle["PartType0/Coordinates"][:]
    rho = handle["PartType0/Densities"][:]
    h = handle["PartType0/SmoothingLengths"][:]
    v = handle["PartType0/Velocities"][:]
    P = handle["PartType0/Pressures"][:]
    u = handle["PartType0/InternalEnergies"][:]
    B = handle["PartType0/MagneticFluxDensities"][:]
    divB = handle["PartType0/MagneticDivergences"][:]
    psi = handle["PartType0/DednerScalars"][:]
    git = handle["Code"].attrs["Git Revision"]
    gitBranch = handle["Code"].attrs["Git Branch"]
    scheme = handle["/HydroScheme"].attrs["Scheme"]
    kernel = handle["/HydroScheme"].attrs["Kernel function"]
    neighbours = handle["/HydroScheme"].attrs["Kernel target N_ngb"]
    mu0 = handle["/PhysicalConstants/InternalUnits"].attrs["vacuum_permeability"]
    dedhyp = handle["/HydroScheme"].attrs["Dedner Hyperbolic Constant"]
    dedpar = handle["/HydroScheme"].attrs["Dedner Parabolic Constant"]
    mhdeta = handle["/HydroScheme"].attrs["Resistive Eta"]
    mhdflavour = handle["/HydroScheme"].attrs["MHD Flavour"]

x = pos[:, 0]
bb = np.sqrt(B[:, 0] * B[:, 0] + B[:, 1] * B[:, 1] + B[:, 2] * B[:, 2])
x -= boxsize * 0.5

# get the exact solution
exx, exsol = exact_RyuJones_1A(t)

# plot everything
fig, ax = pl.subplots(2, 4, sharex=True, figsize=(12, 7))

swift_pts_args = dict(ls="", marker=".", mec="k", ms=0.3)
ex_sol_args = dict(c="r", lw=0.8)

ax[0][0].plot(x, rho, **swift_pts_args)
ax[0][0].plot(exx, exsol["rho"], **ex_sol_args)

ax[0][1].plot(x, v[:, 0], **swift_pts_args)
ax[0][1].plot(exx, exsol["u"], **ex_sol_args)

ax[0][2].plot(x, v[:, 1], **swift_pts_args)
ax[0][2].plot(exx, exsol["v"], **ex_sol_args)

ax[0][3].plot(x, np.log10(abs(divB) * h / bb), **swift_pts_args)
ax03 = ax[0][3].twinx()
ax03.plot(x, psi, "b.", ms=0.3)

ax[1][0].plot(x, B[:, 0], label="SWIFT", **swift_pts_args)
ax[1][0].plot(exx, exsol["Bx"], label="Exact solution", **ex_sol_args)

ax[1][1].plot(x, B[:, 1], **swift_pts_args)
ax[1][1].plot(exx, exsol["By"], **ex_sol_args)

ax[1][2].plot(x, P, **swift_pts_args)
ax[1][2].plot(exx, exsol["P"], **ex_sol_args)

ax[0][0].set_ylabel(r"$\rho$")
ax[0][1].set_ylabel(r"$v_x$")
ax[0][2].set_ylabel(r"$v_y$")
ax[0][3].set_ylabel(r"$err_{\nabla \cdot \mathbf{B}}$")
ax03.set_ylabel(r"$\psi$")
ax[1][0].set_ylabel(r"$B_x$")
ax[1][1].set_ylabel(r"$B_y$")
ax[1][2].set_ylabel(r"$u$")

ax[1][0].legend(loc="best")

# Information -------------------------------------

text_fontsize = 7

ax[1][3].text(
    -0.45,
    0.8,
    "Ryu & Jones 1A shock Tube with $\\gamma=%.3f$ in 3D\nat $t=%.2f$" % (gamma, t),
    fontsize=text_fontsize,
)
ax[1][3].plot([-0.45, 0.45], [0.62, 0.62], "k-", lw=1)
ax[1][3].text(-0.45, 0.5, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=text_fontsize)
"""
ax[1][3].text(
    -0.45, 0.4, "$Branch$ %s" % gitBranch.decode("utf-8"), fontsize=text_fontsize
)
"""
ax[1][3].text(-0.45, 0.3, scheme.decode("utf-8"), fontsize=text_fontsize)
ax[1][3].text(-0.45, 0.2, kernel.decode("utf-8"), fontsize=text_fontsize)
ax[1][3].text(-0.45, 0.1, "$%.2f$ neighbours" % (neighbours), fontsize=text_fontsize)
ax[1][3].plot([-0.45, 0.45], [0.0, 0.0], "k-", lw=1)
ax[1][3].text(-0.45, -0.1, "$\\mu_0:%.2f$ " % (mu0), fontsize=text_fontsize)
ax[1][3].text(-0.45, -0.2, "$Difussion_\\eta:%.2f$ " % (mhdeta), fontsize=text_fontsize)
ax[1][3].text(
    -0.45,
    -0.3,
    "$Dedner: Hyp=%.2f // Par=%.2f$" % (dedhyp, dedpar),
    fontsize=text_fontsize,
)
ax[1][3].plot([-0.45, 0.45], [-0.4, -0.4], "k-", lw=1)
ax[1][3].text(
    -0.45,
    -0.5,
    "$Flavour: $ %s" % mhdflavour.decode("utf-8")[0:30],
    fontsize=text_fontsize,
)

for a in ax[1:-1]:
    a.set_xlabel("x")

"""
# make sure all velocity and magnetic fields plots use the same y axis
# this puts the noise on constant values into context
for a in [ax[1][0], ax[1][2]]:
   a.set_ylim(*ax[1][1].get_ylim())
for a in [ax[2][0], ax[2][2]]:
   a.set_ylim(*ax[2][1].get_ylim())
ax[0][0].set_ylim(0.0, 1.1)
ax[0][1].set_ylim(0.0, 1.1)
ax[0][2].set_ylim(-0.2, 0.2)
ax[1][0].set_ylim(-2.0, 1.0)
ax[1][1].set_ylim(-2.0, 1.0)
ax[1][2].set_ylim(-2.0, 1.0)
ax[2][0].set_ylim(-1.1, 1.1)
ax[2][1].set_ylim(-1.1, 1.1)
ax[2][2].set_ylim(-1.1, 1.1)
"""

ax[1][3].tick_params(
    left=False, right=False, labelleft=False, labelbottom=False, bottom=False
)
ax[1][3].set_xlabel("")
ax[1][3].plot(frameon=False)

# mark the validity area: [0,1] and [3,4] contain the solution of the
# mirrored Riemann problem across the periodic boundary
for a in ax.flatten():
    a.axvline(x=-0.5, linestyle="--", color="k")
    a.axvline(x=0.5, linestyle="--", color="k")
# only plot the relevant part of the solution

ax[0][0].set_xlim(-0.6, 0.6)

# add the time as a title
# ax[0][1].set_title(f"t={t:.2e}")

# save the figure
pl.tight_layout()
pl.savefig(args.output, dpi=100)
