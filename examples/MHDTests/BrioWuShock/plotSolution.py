#!/usr/bin/env python3

"""
plotSolution.py

Plot the density, pressure, divB, velocity components and magnetic field
components for the Brio Wu snapshot input and create a figure with the
given output name.

Usage:
  python3 plotSolution.py SNAPSHOT OUTPUT [--ncell NCELL]
where NCELL is the number of cells to use for the HLL Riemann solver reference
solution (default: 1000).

Also plots the "exact" solution that is hard-coded in exact_solution.py.
"""

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import argparse
from hll_solver import solve_MHD_Riemann_problem
from exact_solution import exact_Brio_Wu

# parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
argparser.add_argument("--ncell", "-n", type=int, default=1000)
args = argparser.parse_args()

# read the variables of interest from the snapshot file
gamma = None
boxsize = None
t = None
x = None
rho = None
v = None
P = None
B = None
divB = None
with h5py.File(args.input, "r") as handle:
    gamma = handle["HydroScheme"].attrs["Adiabatic index"][0]
    boxsize = handle["Header"].attrs["BoxSize"][0]
    t = handle["Header"].attrs["Time"][0]
    x = handle["PartType0/Coordinates"][:, 0]
    rho = handle["PartType0/Densities"][:]
    v = handle["PartType0/Velocities"][:]
    P = handle["PartType0/Pressures"][:]
    B = handle["PartType0/MagneticFluxDensity"][:]
    divB = handle["PartType0/MagneticDivergence"][:]

# recreate the same Riemann problem for the HLL solver
# the initial condition is generated using gamma=2., but since we set the
# internal energy and not the pressure in the IC, this might lead to a different
# initial pressure after SWIFT reads in the internal energy if SWIFT was
# compiled with a different gamma
gamma_ic = 2.0
riemann_problem = {
    "left_state": {
        "rho": 1.0,
        "p": 1.0 * (gamma - 1.0) / (gamma_ic - 1.0),
        "u": 0.0,
        "v": 0.0,
        "w": 0.0,
        "Bx": 0.75,
        "By": 1.0,
        "Bz": 0.0,
    },
    "right_state": {
        "rho": 0.125,
        "p": 0.1 * (gamma - 1.0) / (gamma_ic - 1.0),
        "u": 0.0,
        "v": 0.0,
        "w": 0.0,
        "Bx": 0.75,
        "By": -1.0,
        "Bz": 0.0,
    },
    "time": t,
    "dt": t / args.ncell,
    "xsize": 0.5 * boxsize,
    "gamma": 5.0 / 3.0,
    "ncell": args.ncell,
    "solver": "HLL",
}

# compute the HLL reference solution
print("Computing HLL reference solution...")
refx, solution = solve_MHD_Riemann_problem(riemann_problem)
print("Done.")
refx += 0.5 * boxsize

# get the exact solution
exx, exsol = exact_Brio_Wu(t)
exx += 2.0

# plot everything
fig, ax = pl.subplots(3, 3, sharex=True, figsize=(10, 9))

ax[0][0].plot(x, rho, ".")
ax[0][0].plot(refx, solution["rho"], "-")
ax[0][0].plot(exx, exsol["rho"], "-")
ax[0][1].plot(x, P, ".")
ax[0][1].plot(refx, solution["p"], "-")
ax[0][1].plot(exx, exsol["p"], "-")
ax[0][2].plot(x, divB, ".")
ax[0][2].plot(refx[[0, -1]], [0.0, 0.0], "-")
ax[0][2].plot(exx[[0, -1]], [0.0, 0.0], "-")

for i in range(3):
    ax[1][i].plot(x, v[:, i], ".")
    ax[1][i].plot(refx, solution[["u", "v", "w"][i]], "-")
    ax[1][i].plot(exx, exsol[["u", "v", "w"][i]], "-")
    ax[2][i].plot(x, B[:, i], ".", label="SWIFT")
    ax[2][i].plot(
        refx,
        solution[["Bx", "By", "Bz"][i]],
        "-",
        label=f"HLL solver ({args.ncell} cells)",
    )
    ax[2][i].plot(exx, exsol[["Bx", "By", "Bz"][i]], "-", label="Exact solution")

ax[0][0].set_ylabel("rho")
ax[0][1].set_ylabel("P")
ax[0][2].set_ylabel("divB")
ax[1][0].set_ylabel("vx")
ax[1][1].set_ylabel("vy")
ax[1][2].set_ylabel("vz")
ax[2][0].set_ylabel("Bx")
ax[2][1].set_ylabel("By")
ax[2][2].set_ylabel("Bz")

ax[2][0].legend(loc="best")

for a in ax[2]:
    a.set_xlabel("x")

# make sure all velocity and magnetic fields plots use the same y axis
# this puts the noise on constant values into context
for a in [ax[1][0], ax[1][2]]:
    a.set_ylim(*ax[1][1].get_ylim())
for a in [ax[2][0], ax[2][2]]:
    a.set_ylim(*ax[2][1].get_ylim())

# mark the validity area: [0,1] and [3,4] contain the solution of the
# mirrored Riemann problem across the periodic boundary
for a in ax.flatten():
    a.axvline(x=1.0, linestyle="--", color="k")
    a.axvline(x=3.0, linestyle="--", color="k")
# only plot the relevant part of the solution
ax[0][0].set_xlim(0.9, 3.1)

# add the time as a title
ax[0][1].set_title(f"t={t:.2e}")

# save the figure
pl.tight_layout()
pl.savefig(args.output, dpi=300)
