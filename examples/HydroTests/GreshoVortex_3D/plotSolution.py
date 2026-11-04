################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Computes the analytical solution of the Gresho-Chan vortex and plots the SPH
# answer

# Parameters
gas_gamma = 5.0 / 3.0  # Gas adiabatic index
rho0 = 1  # Gas density
P0 = 0.0  # Constant additional pressure (should have no impact on the
# dynamics)

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

import matplotlib

matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import h5py

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")

snap = int(plt.sys.argv[1])

# Generate the analytic solution at this time
N = 200
R_max = 0.8
solution_r = np.arange(0, R_max, R_max / N)
solution_P = np.zeros(N)
solution_v_phi = np.zeros(N)
solution_v_r = np.zeros(N)

for i in range(N):
    if solution_r[i] < 0.2:
        solution_P[i] = P0 + 5.0 + 12.5 * solution_r[i] ** 2
        solution_v_phi[i] = 5.0 * solution_r[i]
    elif solution_r[i] < 0.4:
        solution_P[i] = (
            P0
            + 9.0
            + 12.5 * solution_r[i] ** 2
            - 20.0 * solution_r[i]
            + 4.0 * np.log(solution_r[i] / 0.2)
        )
        solution_v_phi[i] = 2.0 - 5.0 * solution_r[i]
    else:
        solution_P[i] = P0 + 3.0 + 4.0 * np.log(2.0)
        solution_v_phi[i] = 0.0

solution_rho = np.ones(N) * rho0
solution_s = solution_P / solution_rho**gas_gamma
solution_u = solution_P / ((gas_gamma - 1.0) * solution_rho)

# Read the simulation data
sim = h5py.File("gresho_%04d.hdf5" % snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"][0]
eta = sim["/HydroScheme"].attrs["Kernel eta"][0]
git = sim["Code"].attrs["Git Revision"]

pos = sim["/PartType0/Coordinates"][:, :]
x = pos[:, 0] - boxSize / 2
y = pos[:, 1] - boxSize / 2
vel = sim["/PartType0/Velocities"][:, :]
r = np.sqrt(x**2 + y**2)
v_r = (x * vel[:, 0] + y * vel[:, 1]) / r
v_phi = (-y * vel[:, 0] + x * vel[:, 1]) / r
v_norm = np.sqrt(vel[:, 0] ** 2 + vel[:, 1] ** 2)
rho = sim["/PartType0/Densities"][:]
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]

try:
    diffusion = sim["/PartType0/DiffusionParameters"][:]
    plot_diffusion = True
except:
    plot_diffusion = False

try:
    viscosity = sim["/PartType0/ViscosityParameters"][:]
    plot_viscosity = True
except:
    plot_viscosity = False

# Bin the data
r_bin_edge = np.arange(0.0, 1.0, 0.02)
r_bin = 0.5 * (r_bin_edge[1:] + r_bin_edge[:-1])
rho_bin, _, _ = stats.binned_statistic(r, rho, statistic="mean", bins=r_bin_edge)
v_bin, _, _ = stats.binned_statistic(r, v_phi, statistic="mean", bins=r_bin_edge)
P_bin, _, _ = stats.binned_statistic(r, P, statistic="mean", bins=r_bin_edge)
S_bin, _, _ = stats.binned_statistic(r, S, statistic="mean", bins=r_bin_edge)
u_bin, _, _ = stats.binned_statistic(r, u, statistic="mean", bins=r_bin_edge)
rho2_bin, _, _ = stats.binned_statistic(r, rho**2, statistic="mean", bins=r_bin_edge)
v2_bin, _, _ = stats.binned_statistic(r, v_phi**2, statistic="mean", bins=r_bin_edge)
P2_bin, _, _ = stats.binned_statistic(r, P**2, statistic="mean", bins=r_bin_edge)
S2_bin, _, _ = stats.binned_statistic(r, S**2, statistic="mean", bins=r_bin_edge)
u2_bin, _, _ = stats.binned_statistic(r, u**2, statistic="mean", bins=r_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin**2)
v_sigma_bin = np.sqrt(v2_bin - v_bin**2)
P_sigma_bin = np.sqrt(P2_bin - P_bin**2)
S_sigma_bin = np.sqrt(S2_bin - S_bin**2)
u_sigma_bin = np.sqrt(u2_bin - u_bin**2)

if plot_diffusion:
    alpha_diff_bin, _, _ = stats.binned_statistic(
        r, diffusion, statistic="mean", bins=r_bin_edge
    )
    alpha2_diff_bin, _, _ = stats.binned_statistic(
        r, diffusion**2, statistic="mean", bins=r_bin_edge
    )
    alpha_diff_sigma_bin = np.sqrt(alpha2_diff_bin - alpha_diff_bin**2)

if plot_viscosity:
    alpha_visc_bin, _, _ = stats.binned_statistic(
        r, viscosity, statistic="mean", bins=r_bin_edge
    )
    alpha2_visc_bin, _, _ = stats.binned_statistic(
        r, viscosity**2, statistic="mean", bins=r_bin_edge
    )
    alpha_visc_sigma_bin = np.sqrt(alpha2_visc_bin - alpha_visc_bin**2)

# Plot the interesting quantities
plt.figure(figsize=(7, 7 / 1.6))

line_color = "C4"
binned_color = "C2"
binned_marker_size = 4

scatter_props = dict(
    marker=".",
    ms=1,
    markeredgecolor="none",
    alpha=0.1,
    zorder=-1,
    rasterized=True,
    linestyle="none",
)

errorbar_props = dict(color=binned_color, ms=binned_marker_size, fmt=".", lw=1.2)


# Azimuthal velocity profile -----------------------------
plt.subplot(231)

plt.plot(r, v_phi, **scatter_props)
plt.plot(solution_r, solution_v_phi, "--", color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(r_bin, v_bin, yerr=v_sigma_bin, **errorbar_props)
plt.plot([0.2, 0.2], [-100, 100], ":", color=line_color, alpha=0.4, lw=1.2)
plt.plot([0.4, 0.4], [-100, 100], ":", color=line_color, alpha=0.4, lw=1.2)
plt.xlabel("Radius $r$")
plt.ylabel("Azimuthal velocity $v_\\phi$")
plt.xlim(0, R_max)
plt.ylim(-0.1, 1.2)

# Radial density profile --------------------------------
plt.subplot(232)

plt.plot(r, rho, **scatter_props)
plt.plot(solution_r, solution_rho, "--", color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(r_bin, rho_bin, yerr=rho_sigma_bin, **errorbar_props)
plt.plot([0.2, 0.2], [-100, 100], ":", color=line_color, alpha=0.4, lw=1.2)
plt.plot([0.4, 0.4], [-100, 100], ":", color=line_color, alpha=0.4, lw=1.2)
plt.xlabel("Radius $r$")
plt.ylabel("Density $\\rho$")
plt.xlim(0, R_max)
plt.ylim(rho0 - 0.3, rho0 + 0.3)
# yticks([-0.2, -0.1, 0., 0.1, 0.2])

# Radial pressure profile --------------------------------
plt.subplot(233)

plt.plot(r, P, **scatter_props)
plt.plot(solution_r, solution_P, "--", color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(r_bin, P_bin, yerr=P_sigma_bin, **errorbar_props)
plt.plot([0.2, 0.2], [-100, 100], ":", color=line_color, alpha=0.4, lw=1.2)
plt.plot([0.4, 0.4], [-100, 100], ":", color=line_color, alpha=0.4, lw=1.2)
plt.xlabel("Radius $r$")
plt.ylabel("Pressure $P$")
plt.xlim(0, R_max)
plt.ylim(4.9 + P0, P0 + 6.1)

# Internal energy profile --------------------------------
plt.subplot(234)

plt.plot(r, u, **scatter_props)
plt.plot(solution_r, solution_u, "--", color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(r_bin, u_bin, yerr=u_sigma_bin, **errorbar_props)
plt.plot([0.2, 0.2], [-100, 100], ":", color=line_color, alpha=0.4, lw=1.2)
plt.plot([0.4, 0.4], [-100, 100], ":", color=line_color, alpha=0.4, lw=1.2)
plt.xlabel("Radius $r$")
plt.ylabel("Internal Energy $u$")
plt.xlim(0, R_max)
plt.ylim(7.3, 9.1)


# Radial entropy profile --------------------------------
plt.subplot(235)
plt.xlabel("Radius $r$")
plt.xlim(0, R_max)


if plot_diffusion or plot_viscosity:
    if plot_diffusion:
        plt.plot(r, diffusion, **scatter_props)
        plt.errorbar(
            r_bin,
            alpha_diff_bin,
            yerr=alpha_diff_sigma_bin,
            **errorbar_props,
            label="Diffusion",
        )

    if plot_viscosity:
        plt.plot(r, viscosity, **{**scatter_props, "color": "C3"})
        plt.errorbar(
            r_bin,
            alpha_visc_bin,
            yerr=alpha_visc_sigma_bin,
            **{**errorbar_props, "color": "C4"},
            label="Viscosity",
        )

    plt.ylabel("Rate Coefficient $\\alpha$")
    plt.legend()
else:
    plt.plot(r, S, **scatter_props)
    plt.plot(solution_r, solution_s, "--", color=line_color, alpha=0.8, lw=1.2)
    plt.errorbar(r_bin, S_bin, yerr=S_sigma_bin, **errorbar_props)
    plt.plot([0.2, 0.2], [-100, 100], ":", color=line_color, alpha=0.4, lw=1.2)
    plt.plot([0.4, 0.4], [-100, 100], ":", color=line_color, alpha=0.4, lw=1.2)
    plt.ylabel("Entropy $S$")
    plt.ylim(4.9 + P0, P0 + 6.1)


# Information -------------------------------------
plt.subplot(236, frameon=False)

text_fontsize = 5

plt.text(
    -0.49,
    0.9,
    "Gresho-Chan vortex (3D) with $\\gamma=%.3f$ at $t=%.2f$" % (gas_gamma, time),
    fontsize=text_fontsize,
)
plt.text(-0.49, 0.8, "Background $\\rho_0=%.3f$" % rho0, fontsize=text_fontsize)
plt.text(-0.49, 0.7, "Background $P_0=%.3f$" % P0, fontsize=text_fontsize)
plt.plot([-0.49, 0.1], [0.62, 0.62], "k-", lw=1)
plt.text(-0.49, 0.5, "SWIFT %s" % git.decode("utf-8"), fontsize=text_fontsize)
plt.text(-0.49, 0.4, scheme.decode("utf-8"), fontsize=text_fontsize)
plt.text(-0.49, 0.3, kernel.decode("utf-8"), fontsize=text_fontsize)
plt.text(
    -0.49,
    0.2,
    "$%.2f$ neighbours ($\\eta=%.3f$)" % (neighbours, eta),
    fontsize=text_fontsize,
)
plt.xlim(-0.5, 0.5)
plt.ylim(0, 1)
plt.xticks([])
plt.yticks([])

plt.tight_layout()

plt.savefig("GreshoVortex.png")
