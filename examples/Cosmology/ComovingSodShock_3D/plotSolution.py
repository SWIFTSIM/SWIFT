###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Computes the analytical solution of the Sod shock and plots the SPH answer


# Generates the analytical  solution for the Sod shock test case
# The script works for a given left (x<0) and right (x>0) state and computes the solution at a later time t.
# This follows the solution given in (Toro, 2009)


# Parameters
gas_gamma = 5.0 / 3.0  # Polytropic index
rho_L = 1.0  # Density left state
rho_R = 0.125  # Density right state
v_L = 0.0  # Velocity left state
v_R = 0.0  # Velocity right state
P_L = 1.0  # Pressure left state
P_R = 0.1  # Pressure right state

import sys

sys.path.append("../../HydroTests/")
from riemannSolver import RiemannSolver

import matplotlib

matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py

style.use("../../../tools/stylesheets/mnras.mplstyle")

snap = int(sys.argv[1])


# Read the simulation data
sim = h5py.File("sodShock_%04d.hdf5" % snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
anow = sim["/Header"].attrs["Scale-factor"]
a_i = sim["/Cosmology"].attrs["a_beg"]
H_0 = sim["/Cosmology"].attrs["H0 [internal units]"]
time = 2.0 * (1.0 / np.sqrt(a_i) - 1.0 / np.sqrt(anow)) / H_0
scheme = sim["/HydroScheme"].attrs["Scheme"].decode("utf-8")
kernel = sim["/HydroScheme"].attrs["Kernel function"].decode("utf-8")
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"].decode("utf-8")

x = sim["/PartType0/Coordinates"][:, 0]
v = sim["/PartType0/Velocities"][:, 0] * anow
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]
rho = sim["/PartType0/Densities"][:]

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

x_min = -1.0
x_max = 1.0
x += x_min
N = 1000

# Bin the data
x_bin_edge = np.arange(-0.6, 0.6, 0.02)
x_bin = 0.5 * (x_bin_edge[1:] + x_bin_edge[:-1])
rho_bin, _, _ = stats.binned_statistic(x, rho, statistic="mean", bins=x_bin_edge)
v_bin, _, _ = stats.binned_statistic(x, v, statistic="mean", bins=x_bin_edge)
P_bin, _, _ = stats.binned_statistic(x, P, statistic="mean", bins=x_bin_edge)
S_bin, _, _ = stats.binned_statistic(x, S, statistic="mean", bins=x_bin_edge)
u_bin, _, _ = stats.binned_statistic(x, u, statistic="mean", bins=x_bin_edge)
rho2_bin, _, _ = stats.binned_statistic(x, rho ** 2, statistic="mean", bins=x_bin_edge)
v2_bin, _, _ = stats.binned_statistic(x, v ** 2, statistic="mean", bins=x_bin_edge)
P2_bin, _, _ = stats.binned_statistic(x, P ** 2, statistic="mean", bins=x_bin_edge)
S2_bin, _, _ = stats.binned_statistic(x, S ** 2, statistic="mean", bins=x_bin_edge)
u2_bin, _, _ = stats.binned_statistic(x, u ** 2, statistic="mean", bins=x_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin ** 2)
v_sigma_bin = np.sqrt(v2_bin - v_bin ** 2)
P_sigma_bin = np.sqrt(P2_bin - P_bin ** 2)
S_sigma_bin = np.sqrt(S2_bin - S_bin ** 2)
u_sigma_bin = np.sqrt(u2_bin - u_bin ** 2)

if plot_diffusion:
    alpha_diff_bin, _, _ = stats.binned_statistic(
        x, diffusion, statistic="mean", bins=x_bin_edge
    )
    alpha2_diff_bin, _, _ = stats.binned_statistic(
        x, diffusion ** 2, statistic="mean", bins=x_bin_edge
    )
    alpha_diff_sigma_bin = np.sqrt(alpha2_diff_bin - alpha_diff_bin ** 2)

if plot_viscosity:
    alpha_visc_bin, _, _ = stats.binned_statistic(
        x, viscosity, statistic="mean", bins=x_bin_edge
    )
    alpha2_visc_bin, _, _ = stats.binned_statistic(
        x, viscosity ** 2, statistic="mean", bins=x_bin_edge
    )
    alpha_visc_sigma_bin = np.sqrt(alpha2_visc_bin - alpha_visc_bin ** 2)

# Prepare reference solution
solver = RiemannSolver(gas_gamma)

delta_x = (x_max - x_min) / N
x_s = arange(0.5 * x_min, 0.5 * x_max, delta_x)
rho_s, v_s, P_s, _ = solver.solve(rho_L, v_L, P_L, rho_R, v_R, P_R, x_s / time)

# Additional arrays
u_s = P_s / (rho_s * (gas_gamma - 1.0))  # internal energy
s_s = P_s / rho_s ** gas_gamma  # entropic function

# Plot the interesting quantities
figure(figsize=(7, 7 / 1.6))

# Velocity profile --------------------------------
subplot(231)
plot(x, v, ".", color="r", ms=0.5, alpha=0.2)
plot(x_s, v_s, "--", color="k", alpha=0.8, lw=1.2)
errorbar(x_bin, v_bin, yerr=v_sigma_bin, fmt=".", ms=8.0, color="b", lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Velocity}}~v_x$", labelpad=0)
xlim(-0.5, 0.5)
ylim(-0.1, 0.95)

# Density profile --------------------------------
subplot(232)
plot(x, rho, ".", color="r", ms=0.5, alpha=0.2)
plot(x_s, rho_s, "--", color="k", alpha=0.8, lw=1.2)
errorbar(x_bin, rho_bin, yerr=rho_sigma_bin, fmt=".", ms=8.0, color="b", lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho$", labelpad=0)
xlim(-0.5, 0.5)
ylim(0.05, 1.1)

# Pressure profile --------------------------------
subplot(233)
plot(x, P, ".", color="r", ms=0.5, alpha=0.2)
plot(x_s, P_s, "--", color="k", alpha=0.8, lw=1.2)
errorbar(x_bin, P_bin, yerr=P_sigma_bin, fmt=".", ms=8.0, color="b", lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Pressure}}~P$", labelpad=0)
xlim(-0.5, 0.5)
ylim(0.01, 1.1)

# Internal energy profile -------------------------
subplot(234)
plot(x, u, ".", color="r", ms=0.5, alpha=0.2)
plot(x_s, u_s, "--", color="k", alpha=0.8, lw=1.2)
errorbar(x_bin, u_bin, yerr=u_sigma_bin, fmt=".", ms=8.0, color="b", lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
xlim(-0.5, 0.5)
ylim(0.8, 2.2)

# Entropy profile ---------------------------------
subplot(235)
xlim(-0.5, 0.5)
xlabel("${\\rm{Position}}~x$", labelpad=0)

if plot_diffusion or plot_viscosity:
    if plot_diffusion:
        plot(x, diffusion * 100, ".", color="r", ms=0.5, alpha=0.2)
        errorbar(
            x_bin,
            alpha_diff_bin * 100,
            yerr=alpha_diff_sigma_bin * 100,
            fmt=".",
            ms=8.0,
            color="b",
            lw=1.2,
            label="Diffusion (100x)",
        )

    if plot_viscosity:
        plot(x, viscosity, ".", color="g", ms=0.5, alpha=0.2)
        errorbar(
            x_bin,
            alpha_visc_bin,
            yerr=alpha_visc_sigma_bin,
            fmt=".",
            ms=8.0,
            color="y",
            lw=1.2,
            label="Viscosity",
        )

    ylabel("${\\rm{Rate~Coefficient}}~\\alpha$", labelpad=0)
    legend()
else:
    plot(x, S, ".", color="r", ms=0.5, alpha=0.2)
    plot(x_s, s_s, "--", color="k", alpha=0.8, lw=1.2)
    errorbar(x_bin, S_bin, yerr=S_sigma_bin, fmt=".", ms=8.0, color="b", lw=1.2)
    ylabel("${\\rm{Entropy}}~S$", labelpad=0)
    ylim(0.8, 3.8)

# Information -------------------------------------
subplot(236, frameon=False)

text_fontsize = 5

znow = 1.0 / anow - 1.0
text(
    -0.49,
    0.9,
    "Sod shock with $\\gamma=%.3f$ in 3D at $z=%.2f$" % (gas_gamma, znow),
    fontsize=text_fontsize,
)
text(
    -0.49,
    0.8,
    "Left: $(P_L, \\rho_L, v_L) = (%.3f, %.3f, %.3f)$" % (P_L, rho_L, v_L),
    fontsize=text_fontsize,
)
text(
    -0.49,
    0.7,
    "Right: $(P_R, \\rho_R, v_R) = (%.3f, %.3f, %.3f)$" % (P_R, rho_R, v_R),
    fontsize=text_fontsize,
)
z_i = 1.0 / a_i - 1.0
text(-0.49, 0.6, "Initial redshift: $%.2f$" % z_i, fontsize=text_fontsize)
plot([-0.49, 0.1], [0.52, 0.52], "k-", lw=1)
text(-0.49, 0.4, "SWIFT %s" % git, fontsize=text_fontsize)
text(-0.49, 0.3, scheme, fontsize=text_fontsize)
text(-0.49, 0.2, kernel, fontsize=text_fontsize)
text(
    -0.49,
    0.1,
    "$%.2f$ neighbours ($\\eta=%.3f$)" % (neighbours, eta),
    fontsize=text_fontsize,
)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])

tight_layout()
savefig("SodShock.png", dpi=200)
