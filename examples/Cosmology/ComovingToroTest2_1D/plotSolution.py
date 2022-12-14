###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Computes the analytical solution of the Toro (2009) test 2 and plots the SPH answer

# Generates the analytical solution for the Toro (2009) test case 2
# The script works for a given left (x<0) and right (x>0) state and computes the solution at a later time t.
# This follows the solution given in (Toro, 2009)

# Parameters
gas_gamma = 5.0 / 3.0  # Polytropic index
rho_L = 1.0  # Density left state
rho_R = 1.0  # Density right state
v_L = -2.0  # Velocity left state
v_R = 2.0  # Velocity right state
P_L = 0.4  # Pressure left state
P_R = 0.4  # Pressure right state

import sys

sys.path.append("../../HydroTests/")
from riemannSolver import RiemannSolver

import matplotlib

matplotlib.use("Agg")
from pylab import *
import h5py

style.use("../../../tools/stylesheets/mnras.mplstyle")

snap = int(sys.argv[1])


# Read the simulation data
sim = h5py.File("toroTest2_%04d.hdf5" % snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
anow = sim["/Header"].attrs["Scale-factor"]
a_i = sim["/Cosmology"].attrs["a_beg"]
H_0 = sim["/Cosmology"].attrs["H0 [internal units]"]
time = 2.0 * (1.0 / np.sqrt(a_i) - 1.0 / np.sqrt(anow)) / H_0
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

x = sim["/PartType0/Coordinates"][:, 0]
v = sim["/PartType0/Velocities"][:, 0] * anow
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]
rho = sim["/PartType0/Densities"][:]

N = 1000  # Number of points
x_min = -1.0
x_max = 1.0

x += x_min

# Prepare reference solution
solver = RiemannSolver(gas_gamma)

delta_x = (x_max - x_min) / N
x_s = arange(0.5 * x_min, 0.5 * x_max, delta_x)
rho_s, v_s, P_s, _ = solver.solve(rho_L, v_L, P_L, rho_R, v_R, P_R, x_s / time)
rho_s2, v_s2, P_s2, _ = solver.solve(rho_R, v_R, P_R, rho_L, v_L, P_L, x_s / time)
x_s2 = np.array(x_s)
x_s2 += 1.0
s2neg = x_s2 > 1.0
s2pos = ~s2neg
x_s2[s2neg] -= 2.0

# Additional arrays
u_s = P_s / (rho_s * (gas_gamma - 1.0))  # internal energy
s_s = P_s / rho_s ** gas_gamma  # entropic function
u_s2 = P_s2 / (rho_s2 * (gas_gamma - 1.0))  # internal energy
s_s2 = P_s2 / rho_s2 ** gas_gamma  # entropic function


# Plot the interesting quantities
figure(figsize=(7, 7 / 1.6))

# Velocity profile --------------------------------
subplot(231)
plot(x, v, ".", color="r", ms=4.0)
plot(x_s, v_s, "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2pos], v_s2[s2pos], "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2neg], v_s2[s2neg], "--", color="k", alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Velocity}}~v_x$", labelpad=0)

# Density profile --------------------------------
subplot(232)
plot(x, rho, ".", color="r", ms=4.0)
plot(x_s, rho_s, "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2pos], rho_s2[s2pos], "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2neg], rho_s2[s2neg], "--", color="k", alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho$", labelpad=0)

# Pressure profile --------------------------------
subplot(233)
plot(x, P, ".", color="r", ms=4.0)
plot(x_s, P_s, "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2pos], P_s2[s2pos], "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2neg], P_s2[s2neg], "--", color="k", alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Pressure}}~P$", labelpad=0)

# Internal energy profile -------------------------
subplot(234)
plot(x, u, ".", color="r", ms=4.0)
plot(x_s, u_s, "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2pos], u_s2[s2pos], "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2neg], u_s2[s2neg], "--", color="k", alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)

# Entropy profile ---------------------------------
subplot(235)

plot(x, S, ".", color="r", ms=4.0)
plot(x_s, s_s, "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2pos], s_s2[s2pos], "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2neg], s_s2[s2neg], "--", color="k", alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Entropy}}~S$", labelpad=0)


# Information -------------------------------------
subplot(236, frameon=False)

text_fontsize = 5

z_now = 1.0 / anow - 1.0
text(
    -0.49,
    0.9,
    "Toro (2009) test 2 with  $\\gamma=%.3f$ in 1D at $z=%.2f$" % (gas_gamma, z_now),
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
plot([-0.49, 0.1], [0.62, 0.62], "k-", lw=1)
text(-0.49, 0.5, "SWIFT %s" % git.decode("utf-8"), fontsize=text_fontsize)
text(-0.49, 0.4, scheme.decode("utf-8"), fontsize=text_fontsize)
text(-0.49, 0.3, kernel.decode("utf-8"), fontsize=text_fontsize)
text(
    -0.49,
    0.2,
    "$%.2f$ neighbours ($\\eta=%.3f$)" % (neighbours, eta),
    fontsize=text_fontsize,
)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])

tight_layout()

savefig("ToroTest2.png", dpi=200)
