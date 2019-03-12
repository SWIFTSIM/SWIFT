###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

import numpy as np
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import h5py
import sys
import scipy.stats as stats

# Parameters
gamma = 5. / 3. # Polytropic index
rhoL = 1.       # Initial density in the non vacuum state
vL = 0.         # Initial velocity in the non vacuum state
PL = 1.         # Initial pressure in the non vacuum state
rhoR = 0.       # Initial vacuum density
vR = 0.         # Initial vacuum velocity
PR = 0.         # Initial vacuum pressure

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (9.90,6.45),
'figure.subplot.left'    : 0.045,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.05,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
pl.rcParams.update(params)
pl.rc('font',**{'family':'sans-serif','sans-serif':['Times']})

# Read the snapshot index from the command line argument
snap = int(sys.argv[1])

# Open the file and read the relevant data
file = h5py.File("vacuum_{0:04d}.hdf5".format(snap), "r")
coords = file["/PartType0/Coordinates"]
x = np.sqrt((coords[:,0] - 0.5)**2 + (coords[:,1] - 0.5)**2 + \
            (coords[:,2] - 0.5)**2)
rho = file["/PartType0/Density"][:]
vels = file["/PartType0/Velocities"]
v = np.sqrt(vels[:,0]**2 + vels[:,1]**2 + vels[:,2]**2)
u = file["/PartType0/InternalEnergy"][:]
S = file["/PartType0/Entropy"][:]
P = file["/PartType0/Pressure"][:]
time = file["/Header"].attrs["Time"][0]

scheme = file["/HydroScheme"].attrs["Scheme"]
kernel = file["/HydroScheme"].attrs["Kernel function"]
neighbours = file["/HydroScheme"].attrs["Kernel target N_ngb"][0]
eta = file["/HydroScheme"].attrs["Kernel eta"][0]
git = file["Code"].attrs["Git Revision"]

# Bin the data values
# We let scipy choose the bins and then reuse them for all other quantities
rho_bin, x_bin_edge, _ = \
  stats.binned_statistic(x, rho, statistic = "mean", bins = 50)
rho2_bin, _, _ = \
  stats.binned_statistic(x, rho**2, statistic = "mean", bins = x_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin**2)

v_bin, _, _ = \
  stats.binned_statistic(x, v, statistic = "mean", bins = x_bin_edge)
v2_bin, _, _ = \
  stats.binned_statistic(x, v**2, statistic = "mean", bins = x_bin_edge)
v_sigma_bin = np.sqrt(v2_bin - v_bin**2)

P_bin, _, _ = \
  stats.binned_statistic(x, P, statistic = "mean", bins = x_bin_edge)
P2_bin, _, _ = \
  stats.binned_statistic(x, P**2, statistic = "mean", bins = x_bin_edge)
P_sigma_bin = np.sqrt(P2_bin - P_bin**2)

u_bin, _, _ = \
  stats.binned_statistic(x, u, statistic = "mean", bins = x_bin_edge)
u2_bin, _, _ = \
  stats.binned_statistic(x, u**2, statistic = "mean", bins = x_bin_edge)
u_sigma_bin = np.sqrt(u2_bin - u_bin**2)

S_bin, _, _ = \
  stats.binned_statistic(x, S, statistic = "mean", bins = x_bin_edge)
S2_bin, _, _ = \
  stats.binned_statistic(x, S**2, statistic = "mean", bins = x_bin_edge)
S_sigma_bin = np.sqrt(S2_bin - S_bin**2)

x_bin = 0.5 * (x_bin_edge[1:] + x_bin_edge[:-1])

ref = np.loadtxt("vacuumSpherical3D_exact.txt")

# Plot the interesting quantities
fig, ax = pl.subplots(2, 3)

# Velocity profile
ax[0][0].plot(x, v, "r.", markersize = 0.2)
ax[0][0].plot(ref[:,0], ref[:,2], "k--", alpha = 0.8, linewidth = 1.2)
ax[0][0].errorbar(x_bin, v_bin, yerr = v_sigma_bin, fmt = ".",
  markersize = 8., color = "b", linewidth = 1.2)
ax[0][0].set_xlabel("${\\rm{Radius}}~r$", labelpad = 0)
ax[0][0].set_ylabel("${\\rm{Velocity}}~v_r$", labelpad = 0)
ax[0][0].set_xlim(0., 0.4)
ax[0][0].set_ylim(-0.1, 3.2)

# Density profile
ax[0][1].plot(x, rho, "r.", markersize = 0.2)
ax[0][1].plot(ref[:,0], ref[:,1], "k--", alpha = 0.8, linewidth = 1.2)
ax[0][1].errorbar(x_bin, rho_bin, yerr = rho_sigma_bin, fmt = ".",
  markersize = 8., color = "b", linewidth = 1.2)
ax[0][1].set_xlabel("${\\rm{Radius}}~r$", labelpad = 0)
ax[0][1].set_ylabel("${\\rm{Density}}~\\rho$", labelpad = 0)
ax[0][1].set_xlim(0., 0.4)

# Pressure profile
ax[0][2].plot(x, P, "r.", markersize = 0.2)
ax[0][2].plot(ref[:,0], ref[:,3], "k--", alpha = 0.8, linewidth = 1.2)
ax[0][2].errorbar(x_bin, P_bin, yerr = P_sigma_bin, fmt = ".",
  markersize = 8., color = "b", linewidth = 1.2)
ax[0][2].set_xlabel("${\\rm{Radius}}~r$", labelpad = 0)
ax[0][2].set_ylabel("${\\rm{Pressure}}~P$", labelpad = 0)
ax[0][2].set_xlim(0., 0.4)

# Internal energy profile
ax[1][0].plot(x, u, "r.", markersize = 0.2)
ax[1][0].plot(ref[:,0], ref[:,3] / ref[:,1] / (gamma - 1.), "k--", alpha = 0.8,
              linewidth = 1.2)
ax[1][0].errorbar(x_bin, u_bin, yerr = u_sigma_bin, fmt = ".",
  markersize = 8., color = "b", linewidth = 1.2)
ax[1][0].set_xlabel("${\\rm{Radius}}~r$", labelpad = 0)
ax[1][0].set_ylabel("${\\rm{Internal~Energy}}~u$", labelpad = 0)
ax[1][0].set_xlim(0., 0.4)

# Entropy profile
ax[1][1].plot(x, S, "r.", markersize = 0.2)
ax[1][1].plot(ref[:,0], ref[:,3] / ref[:,1]**gamma, "k--", alpha = 0.8,
              linewidth = 1.2)
ax[1][1].errorbar(x_bin, S_bin, yerr = S_sigma_bin, fmt = ".",
  markersize = 8., color = "b", linewidth = 1.2)
ax[1][1].set_xlabel("${\\rm{Radius}}~r$", labelpad = 0)
ax[1][1].set_ylabel("${\\rm{Entropy}}~S$", labelpad = 0)
ax[1][1].set_xlim(0., 0.4)
ax[1][1].set_ylim(0., 4.)

# Run information
ax[1][2].set_frame_on(False)
ax[1][2].text(-0.49, 0.9,
  "Vacuum test with $\\gamma={0:.3f}$ in 1D at $t = {1:.2f}$".format(
    gamma, time), fontsize = 10)
ax[1][2].text(-0.49, 0.8,
  "Left:~~ $(P_L, \\rho_L, v_L) = ({0:.3f}, {1:.3f}, {2:.3f})$".format(
    PL, rhoL, vL), fontsize = 10)
ax[1][2].text(-0.49, 0.7,
  "Right: $(P_R, \\rho_R, v_R) = ({0:.3f}, {1:.3f}, {2:.3f})$".format(
    PR, rhoR, vR), fontsize = 10)
ax[1][2].plot([-0.49, 0.1], [0.62, 0.62], "k-", lw = 1)
ax[1][2].text(-0.49, 0.5, "$\\textsc{{Swift}}$ {0}".format(git), fontsize = 10)
ax[1][2].text(-0.49, 0.4, scheme, fontsize = 10)
ax[1][2].text(-0.49, 0.3, kernel, fontsize = 10)
ax[1][2].text(-0.49, 0.2,
  "${0:.2f}$ neighbours ($\\eta={1:.3f}$)".format(neighbours, eta),
  fontsize = 10)
ax[1][2].set_xlim(-0.5, 0.5)
ax[1][2].set_ylim(0., 1.)
ax[1][2].set_xticks([])
ax[1][2].set_yticks([])

pl.tight_layout()
pl.savefig("Vacuum.png", dpi = 200)
