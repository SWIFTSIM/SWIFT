################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
gas_gamma = 5./3. # Gas adiabatic index
rho0 = 1          # Gas density
P0 = 0.           # Constant additional pressure (should have no impact on the
                  # dynamics)

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py

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
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


snap = int(sys.argv[1])

# Generate the analytic solution at this time
N = 200
R_max = 0.8
solution_r = arange(0, R_max, R_max / N)
solution_P = zeros(N)
solution_v_phi = zeros(N)
solution_v_r = zeros(N)

for i in range(N):
    if solution_r[i] < 0.2:
        solution_P[i] = P0 + 5. + 12.5*solution_r[i]**2
        solution_v_phi[i] = 5.*solution_r[i]
    elif solution_r[i] < 0.4:
        solution_P[i] = P0 + 9. + 12.5*solution_r[i]**2 - 20.*solution_r[i] + 4.*log(solution_r[i]/0.2)
        solution_v_phi[i] = 2. -5.*solution_r[i]
    else:
        solution_P[i] = P0 + 3. + 4.*log(2.)
        solution_v_phi[i] = 0.

solution_rho = ones(N) * rho0
solution_s = solution_P / solution_rho**gas_gamma
solution_u = solution_P /((gas_gamma - 1.)*solution_rho)

# Read the simulation data
sim = h5py.File("gresho_%04d.hdf5"%snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

pos = sim["/PartType0/Coordinates"][:,:]
x = pos[:,0] - boxSize / 2
y = pos[:,1] - boxSize / 2
vel = sim["/PartType0/Velocities"][:,:]
r = sqrt(x**2 + y**2)
v_r = (x * vel[:,0] + y * vel[:,1]) / r
v_phi = (-y * vel[:,0] + x * vel[:,1]) / r
v_norm = sqrt(vel[:,0]**2 + vel[:,1]**2)
rho = sim["/PartType0/Density"][:]
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]

# Bin te data
r_bin_edge = np.arange(0., 1., 0.02)
r_bin = 0.5*(r_bin_edge[1:] + r_bin_edge[:-1])
rho_bin,_,_ = stats.binned_statistic(r, rho, statistic='mean', bins=r_bin_edge)
v_bin,_,_ = stats.binned_statistic(r, v_phi, statistic='mean', bins=r_bin_edge)
P_bin,_,_ = stats.binned_statistic(r, P, statistic='mean', bins=r_bin_edge)
S_bin,_,_ = stats.binned_statistic(r, S, statistic='mean', bins=r_bin_edge)
u_bin,_,_ = stats.binned_statistic(r, u, statistic='mean', bins=r_bin_edge)
rho2_bin,_,_ = stats.binned_statistic(r, rho**2, statistic='mean', bins=r_bin_edge)
v2_bin,_,_ = stats.binned_statistic(r, v_phi**2, statistic='mean', bins=r_bin_edge)
P2_bin,_,_ = stats.binned_statistic(r, P**2, statistic='mean', bins=r_bin_edge)
S2_bin,_,_ = stats.binned_statistic(r, S**2, statistic='mean', bins=r_bin_edge)
u2_bin,_,_ = stats.binned_statistic(r, u**2, statistic='mean', bins=r_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin**2)
v_sigma_bin = np.sqrt(v2_bin - v_bin**2)
P_sigma_bin = np.sqrt(P2_bin - P_bin**2)
S_sigma_bin = np.sqrt(S2_bin - S_bin**2)
u_sigma_bin = np.sqrt(u2_bin - u_bin**2)


# Plot the interesting quantities
figure()


# Azimuthal velocity profile -----------------------------
subplot(231)

plot(r, v_phi, '.', color='r', ms=0.5)
plot(solution_r, solution_v_phi, '--', color='k', alpha=0.8, lw=1.2)
errorbar(r_bin, v_bin, yerr=v_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
plot([0.2, 0.2], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
plot([0.4, 0.4], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Azimuthal~velocity}}~v_\\phi$", labelpad=0)
xlim(0,R_max)
ylim(-0.1, 1.2)

# Radial density profile --------------------------------
subplot(232)

plot(r, rho, '.', color='r', ms=0.5)
plot(solution_r, solution_rho, '--', color='k', alpha=0.8, lw=1.2)
errorbar(r_bin, rho_bin, yerr=rho_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
plot([0.2, 0.2], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
plot([0.4, 0.4], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho$", labelpad=0)
xlim(0,R_max)
ylim(rho0-0.3, rho0 + 0.3)
#yticks([-0.2, -0.1, 0., 0.1, 0.2])

# Radial pressure profile --------------------------------
subplot(233)

plot(r, P, '.', color='r', ms=0.5)
plot(solution_r, solution_P, '--', color='k', alpha=0.8, lw=1.2)
errorbar(r_bin, P_bin, yerr=P_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
plot([0.2, 0.2], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
plot([0.4, 0.4], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Pressure}}~P$", labelpad=0)
xlim(0, R_max)
ylim(4.9 + P0, P0 + 6.1)

# Internal energy profile --------------------------------
subplot(234)

plot(r, u, '.', color='r', ms=0.5)
plot(solution_r, solution_u, '--', color='k', alpha=0.8, lw=1.2)
errorbar(r_bin, u_bin, yerr=u_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
plot([0.2, 0.2], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
plot([0.4, 0.4], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
xlim(0,R_max)
ylim(7.3, 9.1)


# Radial entropy profile --------------------------------
subplot(235)

plot(r, S, '.', color='r', ms=0.5)
plot(solution_r, solution_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(r_bin, S_bin, yerr=S_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
plot([0.2, 0.2], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
plot([0.4, 0.4], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Entropy}}~S$", labelpad=0)
xlim(0, R_max)
ylim(4.9 + P0, P0 + 6.1)

# Image --------------------------------------------------
#subplot(234)
#scatter(pos[:,0], pos[:,1], c=v_norm, cmap="PuBu", edgecolors='face', s=4, vmin=0, vmax=1)
#text(0.95, 0.95, "$|v|$", ha="right", va="top")
#xlim(0,1)
#ylim(0,1)
#xlabel("$x$", labelpad=0)
#ylabel("$y$", labelpad=0)

# Information -------------------------------------
subplot(236, frameon=False)

text(-0.49, 0.9, "Gresho-Chan vortex with  $\\gamma=%.3f$ at $t=%.2f$"%(gas_gamma,time), fontsize=10)
text(-0.49, 0.8, "Background $\\rho_0=%.3f$"%rho0, fontsize=10)
text(-0.49, 0.7, "Background $P_0=%.3f$"%P0, fontsize=10)
plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])

savefig("GreshoVortex.png", dpi=200)
