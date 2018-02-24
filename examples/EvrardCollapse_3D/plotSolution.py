###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
################################################################################

# Compares the swift result for the 2D spherical Sod shock with a high
# resolution 2D reference result

import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py

# Parameters
gas_gamma = 5./3.      # Polytropic index
rho_L = 1.             # Density left state
rho_R = 0.125          # Density right state
v_L = 0.               # Velocity left state
v_R = 0.               # Velocity right state
P_L = 1.               # Pressure left state
P_R = 0.1              # Pressure right state

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

# Read the simulation data
sim = h5py.File("evrard_%04d.hdf5"%snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

coords = sim["/PartType0/Coordinates"]
x = sqrt((coords[:,0] - 0.5 * boxSize)**2 + (coords[:,1] - 0.5 * boxSize)**2 + \
         (coords[:,2] - 0.5 * boxSize)**2)
vels = sim["/PartType0/Velocities"]
v = sqrt(vels[:,0]**2 + vels[:,1]**2 + vels[:,2]**2)
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
rho = sim["/PartType0/Density"][:]

# Bin the data
x_bin_edge = logspace(-3., log10(2.), 100)
x_bin = 0.5*(x_bin_edge[1:] + x_bin_edge[:-1])
rho_bin,_,_ = stats.binned_statistic(x, rho, statistic='mean', bins=x_bin_edge)
v_bin,_,_ = stats.binned_statistic(x, v, statistic='mean', bins=x_bin_edge)
P_bin,_,_ = stats.binned_statistic(x, P, statistic='mean', bins=x_bin_edge)
S_bin,_,_ = stats.binned_statistic(x, S, statistic='mean', bins=x_bin_edge)
u_bin,_,_ = stats.binned_statistic(x, u, statistic='mean', bins=x_bin_edge)
rho2_bin,_,_ = stats.binned_statistic(x, rho**2, statistic='mean', bins=x_bin_edge)
v2_bin,_,_ = stats.binned_statistic(x, v**2, statistic='mean', bins=x_bin_edge)
P2_bin,_,_ = stats.binned_statistic(x, P**2, statistic='mean', bins=x_bin_edge)
S2_bin,_,_ = stats.binned_statistic(x, S**2, statistic='mean', bins=x_bin_edge)
u2_bin,_,_ = stats.binned_statistic(x, u**2, statistic='mean', bins=x_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin**2)
v_sigma_bin = np.sqrt(v2_bin - v_bin**2)
P_sigma_bin = np.sqrt(P2_bin - P_bin**2)
S_sigma_bin = np.sqrt(S2_bin - S_bin**2)
u_sigma_bin = np.sqrt(u2_bin - u_bin**2)

ref = loadtxt("evrardCollapse3D_exact.txt")

# Plot the interesting quantities
figure()

# Velocity profile --------------------------------
subplot(231)
semilogx(x, -v, '.', color='r', ms=0.2)
semilogx(ref[:,0], ref[:,2], "k--", alpha=0.8, lw=1.2)
errorbar(x_bin, -v_bin, yerr=v_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Velocity}}~v_r$", labelpad=0)
xlim(1.e-3, 2.)
ylim(-1.7, 0.1)

# Density profile --------------------------------
subplot(232)
loglog(x, rho, '.', color='r', ms=0.2)
loglog(ref[:,0], ref[:,1], "k--", alpha=0.8, lw=1.2)
errorbar(x_bin, rho_bin, yerr=rho_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho$", labelpad=0)
xlim(1.e-3, 2.)
ylim(1.e-2, 1.e4)

# Pressure profile --------------------------------
subplot(233)
loglog(x, P, '.', color='r', ms=0.2)
loglog(ref[:,0], ref[:,3], "k--", alpha=0.8, lw=1.2)
errorbar(x_bin, P_bin, yerr=P_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Pressure}}~P$", labelpad=0)
xlim(1.e-3, 2.)
ylim(1.e-4, 1.e3)

# Internal energy profile -------------------------
subplot(234)
loglog(x, u, '.', color='r', ms=0.2)
loglog(ref[:,0], ref[:,3] / ref[:,1] / (gas_gamma - 1.), "k--", alpha=0.8, lw=1.2)
errorbar(x_bin, u_bin, yerr=u_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
xlim(1.e-3, 2.)
ylim(1.e-2, 2.)

# Entropy profile ---------------------------------
subplot(235)
semilogx(x, S, '.', color='r', ms=0.2)
semilogx(ref[:,0], ref[:,3] / ref[:,1]**gas_gamma, "k--", alpha=0.8, lw=1.2)
errorbar(x_bin, S_bin, yerr=S_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Entropy}}~S$", labelpad=0)
xlim(1.e-3, 2.)
ylim(0., 0.25)

# Information -------------------------------------
subplot(236, frameon=False)

text(-0.49, 0.9, "Evrard collapse with $\\gamma=%.3f$ in 3D\nat $t=%.2f$"%(gas_gamma,time), fontsize=10)
plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])

tight_layout()
savefig("EvrardCollapse.png", dpi=200)
