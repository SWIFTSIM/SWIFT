###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Computes the analytical solution of the Noh problem and plots the SPH answer
 

# Parameters
gas_gamma = 5./3.      # Polytropic index
rho0 = 1.              # Background density
P0 = 1.e-6             # Background pressure
v0 = 1

import matplotlib
matplotlib.use("Agg")
from pylab import *
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


# Read the simulation data
sim = h5py.File("noh_%04d.hdf5"%snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

x = sim["/PartType0/Coordinates"][:,0]
v = sim["/PartType0/Velocities"][:,0]
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
rho = sim["/PartType0/Density"][:]

N = 1001  # Number of points
x_min = -1
x_max = 1

x += x_min

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------
x_s = np.arange(0, 2., 2./N) - 1.
rho_s = np.ones(N) * rho0
P_s = np.ones(N) * rho0
v_s = np.ones(N) * v0

# Shock position
u0 = rho0 * P0 * (gas_gamma-1)
us = 0.5 * (gas_gamma - 1) * v0
rs = us * time

# Post-shock values
rho_s[np.abs(x_s) < rs] = rho0 * (gas_gamma + 1) / (gas_gamma - 1)
P_s[np.abs(x_s) < rs] = 0.5 * rho0 * v0**2 * (gas_gamma + 1) 
v_s[np.abs(x_s) < rs] = 0.

# Pre-shock values
rho_s[np.abs(x_s) >= rs] = rho0
P_s[np.abs(x_s) >= rs] = 0
v_s[x_s >= rs] = -v0
v_s[x_s <= -rs] = v0

# Additional arrays
u_s = P_s / (rho_s * (gas_gamma - 1.))  #internal energy
s_s = P_s / rho_s**gas_gamma # entropic function

# Plot the interesting quantities
figure()

# Velocity profile --------------------------------
subplot(231)
plot(x, v, '.', color='r', ms=4.0)
plot(x_s, v_s, '--', color='k', alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Velocity}}~v_x$", labelpad=-4)
xlim(-0.5, 0.5)
ylim(-1.2, 1.2)

# Density profile --------------------------------
subplot(232)
plot(x, rho, '.', color='r', ms=4.0)
plot(x_s, rho_s, '--', color='k', alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho$", labelpad=0)
xlim(-0.5, 0.5)
ylim(0.95, 4.4)

# Pressure profile --------------------------------
subplot(233)
plot(x, P, '.', color='r', ms=4.0)
plot(x_s, P_s, '--', color='k', alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Pressure}}~P$", labelpad=0)
xlim(-0.5, 0.5)
ylim(-0.05, 1.8)

# Internal energy profile -------------------------
subplot(234)
plot(x, u, '.', color='r', ms=4.0)
plot(x_s, u_s, '--', color='k', alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
xlim(-0.5, 0.5)
ylim(-0.05, 0.8)

# Entropy profile ---------------------------------
subplot(235)
plot(x, S, '.', color='r', ms=4.0)
plot(x_s, s_s, '--', color='k', alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Entropy}}~S$", labelpad=-9)
xlim(-0.5, 0.5)
ylim(-0.05, 0.2)

# Information -------------------------------------
subplot(236, frameon=False)

text(-0.49, 0.9, "Noh problem with  $\\gamma=%.3f$ in 1D at $t=%.2f$"%(gas_gamma,time), fontsize=10)
text(-0.49, 0.8, "ICs:~~ $(P_0, \\rho_0, v_0) = (%1.2e, %.3f, %.3f)$"%(1e-6, 1., -1.), fontsize=10)
plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])


savefig("Noh.png", dpi=200)
