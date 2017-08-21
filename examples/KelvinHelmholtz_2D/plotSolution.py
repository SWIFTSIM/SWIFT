###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016  Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Computes the analytical solution of the Gresho-Chan vortex and plots the SPH answer

# Parameters
gas_gamma = 5./3.     # Gas adiabatic index
P1    = 2.5       # Central region pressure
P2    = 2.5       # Outskirts pressure
v1    = 0.5       # Central region velocity
v2    = -0.5      # Outskirts vlocity
rho1  = 2         # Central density
rho2  = 1         # Outskirts density

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

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
sim = h5py.File("kelvinHelmholtz_%04d.hdf5"%snap, "r")
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
v_norm = sqrt(vel[:,0]**2 + vel[:,1]**2)
rho = sim["/PartType0/Density"][:]
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]

# Plot the interesting quantities
figure()


# Azimuthal velocity profile -----------------------------
subplot(231)
scatter(pos[:,0], pos[:,1], c=vel[:,0], cmap="PuBu", edgecolors='face', s=4, vmin=-1., vmax=1.)
text(0.97, 0.97, "${\\rm{Velocity~along}}~x$", ha="right", va="top", backgroundcolor="w")
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Position}}~y$", labelpad=0)
xlim(0, 1)
ylim(0, 1)

# Radial density profile --------------------------------
subplot(232)
scatter(pos[:,0], pos[:,1], c=rho, cmap="PuBu", edgecolors='face', s=4, vmin=0.8, vmax=2.2)
text(0.97, 0.97, "${\\rm{Density}}$", ha="right", va="top", backgroundcolor="w")
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Position}}~y$", labelpad=0)
xlim(0, 1)
ylim(0, 1)

# Radial pressure profile --------------------------------
subplot(233)
scatter(pos[:,0], pos[:,1], c=P, cmap="PuBu", edgecolors='face', s=4, vmin=1, vmax=4)
text(0.97, 0.97, "${\\rm{Pressure}}$", ha="right", va="top", backgroundcolor="w")
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Position}}~y$", labelpad=0)
xlim(0, 1)
ylim(0, 1)

# Internal energy profile --------------------------------
subplot(234)
scatter(pos[:,0], pos[:,1], c=u, cmap="PuBu", edgecolors='face', s=4, vmin=1.5, vmax=5.)
text(0.97, 0.97, "${\\rm{Internal~energy}}$", ha="right", va="top", backgroundcolor="w")
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Position}}~y$", labelpad=0)
xlim(0, 1)
ylim(0, 1)

# Radial entropy profile --------------------------------
subplot(235)
scatter(pos[:,0], pos[:,1], c=S, cmap="PuBu", edgecolors='face', s=4, vmin=0.5, vmax=3.)
text(0.97, 0.97, "${\\rm{Entropy}}$", ha="right", va="top", backgroundcolor="w")
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Position}}~y$", labelpad=0)
xlim(0, 1)
ylim(0, 1)

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

text(-0.49, 0.9, "Kelvin-Helmholtz instability at $t=%.2f$"%(time), fontsize=10)
text(-0.49, 0.8, "Centre:~~~ $(P, \\rho, v) = (%.3f, %.3f, %.3f)$"%(P1, rho1, v1), fontsize=10)
text(-0.49, 0.7, "Outskirts: $(P, \\rho, v) = (%.3f, %.3f, %.3f)$"%(P2, rho2, v2), fontsize=10)
plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])

savefig("KelvinHelmholtz.png", dpi=200)
