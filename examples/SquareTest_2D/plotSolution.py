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

# Computes the analytical solution of the square test

# Parameters
gas_gamma = 5./3.     # Gas adiabatic index
gamma = 5./3.     # Gas adiabatic index
rho0 = 4          # Gas central density
rho1 = 1          # Gas outskirt density
P0 = 2.5          # Gas central pressure
P1 = 2.5          # Gas central pressure
vx = 142.3        # Random velocity for all particles 
vy = -31.4

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
sim = h5py.File("square_%04d.hdf5"%snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

# Analytical soltion
centre_x = 0.5 + time * vx
centre_y = 0.5 + time * vy

while centre_x > 1.:
    centre_x -= 1.
while centre_x < 0.:
    centre_x += 1.
while centre_y > 1.:
    centre_y -= 1.
while centre_y < 0.:
    centre_y += 1.

pos = sim["/PartType0/Coordinates"][:,:]
vel = sim["/PartType0/Velocities"][:,:]
v_norm = sqrt(vel[:,0]**2 + vel[:,1]**2)
rho = sim["/PartType0/Density"][:]
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
x = pos[:,0] - centre_x
y = pos[:,1] - centre_y

# Box wrapping
x[x>0.5] -= 1.
x[x<-0.5] += 1.
y[y>0.5] -= 1.
y[y<-0.5] += 1.

# Azimuthal velocity profile -----------------------------
subplot(231)
scatter(x, y, c=v_norm, cmap="PuBu", edgecolors='face', s=4, vmin=-1., vmax=1.)
text(0.47, 0.47, "${\\rm{Velocity~norm}}$", ha="right", va="top", backgroundcolor="w")
plot([-0.25, 0.25], [0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
plot([-0.25, 0.25], [-0.25, -0.25], '--', color='k', alpha=0.8, lw=2)
plot([0.25, 0.25], [-0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
plot([-0.25, -0.25], [-0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Position}}~y$", labelpad=-7)
xlim(-0.5, 0.5)
ylim(-0.5, 0.5)

# Radial density profile --------------------------------
subplot(232)
scatter(x, y, c=rho, cmap="PuBu", edgecolors='face', s=4, vmin=0., vmax=4.)
text(0.47, 0.47, "${\\rm{Density}}$", ha="right", va="top", backgroundcolor="w")
plot([-0.25, 0.25], [0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
plot([-0.25, 0.25], [-0.25, -0.25], '--', color='k', alpha=0.8, lw=2)
plot([0.25, 0.25], [-0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
plot([-0.25, -0.25], [-0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Position}}~y$", labelpad=-7)
xlim(-0.5, 0.5)
ylim(-0.5, 0.5)

# Radial pressure profile --------------------------------
subplot(233)
scatter(x, y, c=P, cmap="PuBu", edgecolors='face', s=4, vmin=2, vmax=4)
text(0.47, 0.47, "${\\rm{Pressure}}$", ha="right", va="top", backgroundcolor="w")
plot([-0.25, 0.25], [0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
plot([-0.25, 0.25], [-0.25, -0.25], '--', color='k', alpha=0.8, lw=2)
plot([0.25, 0.25], [-0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
plot([-0.25, -0.25], [-0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Position}}~y$", labelpad=-7)
xlim(-0.5, 0.5)
ylim(-0.5, 0.5)

# Internal energy profile --------------------------------
subplot(234)
scatter(x, y, c=u, cmap="PuBu", edgecolors='face', s=4, vmin=0.5, vmax=4.)
text(0.47, 0.47, "${\\rm{Internal~energy}}$", ha="right", va="top", backgroundcolor="w")
plot([-0.25, 0.25], [0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
plot([-0.25, 0.25], [-0.25, -0.25], '--', color='k', alpha=0.8, lw=2)
plot([0.25, 0.25], [-0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
plot([-0.25, -0.25], [-0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Position}}~y$", labelpad=-7)
xlim(-0.5, 0.5)
ylim(-0.5, 0.5)

# Radial entropy profile --------------------------------
subplot(235)
scatter(x, y, c=S, cmap="PuBu", edgecolors='face', s=4, vmin=0., vmax=3.)
text(0.47, 0.47, "${\\rm{Entropy}}$", ha="right", va="top", backgroundcolor="w")
plot([-0.25, 0.25], [0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
plot([-0.25, 0.25], [-0.25, -0.25], '--', color='k', alpha=0.8, lw=2)
plot([0.25, 0.25], [-0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
plot([-0.25, -0.25], [-0.25, 0.25], '--', color='k', alpha=0.8, lw=2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Position}}~y$", labelpad=-7)
xlim(-0.5, 0.5)
ylim(-0.5, 0.5)


# Information -------------------------------------
subplot(236, frameon=False)

text(-0.49, 0.9, "Square test with $\\gamma=%.3f$ at $t=%.2f$"%(gas_gamma,time), fontsize=10)
text(-0.49, 0.8, "Centre:~~~ $(P, \\rho) = (%.3f, %.3f)$"%(P0, rho0), fontsize=10)
text(-0.49, 0.7, "Outskirts: $(P, \\rho) = (%.3f, %.3f)$"%(P1, rho1), fontsize=10)
plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])

savefig("SquareTest.png", dpi=200)
