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
rho0 = 1          # Gas density
P0 = 0.           # Constant additional pressure (should have no impact on the dynamics)

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
'figure.figsize' : (6.45,6.45),
'figure.subplot.left'    : 0.07,
'figure.subplot.right'   : 0.985,
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
sim = h5py.File("gresho_%03d.hdf5"%snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]

pos = sim["/PartType0/Coordinates"][:,:]
x = pos[:,0] - boxSize / 2
y = pos[:, 1] - boxSize / 2
vel = sim["/PartType0/Velocities"][:,:]
r = sqrt(x**2 + y**2)
v_r = (x * vel[:,0] + y * vel[:,1]) / r
v_phi = (-y * vel[:,0] + x * vel[:,1]) / r
v_norm = sqrt(vel[:,0]**2 + vel[:,1]**2)
u = sim["/PartType0/InternalEnergy"][:]

# Plot the interesting quantities
figure()

# Internal energy profile --------------------------------
subplot(221)

plot(r, u, 'x', color='r')
plot(solution_r, solution_u, '--', color='k', alpha=0.8, lw=1.2)
plot([0.2, 0.2], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
plot([0.4, 0.4], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
text(0.02, 8.97, "%s"%scheme, fontsize=9, backgroundcolor="w")
text(0.02, 8.82, "%s"%kernel, fontsize=9, backgroundcolor="w")
xlim(0,R_max)
ylim(7.3, 9.1)

# Azimuthal velocity profile -----------------------------
subplot(222)

plot(r, v_phi, 'x', color='r')
plot(solution_r, solution_v_phi, '--', color='k', alpha=0.8, lw=1.2)
plot([0.2, 0.2], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
plot([0.4, 0.4], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Azimuthal~velocity}}~v_\\phi$", labelpad=0)
xlim(0,R_max)
ylim(-0.1, 1.2)
text(0.75, 1, "$t=%.3f$"%(time), ha="right", va="bottom")

# Radial velocity profile --------------------------------
subplot(223)

plot(r, v_r, 'x', color='r', label="$\\texttt{SWIFT}$")
plot(solution_r, solution_v_r, '--', color='k', alpha=0.8, lw=1.2,  label="$\\rm{Solution}$")
plot([0.2, 0.2], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
plot([0.4, 0.4], [-100, 100], ':', color='k', alpha=0.4, lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Radial~velocity}}~v_r$", labelpad=-5)
legend(loc="upper right", fontsize=10, frameon=False, numpoints=1, handletextpad=0.1, handlelength=1.2)
xlim(0,R_max)
ylim(-0.2, 0.2)
yticks([-0.2, -0.1, 0., 0.1, 0.2])


# Image --------------------------------------------------
subplot(224)
scatter(pos[:,0], pos[:,1], c=v_norm, cmap="PuBu", edgecolors='face', s=4, vmin=0, vmax=1)
text(0.95, 0.95, "$|v|$", ha="right", va="top")
xlim(0,1)
ylim(0,1)
xlabel("$x$", labelpad=0)
ylabel("$y$", labelpad=0)

savefig("greshoVortex_%03d.png"%snap)
