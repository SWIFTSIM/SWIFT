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

# Parameters
gamma = 1.4 # Polytropic index

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
file = h5py.File("interactingBlastWaves_{0:04d}.hdf5".format(snap), "r")
x = file["/PartType0/Coordinates"][:,0]
rho = file["/PartType0/Density"]
v = file["/PartType0/Velocities"][:,0]
u = file["/PartType0/InternalEnergy"]
S = file["/PartType0/Entropy"]
P = file["/PartType0/Pressure"]
time = file["/Header"].attrs["Time"][0]

scheme = file["/HydroScheme"].attrs["Scheme"]
kernel = file["/HydroScheme"].attrs["Kernel function"]
neighbours = file["/HydroScheme"].attrs["Kernel target N_ngb"][0]
eta = file["/HydroScheme"].attrs["Kernel eta"][0]
git = file["Code"].attrs["Git Revision"]

ref = np.loadtxt("interactingBlastWaves1D_exact.txt")

# Plot the interesting quantities
fig, ax = pl.subplots(2, 3)

# Velocity profile
ax[0][0].plot(x, v, "r.", markersize = 4.)
ax[0][0].plot(ref[:,0], ref[:,2], "k--", alpha = 0.8, linewidth = 1.2)
ax[0][0].set_xlabel("${\\rm{Position}}~x$", labelpad = 0)
ax[0][0].set_ylabel("${\\rm{Velocity}}~v_x$", labelpad = 0)
ax[0][0].set_xlim(0., 1.)
ax[0][0].set_ylim(-1., 15.)

# Density profile
ax[0][1].plot(x, rho, "r.", markersize = 4.)
ax[0][1].plot(ref[:,0], ref[:,1], "k--", alpha = 0.8, linewidth = 1.2)
ax[0][1].set_xlabel("${\\rm{Position}}~x$", labelpad = 0)
ax[0][1].set_ylabel("${\\rm{Density}}~\\rho$", labelpad = 0)
ax[0][1].set_xlim(0., 1.)

# Pressure profile
ax[0][2].plot(x, P, "r.", markersize = 4.)
ax[0][2].plot(ref[:,0], ref[:,3], "k--", alpha = 0.8, linewidth = 1.2)
ax[0][2].set_xlabel("${\\rm{Position}}~x$", labelpad = 0)
ax[0][2].set_ylabel("${\\rm{Pressure}}~P$", labelpad = 0)
ax[0][2].set_xlim(0., 1.)

# Internal energy profile
ax[1][0].plot(x, u, "r.", markersize = 4.)
ax[1][0].plot(ref[:,0], ref[:,3] / ref[:,1] / (gamma - 1.), "k--", alpha = 0.8,
              linewidth = 1.2)
ax[1][0].set_xlabel("${\\rm{Position}}~x$", labelpad = 0)
ax[1][0].set_ylabel("${\\rm{Internal~Energy}}~u$", labelpad = 0)
ax[1][0].set_xlim(0., 1.)

# Entropy profile
ax[1][1].plot(x, S, "r.", markersize = 4.)
ax[1][1].plot(ref[:,0], ref[:,3] / ref[:,1]**gamma, "k--", alpha = 0.8,
              linewidth = 1.2)
ax[1][1].set_xlabel("${\\rm{Position}}~x$", labelpad = 0)
ax[1][1].set_ylabel("${\\rm{Entropy}}~S$", labelpad = 0)
ax[1][1].set_xlim(0., 1.)

# Run information
ax[1][2].set_frame_on(False)
ax[1][2].text(-0.49, 0.9,
  "Interacting blast waves test\nwith $\\gamma={0:.3f}$ in 1D at $t = {1:.2f}$".format(
    gamma, time), fontsize = 10)
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
pl.savefig("InteractingBlastWaves.png", dpi = 200)
