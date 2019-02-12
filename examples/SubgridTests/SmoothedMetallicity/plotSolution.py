#!/usr/bin/env python
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Computes the analytical solution of the 3D Smoothed Metallicity example.

import h5py
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Parameters
low_metal = -6     # low metal abundance
high_metal = -5    # High metal abundance
sigma_metal = 0.1  # relative standard deviation for Z

Nelem = 9
# shift all metals in order to obtain nicer plots
low_metal = [low_metal] * Nelem + np.linspace(0, 3, Nelem)
high_metal = [high_metal] * Nelem + np.linspace(0, 3, Nelem)

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

# Plot parameters
params = {
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'font.size': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.figsize': (9.90, 6.45),
    'figure.subplot.left': 0.045,
    'figure.subplot.right': 0.99,
    'figure.subplot.bottom': 0.05,
    'figure.subplot.top': 0.99,
    'figure.subplot.wspace': 0.15,
    'figure.subplot.hspace': 0.12,
    'lines.markersize': 6,
    'lines.linewidth': 3.,
    'text.latex.unicode': True
}

plt.rcParams.update(params)
plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Times']})


snap = int(sys.argv[1])


# Read the simulation data
sim = h5py.File("smoothed_metallicity_%04d.hdf5" % snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
chemistry = sim["/SubgridScheme"].attrs["Chemistry Model"]
git = sim["Code"].attrs["Git Revision"]

pos = sim["/PartType0/Coordinates"][:, :]
d = pos[:, 0] - boxSize / 2
smooth_metal = sim["/PartType0/SmoothedElementAbundance"][:, :]
metal = sim["/PartType0/ElementAbundance"][:, :]
h = sim["/PartType0/SmoothingLength"][:]
h = np.mean(h)

if (Nelem != metal.shape[1]):
    print("Unexpected number of element, please check makeIC.py"
          " and plotSolution.py")
    exit(1)

N = 1000
d_a = np.linspace(-boxSize / 2., boxSize / 2., N)

# Now, work our the solution....


def calc_a(d, high_metal, low_metal, std_dev, h):
    """
    Compute analytical solution:
                        ___ High Metallicity
    Low Metallicity ___/

    where the linear part length is given by the smoothing length.

    Parameters
    ----------

    d: np.array
        Position to compute
    high_metal: float
        Value on the high metallicity side

    low_metal: float
        Value on the low metallicity side

    std_dev: float
        Standard deviation of the noise added

    h: float
        Mean smoothing length
    """

    # solution
    a = np.zeros([len(d), Nelem])
    # function at +- 1 sigma
    sigma = np.zeros([len(d), Nelem, 2])

    for i in range(Nelem):
        # compute low metallicity side
        s = d < -h
        a[s, i] = low_metal[i]
        # compute high metallicity side
        s = d > h
        a[s, i] = high_metal[i]

        # compute non constant parts
        m = (high_metal[i] - low_metal[i]) / (2.0 * h)
        c = (high_metal[i] + low_metal[i]) / 2.0
        # compute left linear part
        s = d < - boxSize / 2.0 + h
        a[s, i] = - m * (d[s] + boxSize / 2.0) + c
        # compute middle linear part
        s = np.logical_and(d >= -h, d <= h)
        a[s, i] = m * d[s] + c
        # compute right linear part
        s = d > boxSize / 2.0 - h
        a[s, i] = - m * (d[s] - boxSize / 2.0) + c

    sigma[:, :, 0] = a * (1 + std_dev)
    sigma[:, :, 1] = a * (1 - std_dev)
    return a, sigma


# compute solution
sol, sig = calc_a(d_a, high_metal, low_metal, sigma_metal, h)

# Plot the interesting quantities
plt.figure()

# Metallicity --------------------------------
plt.subplot(221)
for e in range(Nelem):
    plt.plot(metal[:, e], smooth_metal[:, e], '.', ms=0.5, alpha=0.2)

xmin, xmax = metal.min(), metal.max()
ymin, ymax = smooth_metal.min(), smooth_metal.max()
x = max(xmin, ymin)
y = min(xmax, ymax)
plt.plot([x, y], [x, y], "--k", lw=1.0)
plt.xlabel("${\\rm{Metallicity}}~Z_\\textrm{part}$", labelpad=0)
plt.ylabel("${\\rm{Smoothed~Metallicity}}~Z_\\textrm{sm}$", labelpad=0)

# Metallicity --------------------------------
e = 0
plt.subplot(223)
plt.plot(d, smooth_metal[:, e], '.', color='r', ms=0.5, alpha=0.2)
plt.plot(d_a, sol[:, e], '--', color='b', alpha=0.8, lw=1.2)
plt.fill_between(d_a, sig[:, e, 0], sig[:, e, 1], facecolor="b",
                 interpolate=True, alpha=0.5)
plt.xlabel("${\\rm{Distance}}~r$", labelpad=0)
plt.ylabel("${\\rm{Smoothed~Metallicity}}~Z_\\textrm{sm}$", labelpad=0)
plt.xlim(-0.5, 0.5)
plt.ylim(low_metal[e]-1, high_metal[e]+1)

# Information -------------------------------------
plt.subplot(222, frameon=False)

plt.text(-0.49, 0.9, "Smoothed Metallicity in 3D at $t=%.2f$" % time,
         fontsize=10)
plt.plot([-0.49, 0.1], [0.82, 0.82], 'k-', lw=1)
plt.text(-0.49, 0.7, "$\\textsc{Swift}$ %s" % git, fontsize=10)
plt.text(-0.49, 0.6, scheme, fontsize=10)
plt.text(-0.49, 0.5, kernel, fontsize=10)
plt.text(-0.49, 0.4, chemistry + "'s Chemistry", fontsize=10)
plt.text(-0.49, 0.3, "$%.2f$ neighbours ($\\eta=%.3f$)" % (neighbours, eta),
         fontsize=10)
plt.xlim(-0.5, 0.5)
plt.ylim(0, 1)
plt.xticks([])
plt.yticks([])

plt.tight_layout()
plt.savefig("SmoothedMetallicity.png", dpi=200)
