###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

from h5py import File as HDF5File
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
import makeIC
import sys
sys.path.append("../../../logger/.libs/")
import liblogger as logger

# Plot parameters
plt.style.use("mpl_style")

center = np.array([0.5 * makeIC.boxsize]*3)
id_focus = 0

# Defines the constants
G = 1.189972e-04  # gravitational constant
M = 3.329460e+05  # mass of the sun

# generate the figures
fig_1 = plt.figure()
fig_2 = plt.figure()
fig_3 = plt.figure()


def gravity(t, y):
    """
    Compute the equation of motion.

    Parameters
    ----------

    t: float
      Time of the step
    y: np.array
      Variable to evolve [x, y, vx, vy].
    """
    dy = np.zeros(4)
    dy[:2] = y[2:]

    r = np.sum((y[:2] - center[:2])**2)**0.5

    a = - G * M / r**3
    dy[2] = y[2] * a
    dy[3] = y[3] * a
    return dy


def plotRelative(t, p, *args, **kwargs):
    """
    Wrapper around the function plot from matplotlib, but
    plot the relative evolution of the variable.
    """
    p = (p - p[0]) / p[0]
    plt.plot(t, p, *args, **kwargs)


def doSnapshots():
    """
    Read the snapshots and plot the corresponding variables.
    """

    # get all the filenames
    filenames = glob("simple_orbits_*.hdf5")
    N = len(filenames)
    filenames.sort()

    # generate the output arrays
    E = np.zeros((N, makeIC.num_part))
    t = np.zeros(N)
    p = np.zeros((N, 3))
    v = np.zeros((N, 3))

    for i, f in enumerate(filenames):
        # get the data from the file
        f = HDF5File(f, "r")
        ids = f["PartType1/ParticleIDs"][:]
        sort = np.argsort(ids)
        ids = ids[sort]
        pos = f["PartType1/Coordinates"][sort, :]
        pos -= center
        vel = f["PartType1/Velocities"][sort, :]

        t[i] = f["Header"].attrs["Time"]

        r = np.sum(pos**2, axis=1)**0.5
        v2 = np.sum(vel**2, axis=1)
        E[i, :] = 0.5 * v2 - G * M / r

        # Get the pos / vel of the required particle
        ind = ids == id_focus
        p[i, :] = pos[ind, :]
        v[i, :] = vel[ind, :]

    # Compute the solution
    y0 = np.zeros(4)
    y0[:2] = p[0, :2]
    y0[2:] = v[0, :2]

    # compute the plotting variables
    plt.figure(fig_1.number)
    plotRelative(t, E, ".", label="Snapshot")

    plt.figure(fig_2.number)
    plt.plot(p[:, 0], p[:, 1], "-", label="Snapshot", lw=1.)

    plt.figure(fig_3.number)
    plt.plot(v[:, 0], v[:, 1], "-", label="Snapshot", lw=1.)


def doStatistics():
    """
    Do the plots with the energy output file.
    """
    data = np.genfromtxt("energy.txt", names=True)

    times = data["Time"]
    E = data["E_tot"]
    plt.figure(fig_1.number)
    plotRelative(times, E, "-", label="Statistics")


def doLogger():
    """
    Read the logfile and plot the corresponding variables.
    """
    basename = "index_0000"
    N = 1000

    # Get time limits
    with logger.Reader(basename, verbose=0) as reader:
        t_min, t_max = reader.get_time_limits()
        times = np.linspace(t_min, t_max, N)

        # Create output arrays
        E = np.zeros((N, makeIC.num_part))
        p = np.zeros((N, 3))
        v = np.zeros((N, 3))

        for i, t in enumerate(times):
            # Get the next particles
            pos, vel, ids = reader.get_particle_data(
                ["Coordinates", "Velocities", "ParticleIDs"], t)
            sort = np.argsort(ids)
            ids = ids[sort]
            rel_pos = pos[sort, :] - center
            vel = vel[sort, :]

            # Compute the derived values
            r = np.sum(rel_pos**2, axis=1)**0.5
            v2 = np.sum(vel**2, axis=1)
            E[i, :] = 0.5 * v2 - G * M / r
            ind = ids == id_focus
            p[i, :] = rel_pos[ind, :]
            v[i, :] = vel[ind, :]

    # compute the plotting variables
    plt.figure(fig_1.number)
    plotRelative(times, E, "--", label="Logger (Interpolation)")

    # Compute the solution
    y0 = np.zeros(4)
    y0[:2] = p[0, :2]
    y0[2:] = v[0, :2]
    plt.figure(fig_2.number)
    plt.plot(p[:, 0], p[:, 1], ":r", label="Logger (Interpolation)")

    plt.figure(fig_3.number)
    plt.plot(v[:, 0], v[:, 1], ":r", label="Logger (Interpolation)")


# do all the plots
doStatistics()
doSnapshots()
doLogger()

# add text
plt.figure(fig_1.number)
plt.xlabel("Time [yr]")
plt.ylabel(r"$\frac{E - E(t=0)}{E(t=0)}$")
plt.legend(ncol=2)
plt.savefig("Energy.png")

plt.figure(fig_2.number)
plt.xlabel("Position [AU]")
plt.ylabel("Position [AU]")
plt.axis("equal")
plt.legend()
plt.savefig("Positions.pdf")

plt.figure(fig_3.number)
plt.xlabel("Velocity [AU / yr]")
plt.ylabel("Velocity [AU / yr]")
plt.axis("equal")
plt.legend()
plt.savefig("Velocities.png")

plt.show()
