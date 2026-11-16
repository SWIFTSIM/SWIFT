###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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

import h5py
import numpy as np
import unyt
import swiftsimio as sw
from scipy.optimize import bisect

# Generates a swift IC file for the Rayleigh-Taylor instability in a periodic
# box following the conditions given in Saitoh and Makino 2013: 1202.4277v3

# Parameters
N = [128, 192]  # Particles along one edge
gamma = 7.0 / 5.0  # Gas adiabatic index
dv = 0.025  # velocity perturbation
rho_h = 2  # high density region
rho_l = 1  # Low density region
g = -0.5  # gravitational acceleration
box_size = [1.0, 1.5]  # size of the box

fixed = [0.1, 1.4]  # y-range of non fixed particles
perturbation_box = [0.3, 1.2]  # y-range for the velocity perturbation
fileOutputName = "rayleigh_taylor.hdf5"


# ---------------------------------------------------

if N[0] / box_size[0] != N[1] / box_size[1]:
    raise Exception("Suppose the same ratio for each directions")


def density(y):
    """
    Mass density as function of the position y.
    """
    if isinstance(y, float) or isinstance(y, int):
        y = np.array(y)

    ind = y < 0.5 * box_size[1]
    rho = np.zeros(y.shape)
    tmp = (gamma - 1.0) * g / (gamma * P0)
    alpha = 1.0 / (gamma - 1.0)

    rho[ind] = rho_l * (1 + rho_l * tmp * (y[ind] - 0.5 * box_size[1])) ** alpha
    rho[~ind] = rho_h * (1 + rho_h * tmp * (y[~ind] - 0.5 * box_size[1])) ** alpha

    return rho


def mass(y):
    """
    Integral of the density
    """
    if isinstance(y, float) or isinstance(y, int):
        y = np.array(y)

    ind = y < 0.5 * box_size[1]
    m = np.zeros(y.shape)

    B = (gamma - 1.0) * g / (gamma * P0)
    alpha = 1.0 / (gamma - 1.0)

    m[ind] = (1 + B * rho_l * (y[ind] - 0.5 * box_size[1])) ** (alpha + 1)

    m[~ind] = (1 + B * rho_h * (y[~ind] - 0.5 * box_size[1])) ** (alpha + 1)

    m -= (1 - 0.5 * B * box_size[1] * rho_l) ** (alpha + 1)

    return box_size[0] * m / (B * (alpha + 1))


P0 = rho_h / gamma  # integration constant of the hydrostatic equation
numPart = N[0] * N[1]

m_tot = mass(box_size[1])
m_part = m_tot / numPart


def inv_mass(m):
    """
    Inverse of the function `mass`.
    """
    left = 0
    right = box_size[1]

    def f(y, x):
        return x - mass(y)

    return bisect(f, left, right, args=(m))


def entropy(y):
    """
    Entropy as function of the position y.
    Here we assume isoentropic medium.
    """
    if isinstance(y, float) or isinstance(y, int):
        y = np.array(y)

    ind = y < 0.5 * box_size[1]

    a = np.zeros(y.shape)
    a[ind] = P0 * rho_l ** (-gamma)
    a[~ind] = P0 * rho_h ** (-gamma)
    return a


def smoothing_length(rho, m):
    """
    Compute the smoothing length
    """
    return 1.23 * np.sqrt(m / rho)


def growth_rate():
    """
    Compute the growth rate of the instability.
    Assumes a wavelength equal to the boxsize.
    """
    ymin = 0.0
    ymax = box_size[1]
    A = density(ymax) - density(ymin)
    A /= density(ymax) + density(ymin)
    return np.sqrt(A * np.abs(g) * ymax / (2.0 * np.pi))


def vy(x, y):
    """
    Velocity along the direction y
    """
    ind = y > perturbation_box[0]
    ind = np.logical_and(ind, y < perturbation_box[1])

    v = np.zeros(len(x))

    v[ind] = 1 + np.cos(4 * np.pi * x[ind])
    v[ind] *= 1 + np.cos(5 * np.pi * (y[ind] - 0.5 * box_size[1]))
    v[ind] *= dv
    return v


if __name__ == "__main__":
    # Start by generating grids of particles

    coords = np.zeros((numPart, 3))
    m = np.ones(numPart) * m_part
    u = np.zeros(numPart)
    vel = np.zeros((numPart, 3))
    ids = np.zeros(numPart)

    # generate grid of particles
    y_prev = 0
    uni_id = 0

    # loop over y
    eps = 1e-3 * box_size[1] / N[1]
    for j in range(N[1]):
        m_y = m_part * N[0] + mass(y_prev)
        if m_y > m_tot:
            y_j = box_size[1] - eps * (box_size[1] - y_prev)
        else:
            y_j = inv_mass(m_y)
        y_prev = y_j

        # loop over x
        for i in range(N[0]):

            index = j * N[0] + i

            x = i * box_size[0] / float(N[0])

            coords[index, 0] = x
            coords[index, 1] = y_j
            if y_j < fixed[0] or y_j > fixed[1]:
                uni_id += 1
                ids[index] = uni_id

    print(f"Set 'boundary_particle_max_id: {uni_id}' in rayleigh_taylor.yml")

    ids[ids == 0] = np.linspace(uni_id + 1, numPart, numPart - uni_id)

    if not np.array_equal(np.sort(ids), np.arange(1, numPart + 1)):
        raise ValueError("ids does not contain all unique integers from 1 to numPart")

    # density
    rho = density(coords[:, 1])

    # internal energy
    a = entropy(coords[:, 1])
    u = a * rho ** (gamma - 1) / (gamma - 1.0)

    # smoothing length
    h = smoothing_length(rho, m)

    # Velocity perturbation
    vel[:, 1] = vy(coords[:, 0], coords[:, 1])

    # File
    w = sw.Writer(
        unit_system=unyt.unit_systems.cgs_unit_system,
        boxsize=sw.cosmo_array(
            box_size, unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
        ),
        dimension=2,
    )

    w.gas.coordinates = coords * unyt.cm * w.a.comoving
    w.gas.velocities = vel * unyt.cm / unyt.s * w.a.comoving**0
    w.gas.masses = m * unyt.g * w.a.comoving**0
    w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
    w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
    w.gas.particle_ids = ids.astype(np.uint64)

    w.write(fileOutputName)

    # Add the extra fields not handled by the Writer
    with h5py.File(fileOutputName, "r+") as f:
        f["PartType0"].create_dataset("Density", data=rho, dtype="f")
