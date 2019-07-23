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
from scipy.optimize import bisect

# Generates a swift IC file for the Rayleigh-Taylor instability in a periodic
# box following the conditions given in Saitoh and Makino 2013: 1202.4277v3

# Parameters
N = [128, 192]  # Particles along one edge
gamma = 7./5.  # Gas adiabatic index
dv = 0.025    # velocity perturbation
rho_h = 2    # high density region
rho_l = 1    # Low density region
g = -0.5    # gravitational acceleration
box_size = [1., 1.5]  # size of the box

fixed = [0.1, 1.4]  # y-range of non fixed particles
perturbation_box = [0.3, 1.2]  # y-range for the velocity perturbation
fileOutputName = "rayleigh_taylor.hdf5"


# ---------------------------------------------------

if (N[0] / box_size[0] != N[1] / box_size[1]):
    raise Exception("Suppose the same ratio for each directions")


def density(y):
    """
    Mass density as function of the position y.
    """
    if isinstance(y, float) or isinstance(y, int):
        y = np.array(y)

    ind = y < 0.5 * box_size[1]
    rho = np.zeros(y.shape)
    tmp = (gamma - 1.) * g / (gamma * P0)
    alpha = 1. / (gamma - 1.)

    rho[ind] = rho_l * (1 + rho_l * tmp * (y[ind] - 0.5 * box_size[1]))**alpha
    rho[~ind] = rho_h * (1 + rho_h * tmp * (y[~ind] - 0.5 * box_size[1]))**alpha

    return rho


def mass(y):
    """
    Integral of the density
    """
    if isinstance(y, float) or isinstance(y, int):
        y = np.array(y)

    ind = y < 0.5 * box_size[1]
    m = np.zeros(y.shape)

    B = (gamma - 1.) * g / (gamma * P0)
    alpha = 1. / (gamma - 1.)

    m[ind] = (1 + B * rho_l * (y[ind] - 0.5 * box_size[1]))**(alpha + 1)

    m[~ind] = (1 + B * rho_h * (y[~ind] - 0.5 * box_size[1]))**(alpha + 1)

    m -= (1 - 0.5 * B * box_size[1] * rho_l)**(alpha + 1)

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
    a[ind] = P0 * rho_l**(-gamma)
    a[~ind] = P0 * rho_h**(-gamma)
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
    ymin = 0.
    ymax = box_size[1]
    A = density(ymax) - density(ymin)
    A /= density(ymax) + density(ymin)
    return np.sqrt(A * np.abs(g) * ymax / (2. * np.pi))


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
    uni_id = 1

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
            if (y_j < fixed[0] or y_j > fixed[1]):
                ids[index] = uni_id
                uni_id += 1

    print("You need to compile the code with "
          "--enable-boundary-particles=%i" % uni_id)
    ids[ids == 0] = np.linspace(uni_id, numPart, numPart-uni_id+1)

    # density
    rho = density(coords[:, 1])

    # internal energy
    a = entropy(coords[:, 1])
    u = a * rho**(gamma-1) / (gamma - 1.)

    # smoothing length
    h = smoothing_length(rho, m)

    # Velocity perturbation
    vel[:, 1] = vy(coords[:, 0], coords[:, 1])

    # File
    fileOutput = h5py.File(fileOutputName, 'w')

    # Header
    grp = fileOutput.create_group("/Header")
    grp.attrs["BoxSize"] = box_size
    grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFileOutputsPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = 0
    grp.attrs["Dimension"] = 2

    # Units
    grp = fileOutput.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 1.
    grp.attrs["Unit mass in cgs (U_M)"] = 1.
    grp.attrs["Unit time in cgs (U_t)"] = 1.
    grp.attrs["Unit current in cgs (U_I)"] = 1.
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.

    # Particle group
    grp = fileOutput.create_group("/PartType0")
    ds = grp.create_dataset('Coordinates', (numPart, 3), 'd')
    ds[()] = coords
    ds = grp.create_dataset('Velocities', (numPart, 3), 'f')
    ds[()] = vel
    ds = grp.create_dataset('Masses', (numPart, 1), 'f')
    ds[()] = m.reshape((numPart, 1))
    ds = grp.create_dataset('SmoothingLength', (numPart, 1), 'f')
    ds[()] = h.reshape((numPart, 1))
    ds = grp.create_dataset('InternalEnergy', (numPart, 1), 'f')
    ds[()] = u.reshape((numPart, 1))
    ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
    ds[()] = ids.reshape((numPart, 1))
    ds = grp.create_dataset('Density', (numPart, 1), 'f')
    ds[()] = rho.reshape((numPart, 1))

    fileOutput.close()
