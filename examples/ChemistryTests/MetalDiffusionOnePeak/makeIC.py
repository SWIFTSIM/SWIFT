################################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
#               2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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

import h5py
import numpy as np
import argparse
from astropy import units
from astropy import constants


class store_as_array(argparse._StoreAction):
    """Provides numpy array as argparse arguments."""

    def __call__(self, parser, namespace, values, option_string=None):
        values = np.array(values)
        return super().__call__(parser, namespace, values, option_string)


def parse_options():

    usage = "usage: %prog [options] file"
    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument(
        "--rho",
        action="store",
        dest="rho",
        type=float,
        default=1,
        help="Mean gas density in atom/cm3",
    )

    parser.add_argument(
        "--mass",
        action="store",
        dest="mass",
        type=float,
        default=1e7,
        help="Total gass mass in the box. In solar mass",
    )

    parser.add_argument(
        "--level",
        action="store",
        dest="level",
        type=int,
        default=5,
        help="Resolution level: N = (2**l)**3",
    )

    parser.add_argument(
        "--velocity",
        action=store_as_array,
        nargs=3,
        type=float,
        default=np.zeros(0),
        help="Velocity boost to all particles")

    parser.add_argument(
        "--metal_mass_fraction",
        action="store",
        type=float,
        default=1e-4,
        help="Metal mass fraction of the gas particle placed at the center of the box",
    )

    parser.add_argument(
        "--N_metal",
        action="store",
        type=float,
        default=10,
        help="Default number of metals in GEAR Metal mass fraction of the gas particle placed at the center of the box",
    )

    parser.add_argument(
        "--epsilon",
        action="store",
        type=float,
        default=0.02,
        help="Radius of the homogeneous sphere to seed with metal mass",
    )

    parser.add_argument(
        "--dimension",
        action="store",
        type=int,
        default=3,
        help="Dimensionality of the problem",
    )

    parser.add_argument(
        "-o",
        action="store",
        dest="outputfilename",
        type=str,
        default="box.hdf5",
        help="output filename",
    )

    parser.add_argument(
        '--random_positions',
        default=False,
        action="store_true",
        help="Place the particles randomly in the box.")

    options = parser.parse_args()
    return options


########################################
# main
########################################

opt = parse_options()

N_metal = opt.N_metal
epsilon = opt.epsilon
random_positions = opt.random_positions
level = opt.level
velocity = opt.velocity
dimension = opt.dimension

# Define standard units
UnitMass_in_cgs = 1.988409870698051e43  # 10^10 M_sun in grams
UnitLength_in_cgs = 3.0856775814913673e21  # kpc in centimeters
UnitVelocity_in_cgs = 1e5  # km/s in centimeters per second
UnitCurrent_in_cgs = 1  # Amperes
UnitTemp_in_cgs = 1  # Kelvin
UnitTime_in_cgs = UnitLength_in_cgs / UnitVelocity_in_cgs

UnitMass = UnitMass_in_cgs * units.g
UnitLength = UnitLength_in_cgs * units.cm
UnitTime = UnitTime_in_cgs * units.s
UnitVelocity = UnitVelocity_in_cgs * units.cm / units.s

np.random.seed(1)

# Number of particles
if dimension == 1:
    N = 2**level
elif dimension == 2:
    N = (2**level)**2
elif dimension == 3:
    N = (2**level)**3
else:
    raise ValueError("Only dimensions 1, 2, or 3 are supported.")

# Mean density
rho = opt.rho  # atom/cc
rho = rho * constants.m_p / units.cm ** 3

# Gas mass in the box
M = opt.mass * units.Msun  # in solar mass

# Gas particle mass
m = M / N

# Size of the box
L = (M / rho) ** (1.0 / 3.0)

# Gravitational constant
G = constants.G

print("Dimension                             : {}".format(dimension))
print("Boxsize                               : {}".format(L.to(units.kpc)))
print("Number of particles                   : {}".format(N))

# Convert to code units
m = m.to(UnitMass).value
L = L.to(UnitLength).value
rho = rho.to(UnitMass / UnitLength ** 3).value

#####################
# Generate the particles
#####################

# Generate particles positions
if random_positions:
    print("Sampling random positions in the box")
    pos = np.random.random([N, dimension]) * L
else:
    print("Generating carthesian grid in the box")
    points = np.linspace(0, L, 2**level, endpoint=False)

    if dimension == 1:
        x = points
        pos = np.stack([x, np.zeros_like(x), np.zeros_like(x)], axis=1)

    elif dimension == 2:
        x, y = np.meshgrid(points, points, indexing='ij')
        pos = np.stack(
            [x.ravel(), y.ravel(), np.zeros_like(x.ravel())], axis=1)

    elif dimension == 3:
        x, y, z = np.meshgrid(points, points, points, indexing='ij')
        pos = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)

    else:
        raise ValueError("Only dimensions 1, 2, or 3 are supported.")

# Other properties
vel = np.tile(velocity, (N, 1))
mass = np.ones(N) * m
u = np.zeros(N)
ids = np.arange(N)
h = np.ones(N) * dimension * L / N ** (1 / dimension)
rho = np.ones(N) * rho

print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))
print("Density of the particles (code unit)   : {}".format(rho))


#####################
# Now, take care of the gas particles with the given metallicity
#####################

metallicity = opt.metal_mass_fraction

# Get the particles in a sphere of radius epsilon centered at the box center
box_center = np.full(dimension, L / 2)
if dimension == 1:
    x = pos[:, 0]
    r = np.abs(x-box_center)
elif dimension == 2:
    print("not implemented yet")
else:
    r = np.linalg.norm(pos-box_center, axis=1)

I = np.argwhere(r < epsilon).flatten()

print("Number of particles with given initial metallicity: {}".format(len(I)))

# Give the first metal the desired metallicity
metal_mass_fraction = np.zeros([N, N_metal])
metal_mass_fraction[I, 0] = metallicity

#####################
# Finally write the ICs in the file
#####################

# File
fileOutput = h5py.File(opt.outputfilename, "w")
print("{} saved.".format(opt.outputfilename))

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L] * dimension
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = dimension


# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = UnitLength_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = UnitMass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = UnitTime_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = UnitCurrent_in_cgs
grp.attrs["Unit temperature in cgs (U_T)"] = UnitTemp_in_cgs


# Write Gas particle group
grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=vel, dtype="f")
grp.create_dataset("Masses", data=mass, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("Density", data=rho, dtype="f")
grp.create_dataset("MetalMassFraction", data=metal_mass_fraction, dtype="d")

fileOutput.close()
