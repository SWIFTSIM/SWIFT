################################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
#                         2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
        "--lJ",
        action="store",
        dest="lJ",
        type=float,
        default=0.250,
        help="Jeans wavelength in box size unit",
    )

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
        default=50,
        help="Gas particle mass in solar mass",
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
        "--star_mass",
        action="store",
        type=float,
        default=29.7,
        help="Star particles mass in solar mass",
    )

    parser.add_argument(
        "--n_star",
        action="store",
        type=int,
        default=100,
        help="Number of star particles in the box",
    )

    parser.add_argument(
        "-o",
        action="store",
        dest="outputfilename",
        type=str,
        default="box.hdf5",
        help="output filename",
    )

    parser.add_argument('--random_positions', default=False, action="store_true",
                        help="Place the particles randomly in the box.")

    options = parser.parse_args()
    return options


########################################
# main
########################################

opt = parse_options()

# define standard units
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
N = (2 ** opt.level) ** 3  # number of particles

# Mean density
rho = opt.rho  # atom/cc
rho = rho * constants.m_p / units.cm ** 3

# Gas particle mass
m = opt.mass  # in solar mass
m = m * units.Msun

# Gas mass in the box
M = N * m

# Size of the box
L = (M / rho) ** (1 / 3.0)

# Jeans wavelength in box size unit
lJ = opt.lJ
lJ = lJ * L

# Gravitational constant
G = constants.G

# Jeans wave number
kJ = 2 * np.pi / lJ

# Velocity dispersion
sigma = np.sqrt(4 * np.pi * G * rho) / kJ

print("Boxsize                               : {}".format(L.to(units.kpc)))
print("Number of particles                   : {}".format(N))
print("Equivalent velocity dispertion        : {}".format(sigma.to(units.m / units.s)))

# Convert to code units
m = m.to(UnitMass).value
L = L.to(UnitLength).value
rho = rho.to(UnitMass / UnitLength ** 3).value
sigma = sigma.to(UnitVelocity).value

# Generate the particles

if opt.random_positions:
    print("Sampling random positions in the box")
    pos = np.random.random([N, 3]) * np.array([L, L, L])
else:
    points = np.linspace(0, L, 2**opt.level, endpoint=False)

    # Create a meshgrid of these points
    x, y, z = np.meshgrid(points, points, points)

    # Reshape the grids into a list of particle positions
    pos = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T

vel = np.zeros([N, 3])
mass = np.ones(N) * m
u = np.ones(N) * sigma ** 2
ids = np.arange(N)
h = np.ones(N) * 3 * L / N ** (1 / 3.0)
rho = np.ones(N) * rho

print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))


#####################
# Now, take care of the stars
#####################
N_star = opt.n_star
m_star = opt.star_mass * units.M_sun
m_star = m_star.to(UnitMass).value  # Convert the star mass to internal units

pos_star = np.random.random([N_star, 3]) * np.array([L, L, L])
vel_star = np.zeros([N_star, 3])

mass_star = np.ones(N_star) * m_star
h_star =  np.ones(N_star) * 3 * L / (N + N_star) ** (1 / 3.0)
ids_star = np.arange(N, N + N_star)
star_particle_type =  np.zeros(N_star) # Single star

# Make star explode at regular time, ~ per Myr
star_birth_time = np.linspace(0, 100, 100) * units.Myr
star_birth_time = star_birth_time.to(UnitTime).value

print("Number of stars       : {}".format(N_star))

#####################
# Finally write the ICs in the file
#####################

# File
fileOutput = h5py.File(opt.outputfilename, "w")
print("{} saved.".format(opt.outputfilename))

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, N_star, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, N_star, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3


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
grp.create_dataset("Densities", data=rho, dtype="f")

# Write star particle group
grp = fileOutput.create_group("/PartType4")
grp.create_dataset("Coordinates", data=pos_star, dtype="d")
grp.create_dataset("Velocities", data=vel_star, dtype="f")
grp.create_dataset("Masses", data=mass_star, dtype="f")
grp.create_dataset("ParticleIDs", data=ids_star, dtype="L")
grp.create_dataset("SmoothingLength", data=h_star, dtype="f")
grp.create_dataset("BirthMass", data=mass_star, dtype="f")
grp.create_dataset("BirthTime", data=star_birth_time, dtype="f")
grp.create_dataset("StellarParticleType", data=star_particle_type, dtype="i")
fileOutput.close()
