################################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
#               2025 Darwin Roduit (darwin.roduit@epfl.ch)
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
        default=0.1,
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
        "--boxsize",
        action="store",
        dest="boxsize",
        type=float,
        default=None,
        help="Boxzise in kpc",
    )

    parser.add_argument("--star_mass",
                        action="store",
                        type=float,
                        default=29.7,
                        help="Mass of the star in M_sun")

    parser.add_argument("--star_pos",
                        action=store_as_array,
                        nargs=3,
                        type=float,
                        default=None,
                        help="Position of the star in internal units. By default, places the star at the center of the box")

    parser.add_argument("--star_type",
                        action="store",
                        type=str,
                        default="single_star",
                        choices=["single_star", "continuous_IMF","SSP"],
                        help="Type of the star for the GEAR model.")


    parser.add_argument(
        "-o",
        action="store",
        dest="outputfilename",
        type=str,
        default="box.hdf5",
        help="output filename",
    )


#Ajouter mass etoile, position. Dans le code, dire que c'est une etoile discrete
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

# Size of the box
if opt.boxsize is None:
    # If we don't give the boxsize, then we can determine it
    # Gas mass in the box
    M = N * m
    L = (M / rho) ** (1 / 3.0)
else:
    # If we fixed the boxsize, then the number of particles might change
    L = opt.boxsize*units.kpc
    M = rho * L**3

    print(M.to(units.Msun))
    N = int(np.ceil(M.to(units.Msun) / m))
    print(N)

# Gravitational constant
G = constants.G

print("Boxsize                               : {}".format(L.to(units.kpc)))
print("Total mas                             : {}".format(M.to(units.Msun)))
print("Number of particles                   : {}".format(N))

# Convert to code units
m = m.to(UnitMass).value
L = L.to(UnitLength).value
rho = rho.to(UnitMass / UnitLength ** 3).value

# Generate the particles
pos = np.random.random([N, 3]) * np.array([L, L, L])
vel = np.zeros([N, 3])
mass = np.ones(N) * m
u = np.ones(N)
ids = np.arange(N)
h = np.ones(N) * 3 * L / N ** (1 / 3.0)
rho = np.ones(N) * rho

print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))


#####################
# Now, take care of the star
#####################
N_star = 1
M_star = opt.star_mass*units.M_sun
pos_star = opt.star_pos

# Convert the star mass to internal units
M_star = [M_star.to(UnitMass).value]
print('Mass of the star (internal units)     : {:e}'.format(M_star[0]))

# If no position was given, place the star at the center of the box
if pos_star is None:
    pos_star = np.ones([N_star, 3]) * L/2

#Remaining required data
# vel_star = np.zeros([N_star, 3])
vel_star = np.zeros([N_star, 3])
h_star = np.ones(N_star) * 3 * L / N ** (1 / 3.0) # Same as the gas particles

if opt.star_type == "single_star":
    star_type = 0
elif opt.star_type == "continuous_IMF":
    star_type = 1
elif opt.star_type == "SSP":
    star_type = 2
else:
    raise ValueError(f"Star type {opt.star_type} is not known!")

star_particle_type = np.ones(N_star)*star_type
star_id = [N + N_star]
star_birth_time = np.zeros(N_star)

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
grp.create_dataset("Masses", data=M_star, dtype="f")
grp.create_dataset("ParticleIDs", data=star_id, dtype="L")
grp.create_dataset("SmoothingLength", data=h_star, dtype="f")
grp.create_dataset("BirthMass", data=M_star, dtype="f")
grp.create_dataset("BirthTime", data=star_birth_time, dtype="f")
grp.create_dataset("StellarParticleType", data=star_particle_type, dtype="i")

fileOutput.close()
