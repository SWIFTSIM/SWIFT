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
import unyt
import swiftsimio as sw
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
        "--sink_mass",
        action="store",
        type=float,
        default=50,
        help="Sink particles mass in solar mass",
    )

    parser.add_argument(
        "--sinks_vel",
        action=store_as_array,
        nargs=3,
        type=float,
        default=np.array([10, 0, 0]),
        help="Sink particle velocity. All sinks get the same velocity",
    )

    parser.add_argument(
        "--sink_pos",
        action=store_as_array,
        nargs=3,
        type=float,
        default=np.array([0, 0, 0]),
        help="Sink particle position. Only use it to place one sink particle.",
    )

    parser.add_argument(
        "--n_sink",
        action="store",
        type=int,
        default=10,
        help="Number of sink particles in the box",
    )

    parser.add_argument(
        "-o",
        action="store",
        dest="outputfilename",
        type=str,
        default="box.hdf5",
        help="output filename",
    )

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
UnitTime_in_cgs = UnitLength_in_cgs / UnitVelocity_in_cgs

UnitMass = UnitMass_in_cgs * units.g
UnitLength = UnitLength_in_cgs * units.cm
UnitTime = UnitTime_in_cgs * units.s
UnitVelocity = UnitVelocity_in_cgs * units.cm / units.s

np.random.seed(1)

# Number of particles
N = (2**opt.level) ** 3  # number of particles

# Mean density
rho = opt.rho  # atom/cc
rho = rho * constants.m_p / units.cm**3

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
rho = rho.to(UnitMass / UnitLength**3).value
sigma = sigma.to(UnitVelocity).value

# Generate the particles
pos = np.random.random([N, 3]) * np.array([L, L, L])
vel = np.zeros([N, 3])
mass = np.ones(N) * m
u = np.ones(N) * sigma**2
ids = np.arange(1, N + 1)
h = np.ones(N) * 3 * L / N ** (1 / 3.0)
rho = np.ones(N) * rho

print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))


#####################
# Now, take care of the sink
#####################
N_sink = opt.n_sink
m_sink = opt.sink_mass * units.M_sun
m_sink = m_sink.to(UnitMass).value  # Convert the sink mass to internal units

if N_sink == 1:
    pos_sink = np.reshape(opt.sinks_vel, (N_sink, 3))
else:
    pos_sink = np.random.random([N_sink, 3]) * np.array([L, L, L])

if opt.sinks_vel is not None:
    vel_sink = np.tile(
        opt.sinks_vel, (N_sink, 1)
    )  # Duplicate the velocity for all sinks
else:
    np.zeros([N_sink, 3])

mass_sink = np.ones(N_sink) * m_sink
h_sink = np.ones(N_sink) * 3 * L / (N + N_sink) ** (1 / 3.0)
ids_sink = np.arange(N + 1, N + N_sink + 1)

#####################
# Finally write the ICs in the file
#####################

# File
UL = unyt.Unit(UnitLength_in_cgs * unyt.cm)
UM = unyt.Unit(UnitMass_in_cgs * unyt.g)
Ut = unyt.Unit(UnitTime_in_cgs * unyt.s)
unit_system = unyt.UnitSystem("homogeneous_box_sink_particles", UL, UM, Ut)

w = sw.Writer(
    unit_system=unit_system,
    boxsize=sw.cosmo_array(
        [L, L, L], UL, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

# Write Gas particle group
w.gas.coordinates = pos * UL * w.a.comoving
w.gas.velocities = vel * UL / Ut * w.a.comoving**0
w.gas.masses = mass * UM * w.a.comoving**0
w.gas.smoothing_lengths = h * UL * w.a.comoving
w.gas.internal_energy = u * UL**2 / Ut**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

# Write sink particle group
w.sinks.coordinates = pos_sink * UL * w.a.comoving
w.sinks.velocities = vel_sink * UL / Ut * w.a.comoving**0
w.sinks.masses = mass_sink * UM * w.a.comoving**0
w.sinks.particle_ids = ids_sink.astype(np.uint64)

w.write(opt.outputfilename)

with h5py.File(opt.outputfilename, "r+") as fileOutput:
    fileOutput["PartType0"].create_dataset("Densities", data=rho, dtype="f")
    fileOutput["PartType3"].create_dataset("SmoothingLength", data=h_sink, dtype="f")

print("{} saved.".format(opt.outputfilename))
