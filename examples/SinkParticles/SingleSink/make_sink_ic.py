################################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
# Copyright (c) 2024 Jonathan Davies (j.j.davies@ljmu.ac.uk)
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
from optparse import OptionParser
from astropy import units
from astropy import constants


def parse_options():

    usage = "usage: %prog [options] file"
    parser = OptionParser(usage=usage)

    parser.add_option(
        "--rho",
        action="store",
        dest="rho",
        type="float",
        default=0.1,
        help="Mean gas density in atom/cm3"
    )

    parser.add_option(
        "--sigma",
        action="store",
        dest="sigma",
        type="float",
        default=10.,
        help="Velocity dispersion of the gas in km/s (sets internal energy)."
    )

    parser.add_option(
        "--gas-mass",
        action="store",
        dest="gas_mass",
        type="float",
        default=1e5,
        help="Gas particle mass in solar masses"
    )

    parser.add_option(
        "--sink-mass",
        action="store",
        dest="sink_mass",
        type="float",
        default=5e6,
        help="Sink particle mass in solar masses"
    )

    parser.add_option(
        "--sink-velocity",
        action="store",
        dest="sink_mach",
        type="float",
        default=0.,
        help="Sink velocity as a multiple of the sound speed (i.e. Mach number)"
    )

    parser.add_option(
        "--level",
        action="store",
        dest="level",
        type="int",
        default=6,
        help="Resolution level: N = (2**l)**3"
    )

    parser.add_option(
        "-o",
        action="store",
        dest="outputfilename",
        type="string",
        default="ics.hdf5",
        help="output filename",
    )

    (options, args) = parser.parse_args()

    files = args

    return files, options


########################################
# main
########################################

files, opt = parse_options()

# define standard units
UnitMass_in_cgs = 1.988409870698051e43  # 10^10 M_sun in grams
UnitLength_in_cgs = 3.0856775814913673e21 # kpc in centimeters
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
m = opt.gas_mass * units.Msun  # in solar mass

# Sink particle mass
sm = opt.sink_mass * units.Msun # in solar mass

# Gas mass in the box
M = N * m

# Size of the box
L = (M / rho) ** (1 / 3.0)

# Gravitational constant
G = constants.G

print("Number of particles                   : {}".format(N))

# Convert to code units
m = m.to(UnitMass).value
sm = sm.to(UnitMass).value
L = L.to(UnitLength).value
rho = rho.to(UnitMass / UnitLength ** 3).value
sigma = (opt.sigma * units.km/units.s).to(UnitVelocity).value

# Generate the particles
pos = np.random.random([N, 3]) * np.array([L, L, L])
vel = np.zeros([N, 3])
mass = np.ones(N) * m
u = np.ones(N) * sigma**2
ids = np.arange(N)
h = np.ones(N) * 3 * L / N ** (1 / 3.0)
rho = np.ones(N) * rho

print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))

# Generate the sink
# Always put it 10% of the way through the box to give it room to move (if it's going to)

NSINK = 1

sink_pos = np.ones([NSINK, 3])
sink_pos[:,0] = L/10.
sink_pos[:,1] = L/2.
sink_pos[:,2] = L/2.

sink_mass = np.array([sm,])
sink_ids = np.array([2*ids[-1]])
sink_h = np.array([3 * L / N ** (1 / 3.0),])

gas_cs = np.sqrt(sigma**2 * 5./3. * ((5./3.)-1))

sink_vel = np.zeros([NSINK, 3])
sink_vel[:,0] += gas_cs * opt.sink_mach

print(f"Sink velocity: {gas_cs * opt.sink_mach}")

if gas_cs * opt.sink_mach > 0:
    sink_time_in_box = L*0.9 / (gas_cs * opt.sink_mach)
    print(f"Sink will leave box (neglecting dynamical friction) after time: {sink_time_in_box}")


# File
fileOutput = h5py.File(opt.outputfilename, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]
grp.attrs["NumPart_Total"] = [N, 0, 0, NSINK, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, NSINK, 0, 0]
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

# Particle groups
grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=vel, dtype="f")
grp.create_dataset("Masses", data=mass, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("Densities", data=rho, dtype="f")

grp = fileOutput.create_group("/PartType3")
grp.create_dataset("Coordinates", data=sink_pos, dtype="d")
grp.create_dataset("Velocities", data=sink_vel, dtype="f")
grp.create_dataset("Masses", data=sink_mass, dtype="f")
grp.create_dataset("SmoothingLength", data=sink_h, dtype="f")
grp.create_dataset("ParticleIDs", data=sink_ids, dtype="L")

fileOutput.close()

print(f"{opt.outputfilename} saved.")