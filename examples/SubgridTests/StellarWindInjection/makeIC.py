################################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
#               2025 Darwin Roduit (darwin.roduit@epfl.ch)
#               2025 Zachary Rey (zachary.rey@epfl.ch)
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
import unyt
import swiftsimio as sw


class store_as_array(argparse._StoreAction):
    """Provides numpy array as argparse arguments."""

    def __call__(self, parser, namespace, values, option_string=None):
        values = np.array(values)
        return super().__call__(parser, namespace, values, option_string)


def parse_options():

    usage = "usage: %prog [options] file"
    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument(
        "--T",
        action="store",
        dest="T",
        type=float,
        default=0.250,
        help="Target gas temperature in Kelvin",
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
        "--star_mass",
        action="store",
        type=float,
        default=29.7,
        help="Mass of the star in M_sun",
    )

    parser.add_argument(
        "--star_pos",
        action=store_as_array,
        nargs=3,
        type=float,
        default=None,
        help="Position of the star in internal units (works only for one star). By default, places the star at the center of the box",
    )

    parser.add_argument(
        "--star_type",
        action="store",
        type=int,
        default=0,
        help="The type of the star particle : 0=single star; 1=continuous IMF; 2=old single population star",
    )

    parser.add_argument(
        "--random_positions",
        default=False,
        action="store_true",
        help="Place the particles randomly in the box.",
    )

    parser.add_argument(
        "--num_stars",
        action="store",
        type=int,
        default=1,
        help="The number of stars in the simulation",
    )

    parser.add_argument(
        "-o",
        action="store",
        dest="outputfilename",
        type=str,
        default="box.hdf5",
        help="output filename",
    )

    # Ajouter mass etoile, position. Dans le code, dire que c'est une etoile discrete
    options = parser.parse_args()
    return options


def calculate_mu(T_k, H_frac):
    """
    Calculate mean molecular weight based on a temperature threshold.
    """
    T_transition = 1.1e4  # Kelvin

    if T_k > T_transition:
        # Fully ionized primordial gas
        return 4.0 / (8.0 - 5.0 * (1.0 - H_frac))
    else:
        # Neutral primordial gas
        return 4.0 / (1.0 + 3.0 * H_frac)


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

UnitMass = UnitMass_in_cgs * unyt.g
UnitLength = UnitLength_in_cgs * unyt.cm
UnitTime = UnitTime_in_cgs * unyt.s
UnitVelocity = UnitVelocity_in_cgs * unyt.cm / unyt.s

np.random.seed(1)

# Number of particles
N = (2**opt.level) ** 3  # number of gas particles

# Mean density
rho = opt.rho  # atom/cc
rho = rho * unyt.mp / unyt.cm**3

# Gas particle mass
m = opt.mass  # in solar mass
m = m * unyt.Msun

# Gas mass in the box
M = N * m

# Size of the box
L = (M / rho) ** (1 / 3.0)

# Thermodynamics properties
gamma = 5.0 / 3.0  # Ideal gas
H_frac = 0.76  # Standard primordial hydrogen mass fraction
mu = calculate_mu(opt.T, H_frac)
T = opt.T * unyt.K
u_unyt = (unyt.kb * T) / ((gamma - 1.0) * mu * unyt.mp)

print("Boxsize                               : {}".format(L.to(unyt.kpc)))
print("Number of particles                   : {}".format(N))
print("Target Temperature                    : {}".format(T))
print("Mean molecular weight (mu)            : {}".format(mu))
print("Specific Internal Energy (cgs)        : {}".format(u_unyt.to("erg/g")))

# Convert to code units
m = m.to(UnitMass).value
L = L.to(UnitLength).value
rho = rho.to(UnitMass / UnitLength**3).value
u_code = u_unyt.to(UnitVelocity**2).value

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
u = np.ones(N) * u_code
ids = np.arange(1, N + 1)
h = np.ones(N) * 3 * L / N ** (1 / 3.0)
rho = np.ones(N) * rho

print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))

############################
# Now, take care of the star
############################
N_star = opt.num_stars
M_star = opt.star_mass * unyt.Msun
pos_star = opt.star_pos

# Convert the star mass to internal units
M_star = [M_star.to(UnitMass).value]
print("Mass of the star (internal units)     : {:e}".format(M_star[0]))

# If no position was given, place the star at the center of the box
if pos_star is None and N_star != 0:
    pos_star = np.ones([N_star, 3]) * L / (1 + N_star)

# Remaining required data
vel_star = np.zeros([N_star, 3])
h_star = np.ones(N_star) * 3 * L / N ** (1 / 3.0)  # Same as the gas particles
star_particle_type = np.ones(N_star) * opt.star_type  # Single star by default
star_id = [N + N_star]
star_birth_time = np.zeros(N_star)  # birth time of the star in Myr

#####################
# Finally write the ICs in the file
#####################

# File
units = unyt.UnitSystem(
    "stellar_wind", unyt.Unit(UnitLength), unyt.Unit(UnitMass), unyt.Unit(UnitTime)
)

w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [L, L, L], UnitLength, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

# Gas particle group
w.gas.coordinates = pos * UnitLength * w.a.comoving
w.gas.velocities = vel * UnitVelocity * w.a.comoving**0
w.gas.masses = mass * UnitMass * w.a.comoving**0
w.gas.smoothing_lengths = h * UnitLength * w.a.comoving
w.gas.internal_energy = u * UnitVelocity**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

if N_star != 0:
    # Star particle group
    w.stars.coordinates = np.atleast_2d(pos_star) * UnitLength * w.a.comoving
    w.stars.velocities = vel_star * UnitVelocity * w.a.comoving**0
    w.stars.masses = np.array(M_star) * UnitMass * w.a.comoving**0
    w.stars.smoothing_lengths = h_star * UnitLength * w.a.comoving
    w.stars.particle_ids = np.array(star_id, dtype=np.uint64)

w.write(opt.outputfilename)
print("{} saved.".format(opt.outputfilename))

# Add the extra fields
with h5py.File(opt.outputfilename, "r+") as f:
    f["PartType0"].create_dataset("Densities", data=rho, dtype="f")

    if N_star != 0:
        grp = f["PartType4"]
        grp.create_dataset("BirthMass", data=M_star, dtype="f")
        grp.create_dataset("BirthTime", data=star_birth_time, dtype="f")
        grp.create_dataset("StellarParticleType", data=star_particle_type, dtype="i")
