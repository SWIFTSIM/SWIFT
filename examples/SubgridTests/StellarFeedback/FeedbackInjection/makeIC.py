################################################################################
# This file is part of SWIFT.
# Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
    parser = argparse.ArgumentParser(
        description="Generate a homogeneous box with a Jeans unstable configuration and an exponentially decreasing density profile in the z-axis."
    )

    # Arguments
    parser.add_argument(
        "--rho",
        type=float,
        default=0.1,
        help="Mean gas density in atom/cm^3 (default: 0.1)",
    )

    parser.add_argument(
        "--mass",
        type=float,
        default=50.0,
        help="Gas particle mass in solar masses (default: 50.0)",
    )

    parser.add_argument(
        "--level",
        type=int,
        default=5,
        help="Resolution level: N = (2**level)^3 (default: 5)",
    )

    parser.add_argument(
        "--z_scale",
        type=float,
        default=0.01,
        help="Exponential scale height for density variation in the z-axis (default: 0.1)",
    )

    parser.add_argument(
        "--star_mass",
        type=float,
        default=29.7,
        help="Mass of the star in solar masses (default: 29.7)",
    )

    parser.add_argument(
        "--star_pos",
        type=store_as_array,
        default=None,
        help="Position of the star in internal units as a comma-separated list of three floats (e.g., '0.5,0.5,0.5'). Defaults to the center of the box.",
    )

    parser.add_argument(
        "--distribution",
        type=str,
        choices=["gaussian", "exponential"],
        default="gaussian",
        help="Choose the z-axis density profile: 'gaussian' or 'exponential' (default: exponential).",
    )

    parser.add_argument("--seed", type=int, default=1, help="Random seed")

    parser.add_argument(
        "-o",
        "--outputfilename",
        type=str,
        default="box.dat",
        help="Output filename (default: 'box.dat')",
    )

    # Parse the arguments
    args = parser.parse_args()
    return args


# Define the z-axis distribution
def generate_z_coordinates(N, L, z_scale, distribution="exponential"):
    z_center = L / 2  # Center of the box
    h_mean = L / N ** (1 / 3.0)  # Mean inter-particle separation in the midplane

    if distribution == "exponential":
        # Exponential distribution
        random_signs = np.random.choice([-1, 1], N)  # Random direction for Z offset
        z_offsets = -z_scale * np.log(np.random.uniform(0, 1, N))
        z_coords = (z_center + random_signs * z_offsets) % L
    elif distribution == "gaussian":
        # Gaussian distribution
        z_offsets = np.random.normal(
            loc=0.0, scale=h_mean, size=N
        )  # Gaussian distribution
        z_coords = (z_center + z_offsets) % L
    else:
        raise ValueError(f"Unknown distribution type: {distribution}")

    return z_coords


########################################
# main
########################################

opt = parse_options()
seed = opt.seed


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

np.random.seed(seed)

# Number of particles
N = (2 ** opt.level) ** 3  # number of particles

# Mean density
rho_mean = opt.rho  # atom/cc
rho_mean = rho_mean * constants.m_p / units.cm ** 3

# Gas particle mass
m = opt.mass  # in solar mass
m = m * units.Msun

# Gas mass in the box
M = N * m

# Size of the box
L = (M / rho_mean) ** (1 / 3.0)

# Gravitational constant
G = constants.G

print("Number of particles                   : {}".format(N))

# Convert to code units
m = m.to(UnitMass).value
L = L.to(UnitLength).value
rho_mean = rho_mean.to(UnitMass / UnitLength ** 3).value

# Generate the particles
pos = np.zeros([N, 3])

# Uniform distribution in X and Y
pos[:, 0] = np.random.random(N) * L
pos[:, 1] = np.random.random(N) * L

# Generate Z-coordinates based on the chosen distribution
pos[:, 2] = generate_z_coordinates(N, L, opt.z_scale * L, opt.distribution)

# Other particle properties
vel = np.zeros([N, 3])
mass = np.ones(N) * m
u = np.ones(N)
ids = np.arange(N)
h = np.ones(N) * 3 * L / N ** (1 / 3.0)
rho = np.ones(N) * rho_mean

print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))

#####################
# Add Background particles. These avoid Swift complaining that there are not
# enough particles in some top cells
#####################

# Number of background particles
level_background = opt.level - 1
N_background = (2 ** level_background) ** 3
print("Number of background particles         : {}".format(N_background))

# Generate background Cartesian grid
grid_spacing = L / (2 ** level_background)
grid_coords = np.linspace(grid_spacing / 2, L - grid_spacing / 2, 2 ** level_background)
x, y, z = np.meshgrid(grid_coords, grid_coords, grid_coords)
pos_background = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T

# Now randomly shift the whole grid so that it's not at the same location for
# different seeds.
shift_x = np.random.uniform(0.0, grid_spacing)
shift_y = np.random.uniform(0.0, grid_spacing)
shift_z = np.random.uniform(0.0, grid_spacing)

pos_background[:, 0] = (pos_background[:, 0] + shift_x) % L
pos_background[:, 1] = (pos_background[:, 1] + shift_y) % L
pos_background[:, 2] = (pos_background[:, 2] + shift_z) % L

# Assign properties to background particles
vel_background = np.zeros_like(pos_background)  # Zero velocity
mass_background = np.ones(N_background) * m  # Same mass
u_background = np.ones(N_background)  # Same internal energy
ids_background = np.arange(N, N + N_background)
h_background = np.ones(N_background) * (3 * L / (N + N_background) ** (1 / 3.0))
rho_background = np.ones(N_background) * (rho_mean / N * N_background)  # Adjust
# density for background

# Merge the properties
pos = np.vstack([pos, pos_background])
vel = np.vstack([vel, vel_background])
mass = np.concatenate([mass, mass_background])
u = np.concatenate([u, u_background])
ids = np.concatenate([ids, ids_background])
h = np.concatenate([h, h_background])
rho = np.concatenate([rho, rho_background])

N_gas_tot = N_background + N

#####################
# Now, take care of the star
#####################
N_star = 1
M_star = opt.star_mass * units.M_sun
pos_star = opt.star_pos

# Convert the star mass to internal units
M_star = [M_star.to(UnitMass).value]
print("Mass of the star (internal units)     : {:e}".format(M_star[0]))

# If no position was given, place the star at the center of the box
if pos_star is None:
    pos_star = np.ones([N_star, 3]) * L / 2

# Remaining required data
vel_star = np.zeros([N_star, 3])
h_star = np.ones(N_star) * 3 * L / N ** (1 / 3.0)  # Same as the gas particles
star_particle_type = np.ones(N_star) * 0  # Single star
star_id = [N + N_star]
star_birth_time = np.zeros(N_star)

print("Smoothing length of the star (internal units)     : {:e}".format(h_star[0]))


#####################
# Finally write the ICs in the file
#####################

# File
fileOutput = h5py.File(opt.outputfilename, "w")
print("{} saved.".format(opt.outputfilename))

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]
grp.attrs["NumPart_Total"] = [N_gas_tot, 0, 0, 0, N_star, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N_gas_tot, 0, 0, 0, N_star, 0]
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

# Particle group
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
