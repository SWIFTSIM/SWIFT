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

    # Do we want to generate positions randomly instead of a grid.
    parser.add_argument('--random_positions', default=False, action="store_true",
                        help="Place the particles randomly in the box. If false, the particles are placed on a regular carthesian grid")

    # Do we want to add a sheering effect?
    parser.add_argument('--add_shear', default=False, action="store_true",
                        help="Particle velocity is v_y = x/2 - boxsize/2")

    # Add a velocity boost
    parser.add_argument("--velocity",
                        action=store_as_array,
                        nargs=3,
                        type=float,
                        default=np.zeros(0),
                        help="Velocity boost to all particles")

    parser.add_argument(
        "--Z_R",
        action="store",
        type=float,
        default=1e-4,
        help="Metal mass fraction of the gas particles at x >= L/2",
    )

    parser.add_argument(
        "--Z_L",
        action="store",
        type=float,
        default=0,
        help="Metal mass fraction of the gas particles at x < L/2",
    )

    parser.add_argument(
        "--rho_R",
        action="store",
        type=float,
        default=1.0,
        help="Gas density if x >= L/2. In atom/cm3",
    )

    parser.add_argument(
        "--rho_L",
        action="store",
        type=float,
        default=1.0,
        help="Gas density if x < L/2. In atom/cm3",
    )

    parser.add_argument(
        "--L",
        action="store",
        type=float,
        default=1.0,
        help="Boxsize in kpc",
    )

    parser.add_argument(
        "--N_metal",
        action="store",
        type=float,
        default=10,
        help="Default number of metals in GEAR Metal mass fraction of the gas particle placed at the center of the box",
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

    options = parser.parse_args()
    return options

# %%
########################################
# main
########################################


opt = parse_options()

N_metal = opt.N_metal
random_positions = opt.random_positions
add_shear = opt.add_shear
level = opt.level
velocity = opt.velocity
dimension = opt.dimension

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
if dimension == 1:
    N = 2**level
elif dimension == 2:
    N = (2**level)**2
elif dimension == 3:
    N = (2**level)**3
else:
    raise ValueError("Only dimensions 1, 2, or 3 are supported.")

# Gas mass in the box
M = opt.mass * units.Msun  # in solar mass

# Gas particle mass
m = M / N

# Size of the box
L = opt.L*units.kpc

# Gravitational constant
G = constants.G

print("Dimension                             : {}".format(dimension))
print("Boxsize                               : {}".format(L.to(units.kpc)))
print("Number of particles                   : {}".format(N))

# Convert to code units
m = m.to(UnitMass).value
L = L.to(UnitLength).value

#####################
# Generate the particles
#####################

# %% Generate particles positions

# Calculate the inter-particle spacing (size of one grid cell)
if dimension == 1:
    grid_size = L / (2**level)
else:
    grid_size = L / (2**level) # This is the spacing along one dimension

# Define a small half-shift (s): this ensures particles are centered in their cell,
# moving them away from the 0 boundary.
shift = 0.5 * grid_size

if random_positions:
    print("Sampling random positions in the box (Range: [0, L))")
    # Add a minimal offset to ensure strict (0, L) if random() returns exactly 0
    epsilon = 1e-12 * L
    pos = np.random.random([N, dimension]) * (L - 2*epsilon) + epsilon

else:
    print("Generating carthesian grid in the box (Range: (0, L))")

    # Start the grid points at a small offset (half a grid cell)
    # This ensures the first particle is at 'shift' (e.g., L/2N) and the last is
    # at L - shift, moving them all strictly off the boundary.
    points = np.linspace(shift, L - shift, 2**level, endpoint=True)

    # Note on np.linspace:
    # We use endpoint=True here because we defined the range [shift, L-shift] for 2**level points.
    # The step size will now be exactly 'grid_size', centered correctly.

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

print("Inter-particle distance (code unit)   : {}".format(grid_size))

# %% Velocity

# Ensure 'x' exists and is the 0-th axis of pos
x = pos[:, 0]
xmid = L/2.0

# Prepare the velocity array with 0s
vel = np.tile(velocity, (N, 1))

# Add shearing velocity
if add_shear:
    print("Adding shear...")
    vy_shear = (x / 2.0 - L / 2.0) * units.kpc / units.Gyr
    vy_shear = vy_shear.to(UnitVelocity).value

    # Create the shear velocity array and apply it only to the y-axis (index 1)
    v_shear = np.zeros((N, dimension))
    if dimension > 1:
        v_shear[:, 1] = vy_shear
    else:
        raise ValueError("Shear along y-direction requires dimension > 1")

    # Sum contributions
    vel = vel + v_shear

# %% Add density jump
rho_L = opt.rho_L
rho_R = opt.rho_R
rho = np.zeros(N)  # Init rho

if rho_L != rho_R:
    I_L = np.argwhere(x < xmid).flatten()
    I_R = np.argwhere(x >= xmid).flatten()
    rho[I_L] = rho_L
    rho[I_R] = rho_R
else:
    rho = rho_L  # atom/cc

# Unit conversion (Do we need to change to dimension ?)
rho = rho * constants.m_p / units.cm ** dimension
rho = rho.to(UnitMass / UnitLength ** dimension).value  # Code units

print("Density of the particles (code unit)   : {}".format(rho))

# %% Init the rest of the variables
mass = np.ones(N) * m
u = np.zeros(N)
ids = np.arange(N)
h = np.ones(N) * dimension * L / N ** (1 / dimension)

# %% Add metallicity
Z_L = opt.Z_L
Z_R = opt.Z_R
metal_mass_fraction = np.zeros([N, N_metal])  # Init Fe metal mass fraction

# Get the particles with x >= xmid
I_R = np.argwhere(x >= xmid).flatten()
I_L = np.argwhere(x < xmid).flatten()

metal_mass_fraction[I_L, 0] = Z_L
metal_mass_fraction[I_R, 0] = Z_R

print(f"Particles with x>= {xmid} are given a Z_Fe = {
      Z_R} and particle with z < {xmid} are given Z_Fe = {Z_L}")

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
