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
    
    # Do we want to generate positions randomly instead of a grid.
    parser.add_argument('--random_position', default=False, action="store_true",
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
        default=10.0,
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
        "-o",
        action="store",
        dest="outputfilename",
        type=str,
        default="box.hdf5",
        help="output filename",
    )

    options = parser.parse_args()
    return options

#%%
########################################
# main
########################################

opt = parse_options()

N_metal = opt.N_metal
random_position = opt.random_position
add_shear = opt.add_shear
level = opt.level
velocity = opt.velocity

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
N = (2 ** level) ** 3  # number of particles

# Gas mass in the box
M = opt.mass * units.Msun  # in solar mass

# Gas particle mass
m = M / N

# Size of the box
L = opt.L*units.kpc

# Gravitational constant
G = constants.G

print("Boxsize                               : {}".format(L.to(units.kpc)))
print("Number of particles                   : {}".format(N))

# Convert to code units
m = m.to(UnitMass).value
L = L.to(UnitLength).value

#%% Generate particles positions

if random_position:
    print("Sampling random positions in the box")
    pos = np.random.random([N, 3]) * np.array([L, L, L])
else:
    points = np.linspace(0, L, 2**level, endpoint=False)

    # Create a meshgrid of these points
    x, y, z = np.meshgrid(points, points, points)

    # Reshape the grids into a list of particle positions
    pos = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T

xmid = L/2.0
x = pos[:, 0]
print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))

#%% Velocity
vel = np.tile(velocity, (N, 1)) # np.zeros([N, 3])

# Add shear
if add_shear:
    vy_shear = x/2.0 - L/2.0
    
    # Create the 3D shear velocity array [0, vy_shear, 0]
    v_shear = np.zeros_like(vel)  # Initialize with zeros
    v_shear[:, 1] = vy_shear      # Assign vy_shear to the y-component
    
    # Sum all contributions
    vel = vel + v_shear

#%% Add density jump
rho_L = opt.rho_L
rho_R = opt.rho_R
rho = np.zeros(N) # Init rho

if rho_L != rho_R:
    I_L = np.argwhere(x >= xmid).flatten()
    I_R = np.argwhere(x < xmid).flatten()
    rho[I_L] = rho_L
    rho[I_R] = rho_R
else:
    rho = rho_L  # atom/cc
    
# Unit conversion
rho = rho * constants.m_p / units.cm ** 3
rho = rho.to(UnitMass / UnitLength ** 3).value # Code units

print("Density of the particles (code unit)   : {}".format(rho))

#%% Init the rest of the variables
mass = np.ones(N) * m
u = np.zeros(N)
ids = np.arange(N)
h = np.ones(N) * 3 * L / N ** (1 / 3.0)

#%% Add metallicity
Z_L = opt.Z_L
Z_R = opt.Z_R
metal_mass_fraction = np.zeros([N, N_metal]) # Init Fe metal mass fraction

# Get the particles with x >= xmid
I_R = np.argwhere(x >= xmid).flatten()
I_L = np.argwhere(x < xmid).flatten()

print(I_L)

metal_mass_fraction[I_L, 0] = Z_L
metal_mass_fraction[I_R, 0] = Z_R

print(f"Particles with x>= {xmid} are given a Z_Fe = {Z_R} and particle with z < {xmid} are given Z_Fe = {Z_L}")

#####################
# Finally write the ICs in the file
#####################

# File
fileOutput = h5py.File(opt.outputfilename, "w")
print("{} saved.".format(opt.outputfilename))

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
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
grp.create_dataset("Density", data=rho, dtype="f")
grp.create_dataset("MetalMassFraction", data=metal_mass_fraction, dtype="d")

fileOutput.close()
