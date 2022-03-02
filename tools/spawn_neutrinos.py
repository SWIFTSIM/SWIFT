#!/usr/bin/env python3

# Utility to spawn placeholder neutrinos in an existing initial conditions file.
# The neutrinos will be distributed uniformly in space and with with zero mass
# and velocity.

import h5py
import numpy as np
import sys

# This script does not need to be changed for different particle numbers.
# Usage: ./spawn_neutrinos.py filename

# Constants
Mpc_cgs = 3.08567758e24
Default_N_nu_per100Mpc = 72  # 72^3 neutrinos for a (100 Mpc)^3 box
Default_nr_neutrinos_per_Mpc3 = (Default_N_nu_per100Mpc / 100.0) ** 3

# Read command line arguments
if len(sys.argv) <= 1 or sys.argv[1] == "--help" or sys.argv[1] == "-h":
    print("Usage: ./spawn_neutrinos.py filename (use -h to show this message)")
    exit(0)

# Open the hdf5 file
fname = sys.argv[1]
f = h5py.File(fname, "r+")
print("Operating on '" + fname + "'")
print("")

# Check the unit system
if "Units" in f.keys() and "Unit length in cgs (U_L)" in f["Units"].attrs.keys():
    Length_Unit = f["Units"].attrs["Unit length in cgs (U_L)"] / Mpc_cgs  # Mpc
else:
    Length_Unit = 1.0  # Mpc

# Extract the box dimensions and volume
L = f["Header"].attrs["BoxSize"] / Length_Unit  # Mpc
V = L ** 3 if np.isscalar(L) else np.product(L)  # Mpc^3
if not np.isscalar(L) and len(L) != 3:
    raise ValueError("Box dimensions are not cubic")

# Check that the file does not have any neutrinos
nparts = f["Header"].attrs["NumPart_Total"]
while len(nparts) < 6:
    nparts = np.append(nparts, 0)
if nparts[6] != 0 or "PartType6" in f.keys():
    raise IOError("This file already has neutrinos.")

# Compute the default number of neutrinos (round to nearest cubic number)
Default_N_nu = round((Default_nr_neutrinos_per_Mpc3 * V) ** (1.0 / 3.0))
Default_Nr_neutrinos = int(Default_N_nu ** 3)

print("The box dimensions are " + str(L) + " Mpc.")
print(
    "The default number of neutrinos is "
    + "%g" % Default_N_nu_per100Mpc
    + "^3 per (100 Mpc)^3."
)
print(
    "The default number of neutrinos is "
    + "%g" % Default_N_nu
    + "^3 = "
    + str(Default_Nr_neutrinos)
    + "."
)
print("")

# Request the number of neutrino particles to be spawned (with default option)
Nr_neutrinos = int(
    input(
        "Enter the number of neutrinos (default " + "%d" % Default_Nr_neutrinos + "): "
    )
    or "%d" % Default_Nr_neutrinos
)

nparts[6] = Nr_neutrinos

print("")
print("The number of particles per type will be:")
print("{:25s}: {:12d}".format("Gas", nparts[0]))
print("{:25s}: {:12d}".format("Dark Matter", nparts[1]))
print("{:25s}: {:12d}".format("Background Dark Matter", nparts[2]))
print("{:25s}: {:12d}".format("Sinks", nparts[3]))
print("{:25s}: {:12d}".format("Stars", nparts[4]))
print("{:25s}: {:12d}".format("Black Holes", nparts[5]))
print("{:25s}: {:12d}".format("Neutrinos", nparts[6]))
print("")

firstID = int(nparts[0:6].sum() + 1)
print("The first particle ID of the first neutrino will be: " + str(firstID))
print("")

confirm = input("Enter y to confirm: ")
if confirm != "y":
    print("Not confirmed. Done for now.")
    exit(0)

print("")
print("Generating particle data...")

# Generate particle data
x = np.random.uniform(0, L, (Nr_neutrinos, 3)) * Length_Unit
v = np.zeros((Nr_neutrinos, 3))
m = np.zeros(Nr_neutrinos)
ids = np.arange(firstID, firstID + Nr_neutrinos)

print("Writing particle data...")

# Store the particle data
f.create_group("PartType6")
f["PartType6/Coordinates"] = x
f["PartType6/Velocities"] = v
f["PartType6/Masses"] = m
f["PartType6/ParticleIDs"] = ids

# Update the header
f["Header"].attrs["NumPart_Total"] = nparts

print("All done.")

# And close
f.close()
