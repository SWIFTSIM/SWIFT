#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plots the evolution of the DM particles with the CSDS.
"""

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from velociraptor import load
from velociraptor.particles import load_groups

# Include the CSDS library to the path before importing the CSDS python module
sys.path.append("../../../csds/src/.libs/")

# Now we can load the csds
import libcsds as csds

#%%
def get_redshift(a):
    return 1 / a - 1


def get_scale_factor(z):
    return 1 / (z + 1)


#%% Argparse options
def parse_option():
    description = """"Plots the evolution of the DM particles with the CSDS. """
    epilog = """
Examples:
--------

python3 csds_analysis.py csds_index_0000.dump
python3 csds_analysis.py csds_index_0000.dump --halo 2
python3 csds_analysis.py csds_index_0000.dump --halo 2 -n 250
"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("csds_file", type=str, help="CSDS .dump file name.")

    parser.add_argument(
        "--halo",
        action="store",
        type=int,
        dest="halo_index",
        default=1,
        help="Halo id of which we want trace the evolution. The halo id are given by velociraptor.",
    )

    parser.add_argument(
        "-n",
        action="store",
        type=int,
        default=500,
        help="Number of timestamps to fetch the CSDS.",
    )

    parser.add_argument(
        "--stf_file",
        action="store",
        type=str,
        dest="stf_file",
        default="stf_output/snap_0031",
        help="stf file basename.",
    )

    args = parser.parse_args()

    return args


#%% Parse the arguments
args = parse_option()

stf_file = args.stf_file
halo_index = args.halo_index
filename = args.csds_file
n_files = 1  # Number of csds files
n = args.n  # Number of timestamps to fetch the CSDS
dt = 1e-3  # small timestep for csds start (starts at t0+dt, see below)

# These are the number of timestamps we want for z > 2 (n_high_z) and z <=2 (n_low_z). This allows to provide more images outputs for 0 < z < 0 for the movie.
n_high_z = int(n / 10)
n_low_z = int(9 / 10 * n)

print("\n-------------------------------------")
print("Welcome to csds_anlysis.py !\n")

#%% Extract particles from halo
# Load Velociraptor data
print("Extracting halo data...")

catalogue = load(stf_file + ".properties.0")
groups = load_groups(stf_file + ".catalog_groups.0", catalogue=catalogue)

# Get the particles in halo_index
particles, unbound_particles = groups.extract_halo(halo_index=halo_index)

# Get the IDs
IDs = particles.particle_ids
N_particles = len(IDs)

#%% Now, we can open the CSDS, get the evolution of the particles though time with their IDs
cm_positions = np.zeros((n_files, n, 3))
positions = np.empty((n, N_particles, 3))
velocities = np.empty((n, N_particles, 3))
n_particles = np.zeros((n_files, n), dtype=int)
m_tot = np.zeros((n_files, n))

file_index = 0

if filename.endswith(".dump"):
    filename = filename[:-5]

print("Openening the CSDS...")

#%%
with csds.Reader(
    filename, verbose=0, number_index=10, restart_init=False, use_cache=True
) as reader:

    print("The CSDS is opened.")

    # Check the time limits (scale-factors in this example)
    t0, t1 = reader.get_time_limits()

    print("Scale factor limits: [{:e}, {:e}]".format(t0, t1))
    print(np.log10(t0), np.log10(t1))
    print("Redshift limits: [{:e}, {:e}]".format(get_redshift(t0), get_redshift(t1)))

    # Ensure that the fields are present
    fields = ["Coordinates", "Masses", "ParticleIDs", "Velocities"]
    missing = set(fields).difference(
        set(reader.get_list_fields(part_type=csds.particle_types.gas))
    )

    if missing:
        raise Exception("Fields %s not found in the logfile." % missing)

    # Create the list of IDs for *all* particle types
    IDs_list = [None] * csds.particle_types.count
    IDs_list[csds.particle_types.dark_matter] = IDs

    # Read the particles by IDs
    out = reader.get_data(fields=fields, time=t0, filter_by_ids=IDs_list)

    # Print the missing ids
    dm_ids, ids_found = set(IDs), set(out["ParticleIDs"])
    diff_ids = list(dm_ids.difference(ids_found))
    diff_found = list(ids_found.difference(dm_ids))
    print("The following ids were not found: ", np.array(diff_ids))
    print("The following ids are wrongly missing: ", np.array(diff_found))

    # Reverse the time: it is better to start from the end and rewind to the beginning
    high_z = np.linspace(t0 + dt, get_scale_factor(2), num=n_high_z, endpoint=False)
    low_z = np.linspace(get_scale_factor(2), t1, num=n_low_z)
    times = np.concatenate((high_z, low_z))[::-1]

    print("Starting to move through time")
    for i, t in enumerate(times):
        out = reader.get_data(fields=fields, time=t, filter_by_ids=IDs_list)

        # Get positions and velocities
        positions[i] = out["Coordinates"]
        velocities[i] = out["Velocities"]

        # Compute CM of halo
        n_particles[file_index, i] = len(out["ParticleIDs"])
        m_tot[file_index, i] = np.sum(out["Masses"])
        pos_cm = np.sum(out["Masses"][:, np.newaxis] * out["Coordinates"], axis=0)
        pos_cm /= m_tot[file_index, i]
        cm_positions[file_index, i, :] = pos_cm

print("The CSDS is now closed.")
#%% Do a movie showing the evolution of the halo (positions)
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots(num=1, figsize=(5, 5), dpi=300)
fig.tight_layout()


def animate(index):
    # Read the positions in the right order
    i = -index - 1
    x = positions[i, :, 0]
    y = positions[i, :, 1]
    scale_factor = times[i]
    redshift = 1 / scale_factor - 1

    # Clear the previous data
    ax.cla()
    ax.text(
        0.975,
        0.975,
        "z = {:.2f}".format(redshift),
        ha="right",
        va="top",
        transform=ax.transAxes,
        color="k",
    )
    ax.scatter(x, y, c="k", zorder=1, marker=".")
    ax.set_aspect("equal", "box")
    return ax


print("Making movie...")
animation = FuncAnimation(
    fig, animate, range(positions.shape[0]), fargs=[], interval=n / 24
)
animation.save("halo_evolution.mp4")
print("End :)")
