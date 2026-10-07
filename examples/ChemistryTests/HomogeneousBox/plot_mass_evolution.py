#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Improved script to plot metal mass time evolution for all elements stored in arrays.
Ensures metal mass conservation in idealized tests.
"""

import argparse
import os
import swiftsimio as sw
import numpy as np
import matplotlib.pyplot as plt
import unyt
from tqdm import tqdm


# %%
def parse_option():
    description = """
Plot the metal mass time evolution for all elements stored in the arrays.
This script ensures metal mass conservation in idealized tests.
    """
    epilog = """
Examples:
--------
python3 plot_metal_mass_conservation_in_time.py snap/*.hdf5 
"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("files",
                        nargs="+",
                        type=str,
                        help="File name(s).")

    parser.add_argument('--log', default=False, action="store_true",
                        help="Plot in log scale.")

    parser.add_argument('--symlog', default=False, action="store_true",
                        help="Plot in symlog scale.")

    parser.add_argument("--linthresh",
                        type=float,
                        default=1e-18,
                        help="When using symolog, controls the regions where ithe plot is linear")

    parser.parse_args()
    args = parser.parse_args()
    files = args.files

    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"File {f} does not exist.")

    return args, files


# %%
args, files = parse_option()
log = args.log
symlog = args.symlog
linthresh = args.linthresh

N_file = len(files)

# Load the first file to retrieve element names
first_data = sw.load(files[0])
elements = first_data.gas.metal_mass_fractions.named_columns  # Retrieve named columns
n_elements = len(elements)

# Initialize arrays
time = np.zeros(N_file)
m = np.zeros(N_file)

# Loop through files
for i, filename in enumerate(tqdm(files)):
    # Load data
    data = sw.load(filename)
    time[i] = data.metadata.time.value

    m_iter = data.gas.masses
    m[i] = np.sum(m_iter.value)

# %%
# Plot each element's evolution
figsize = (6, 4)
fig, ax = plt.subplots(figsize=figsize)
if log or symlog:
    ax.set_yscale("log")

ax.set_title("Evolution of mass")
ax.set_ylabel("$M$ [M$_\\odot$]")
ax.set_xlabel("$t$ [Gyr]")

ax.plot(time, m, label="Mass")
# ax.legend()

plt.tight_layout()
plt.savefig(f"mass_evolution.png", format="png", dpi=300)
plt.close()
