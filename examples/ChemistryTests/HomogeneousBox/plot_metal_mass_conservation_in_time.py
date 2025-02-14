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
m_elements = np.zeros((N_file, n_elements))
m_feedback = np.zeros((N_file, n_elements))
m_diff = np.zeros((N_file, n_elements))
m_diff_abs = np.zeros((N_file, n_elements))

# Feedback and diffused metal masses
data = sw.load(files[0])
has_feedback_masses = hasattr(data.gas, "feedback_metal_masses")
has_diffused_masses = hasattr(data.gas, "diffused_metal_masses")

# Loop through files
for i, filename in enumerate(tqdm(files)):
    # Load data
    data = sw.load(filename)
    time[i] = data.metadata.time.value

    # Metal mass fractions
    for j, element in enumerate(elements):
        m_iter = getattr(data.gas.metal_mass_fractions, element) * data.gas.masses
        m_elements[i, j] = np.sum(m_iter.value)


    if has_feedback_masses:
        for j in range(n_elements):
            m_feedback[i, j] = np.sum(data.gas.feedback_metal_masses[:, j].to(unyt.Msun).value)

    if has_diffused_masses:
        for j in range(n_elements):
            m_diff_abs[i, j] = np.sum(np.abs(data.gas.diffused_metal_masses[:, j]).to(unyt.Msun).value)
            m_diff[i, j] = np.sum(data.gas.diffused_metal_masses[:, j].to(unyt.Msun).value)

# %%
# Plot each element's evolution
figsize = (6, 4)

for j, element in enumerate(elements):
    
    if element != "fe":
        continue
    
    fig, ax = plt.subplots(figsize=figsize)
    if log or symlog:
        ax.set_yscale("log")

    ax.set_title(f"Evolution of {element}")
    ax.set_ylabel(f"$M_{{\\mathrm{{{element}}}}}$ [M$_\\odot$]")
    ax.set_xlabel("$t$ [Gyr]")

    ax.plot(time, m_elements[:, j], label=f"{element} mass")

    if has_feedback_masses:
        ax.plot(time, m_feedback[:, j], label=f"{element} injected (feedback)")
    # if has_diffused_masses:
        # ax.plot(time, m_diff_abs[:, j], label=f"Sum in |.| of {element} diffused")
    ax.legend()

    plt.tight_layout()
    plt.savefig(f"metal_mass_evolution_{element}.png", format="png", dpi=300)
    plt.close()

#%%
# Plot each element's evolution

if has_diffused_masses:
    figsize = (6, 4)

    for j, element in enumerate(elements):
        
        if element != "fe":
            continue
    
        fig, ax = plt.subplots(figsize=figsize)
        if log:
            ax.set_yscale("log")
            m_diff[:, j] = np.abs(m_diff[:, j])

        if symlog:
            ax.set_yscale("symlog", linthresh=linthresh)

        ax.set_title(f"Evolution of {element}")
        ax.set_ylabel(f"$M_{{\\mathrm{{{element}}}}}$ [M$_\\odot$]")
        ax.set_xlabel("$t$ [Gyr]")

        ax.plot(time, m_diff[:, j], label=f"{element} diffused")
        ax.legend()

        plt.tight_layout()
        plt.savefig(f"metal_mass_evolution_diffused_{element}.png", format="png", dpi=300)
        plt.close()
