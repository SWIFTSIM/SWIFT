#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 10:39:41 2024

@author: darwinr
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
Plot the Fe mass time evolution. This script ensure metal mass conservation in idealized tests.
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
                        help="Plot in log.")

    parser.parse_args()
    args = parser.parse_args()
    files = args.files

    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError("You need to provide one file")

    return args, files

#%%

args, files = parse_option()
log = args.log

N_file = len(files)

time = np.zeros(N_file)
m_fe = np.zeros(N_file)

for i, filename in enumerate(tqdm(files)):
    # Load data
    data = sw.load(filename)
    boxsize = data.metadata.boxsize

    if hasattr(data.gas.metal_mass_fractions, 'fe'):
        m_fe_iter = data.gas.metal_mass_fractions.fe * data.gas.masses
    else:
        m_fe_iter = data.gas.metal_mass_fractions[:, 0] * data.gas.masses

    m_fe[i] = np.sum(m_fe_iter.value)
    time[i] = data.metadata.time.value


figsize = (4.6, 6.2)
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figsize, num=1)

if log:
    ax.set_yscale("log")

ax.set_ylabel("$M_{\mathrm{Fe}}$ [M$_\odot$]")
ax.set_xlabel("$t$ [Gyr]")

ax.plot(time, m_fe)
plt.savefig("metal_mass_evolution.png", format='png',
            bbox_inches='tight', dpi=300)
plt.close()
