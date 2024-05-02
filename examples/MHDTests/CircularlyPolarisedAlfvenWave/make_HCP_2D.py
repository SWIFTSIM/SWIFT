#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt

Nside = 8

# Generate a 2D HCP lattice

# Generate a regular lattice
N = Nside * Nside
grid_x = np.linspace(1./(2 * Nside), 1 - 1./(2 * Nside), Nside)
grid_y = np.linspace(1./(2 * Nside), 1 - 1./(2 * Nside), Nside)
xv, yv = np.meshgrid(grid_x, grid_y)
pos = np.array((xv.flatten(), yv.flatten(), np.zeros(N))).T

# Shift every second row
mask = pos[:,1] % (2. / Nside) == (3. / (2 * Nside))
pos[mask, 0] += 1. / (2 * Nside)

# Shift everything to avoid particles on the very edge
pos[:, 0] -= 1. / (4 * Nside)

# Compute sensible values for h
h = np.ones(N) * (1. / (2 * Nside))
