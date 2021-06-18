#!/usr/bin/env python3

# ---------------------------------------------------------------------
# Write initial conditions for gas and stars where there are
# background particles on a uniform grid, and another layer of
# gas and stars that are sampled from a sine wave
# ---------------------------------------------------------------------

from swiftsimio import Writer
from swiftsimio.units import cosmo_units

import unyt
import numpy as np
from scipy import special as sp
from matplotlib import pyplot as plt

np.random.seed(666)

# Box is 1 Mpc
boxsize = 1000 * unyt.kpc


n_p = 20  # number of background particles in each dimension
n_s = 20  # number of background star particles in each dimension
n_sample_p = 10000  # how many hydro particles to sample the distribution
n_sample_s = 5000  # how many hydro particles to sample the distribution


# Get uniform distribution

xp = []
xs = []

dx = boxsize / n_p
ds = boxsize / n_s

for i in range(n_p):
    x = (i + 0.501) * dx
    for j in range(n_p):
        y = (j + 0.501) * dx
        for k in range(n_p):
            z = (k + 0.501) * dx
            xp.append(np.array([x, y, z], dtype=np.float))
xp = np.array(xp)
velp = np.zeros((n_p ** 3, 3), dtype=np.float)

for i in range(n_s):
    x = (i + 0.001) * ds
    for j in range(n_s):
        y = (j + 0.001) * ds
        for k in range(n_s):
            z = (k + 0.001) * ds
            xs.append(np.array([x, y, z], dtype=np.float))
xs = np.array(xs)
vels = np.zeros((n_s ** 3, 3), dtype=np.float)


amplitude = 0.5
velp_max = 20.0
vels_max = 200.0


def sine(x, amplitude=1.0):
    # raise the sine by 1.01 so you don't get negative probablities
    return amplitude * np.sin(x / boxsize.value * 2.0 * np.pi) + 1.01


def sample(n):
    samples = 0
    keep = np.zeros((n, 3), dtype=np.float)
    while samples < n:
        sample = np.zeros(3, dtype=np.float)

        found = False
        while not found:
            pick = np.random.uniform(0.0, boxsize.value, 1)
            prob_x = sine(pick, amplitude)
            confirm = np.random.uniform(0.0, 1.0) * (amplitude + 1.01)
            if confirm <= prob_x:
                sample[0] = pick  # * np.random.choice([1, -1])
                found = True
        found = False
        while not found:
            pick = np.random.uniform(0.0, boxsize.value, 1)
            prob_y = sine(pick, amplitude)
            confirm = np.random.uniform(0.0, 1.0) * (amplitude + 1.01)
            if confirm <= prob_y:
                sample[1] = pick
                found = True
        found = False
        while not found:
            pick = np.random.uniform(0.0, boxsize.value, 1)
            prob_z = sine(pick, amplitude)
            confirm = np.random.uniform(0.0, 1.0) * (amplitude + 1.01)
            if confirm <= prob_z:
                sample[2] = pick
                found = True

        keep[samples] = sample
        samples += 1

        if samples % 1000 == 0:
            print("sampled", samples, "/", n)

    return keep


xp_sampled = sample(n_sample_p)
velp_sampled = np.random.uniform(-1, 1, xp_sampled.shape) * velp_max

xs_sampled = sample(n_sample_s)
vels_sampled = np.random.uniform(-1, 1, xs_sampled.shape) * vels_max

vels_norm = np.sqrt(np.sum(vels_sampled ** 2, axis=1))
velp_norm = np.sqrt(np.sum(velp_sampled ** 2, axis=1))
#  print("min/max vels:", velp_norm.min(), velp_norm.max())
#  print("min/max vels:", vels_norm.min(), vels_norm.max())


xp_tot = np.vstack((xp, xp_sampled))
xp = unyt.unyt_array(xp_tot, boxsize.units)

xs_tot = np.vstack((xs, xs_sampled))
xs = unyt.unyt_array(xs_tot, boxsize.units)

vp_tot = np.vstack((velp, velp_sampled))
vp = unyt.unyt_array(vp_tot, unyt.km / unyt.s)

vs_tot = np.vstack((vels, vels_sampled))
vs = unyt.unyt_array(vs_tot, unyt.km / unyt.s)


w = Writer(cosmo_units, boxsize)
w.gas.coordinates = xp
w.stars.coordinates = xs
w.gas.velocities = vp
w.stars.velocities = vs
w.gas.masses = np.ones(xp.shape[0], dtype=np.float) * 1e6 * unyt.msun
w.stars.masses = np.random.uniform(1e8, 1e10, size=xs.shape[0]) * unyt.msun
# Generate internal energy corresponding to 10^4 K
w.gas.internal_energy = (
    np.ones(xp.shape[0], dtype=float) * (1e4 * unyt.kb * unyt.K) / (1e6 * unyt.msun)
)


# Generate initial guess for smoothing lengths based on MIPS
w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)
w.stars.generate_smoothing_lengths(boxsize=boxsize, dimension=3)

# If IDs are not present, this automatically generates
w.write("randomized-sine.hdf5")
