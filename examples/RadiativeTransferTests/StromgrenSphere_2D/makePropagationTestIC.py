#!/usr/bin/env python3

# ---------------------------------------------------------------------
# Add a single star in the center of a glass distribution
# Intended for the photon propagation test.
# ---------------------------------------------------------------------

from swiftsimio import Writer
import unyt
import numpy as np
import h5py


glass = h5py.File("glassPlane_128.hdf5", "r")
parts = glass["PartType0"]
xp = parts["Coordinates"][:]
h = parts["SmoothingLength"][:]
glass.close()

# find particles closest to the center
# and select a couple of them to put the star in their middle
r = np.sqrt(np.sum((0.5 - xp) ** 2, axis=1))
mininds = np.argsort(r)
center_parts = xp[mininds[:4]]
xs = center_parts.sum(axis=0) / center_parts.shape[0]

# Double-check all particles for boundaries
for i in range(2):
    mask = xp[:, i] < 0.0
    xp[mask, i] += 1.0
    mask = xp[:, i] > 1.0
    xp[mask, i] -= 1.0


unitL = unyt.cm
t_end = 1e-3 * unyt.s
edgelen = unyt.c.to("cm/s") * t_end * 2.0
edgelen = edgelen.to(unitL)
boxsize = unyt.unyt_array([edgelen.v, edgelen.v, 0.0], unitL)

xs = unyt.unyt_array(
    [np.array([xs[0] * edgelen, xs[1] * edgelen, 0.0 * edgelen])], unitL
)
xp *= edgelen
h *= edgelen


w = Writer(unit_system=unyt.unit_systems.cgs_unit_system, box_size=boxsize, dimension=2)

w.gas.coordinates = xp
w.stars.coordinates = xs
w.gas.velocities = np.zeros(xp.shape) * (unyt.cm / unyt.s)
w.stars.velocities = np.zeros(xs.shape) * (unyt.cm / unyt.s)
w.gas.masses = np.ones(xp.shape[0], dtype=np.float64) * 100 * unyt.g
w.stars.masses = np.ones(xs.shape[0], dtype=np.float64) * 100 * unyt.g
w.gas.internal_energy = (
    np.ones(xp.shape[0], dtype=np.float64) * (300.0 * unyt.kb * unyt.K) / (unyt.g)
)

w.gas.smoothing_length = h
w.stars.smoothing_length = w.gas.smoothing_length[:1]

w.write("propagationTest-2D.hdf5")
