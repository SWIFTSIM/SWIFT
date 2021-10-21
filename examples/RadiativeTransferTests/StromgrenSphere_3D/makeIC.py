#!/usr/bin/env python3

# ---------------------------------------------------------------------
# Add a single star in the center of a glass distribution
# ---------------------------------------------------------------------

from swiftsimio import Writer
import unyt
import numpy as np
import h5py


glass = h5py.File("glassCube_64.hdf5", "r")
parts = glass["PartType0"]
xp = parts["Coordinates"][:]
h = parts["SmoothingLength"][:]
glass.close()

# replace the particle closest to the center
# by the star
r = np.sqrt(np.sum((0.5 - xp) ** 2, axis=1))
rmin = np.argmin(r)
xs = xp[rmin]
xp = np.delete(xp, rmin, axis=0)
h = np.delete(h, rmin)


unitL = unyt.cm
t_end = 1e-3 * unyt.s
edgelen = unyt.c.to("cm/s") * t_end * 2.0
edgelen = edgelen.to(unitL)
boxsize = unyt.unyt_array([edgelen.v, edgelen.v, edgelen.v], unitL)

xs = unyt.unyt_array(
    [np.array([xs[0] * edgelen, xs[1] * edgelen, xs[2] * edgelen])], unitL
)
xp *= edgelen
h *= edgelen


w = Writer(unit_system=unyt.unit_systems.cgs_unit_system, box_size=boxsize, dimension=3)

w.gas.coordinates = xp
w.stars.coordinates = xs
w.gas.velocities = np.zeros(xp.shape) * (unyt.cm / unyt.s)
w.stars.velocities = np.zeros(xs.shape) * (unyt.cm / unyt.s)
w.gas.masses = np.ones(xp.shape[0], dtype=np.float) * 100 * unyt.g
w.stars.masses = np.ones(xs.shape[0], dtype=np.float) * 100 * unyt.g
w.gas.internal_energy = (
    np.ones(xp.shape[0], dtype=np.float) * (300.0 * unyt.kb * unyt.K) / (unyt.g)
)

w.gas.smoothing_length = h
w.stars.smoothing_length = w.gas.smoothing_length[:1]

w.write("stromgrenSphere-3D.hdf5")
