from __future__ import print_function

from numpy import pi, sum, max, min, log, exp
from numpy import diff, where, unique, argsort
from numpy import array, full, linspace, mgrid
from numpy import transpose, append, digitize
from numpy import concatenate, ones, meshgrid
from numpy import zeros
from numpy.linalg import norm


class Cube:
    """ Class for creating a distribution of particles in a Cube.

        Arguments:
            n     : total number of desired points to represent the Cube
            side  : side (with units)
            mass  : mass (with units)
    """

    def __init__(self, n=32, side=1.0, mass=1.0):

        nside = int(n)

        naux = nside ** 3
        h = ones(naux) / 2.0 * side / nside  # min particle separation

        # z3     = mgrid[0:nside, 0:nside, 0:nside].T.reshape(naux, 3)
        ##grid   = (concatenate((z3, z3+[0.5, 0.5, 0], z3+[0, 0.5, 0.5],
        ##          z3+[0.5, 0, 0.5])) + 0.25)/nside * side
        # grid   = (concatenate((z3+[0.5, 0.5, 0], z3+[0, 0.5, 0.5],
        #          z3+[0.5, 0, 0.5])))/nside * side

        # pos    = grid
        x = linspace(0.0, 1.0, nside + 1)
        x = 0.5 * (x[1:] + x[:-1])
        y = x
        z = x
        xx, yy, zz = meshgrid(x, y, z)
        pos = zeros((naux, 3))
        pos[:, 0] = xx.reshape((naux))
        pos[:, 1] = yy.reshape((naux))
        pos[:, 2] = zz.reshape((naux))

        npart = len(pos)
        masses = ones(npart) * mass / float(npart)  # uniform masses
        print("We placed {:d} gas cells in a cube.".format(npart))

        self.npart = npart
        self.dx = h
        self.pos = pos
        self.lside = side
        self.mass = masses
