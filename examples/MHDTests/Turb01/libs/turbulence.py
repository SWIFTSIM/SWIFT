from __future__ import print_function

from sys import exit
from numpy import meshgrid, sqrt, log, exp, zeros, linspace, array, cross
from numpy import fft, random
from numpy import pi
from time import time
from scipy.interpolate import RegularGridInterpolator


class VelocityGrid:
    """ Class for creating a 3-D grid with turbulent velocity field.
        The velocities are produced from a Gaussian random distribution
        with a given power spectrum.
        Each velocity component is sampled in Fourier space, and then
        transformed back to real space.
        Based on Dubinski et al. (1995).

        Arguments:
           npow : power index of spectrum.
           ngrid: number of grid points per dimension (must be even).
           xmax : outer scale of turbulence.
           dx   : physical separation between neighboring points.
           seed : number that determines the random realization.
    """

    def __init__(self, npow=-4.0, ngrid=256, xmax=1.0, dx=0.01, seed=2021987):

        start = time()
        print("Creating 3-D velocity grid with power spectrum P_k~k**{}".format(npow))

        if ngrid % 2 != 0:
            print("Grid points must be an even number. Exiting.")
            exit()
        nc = int(ngrid / 2) + 1

        kmax = 2 * pi / dx
        kmin = 2 * pi / xmax

        kx = fft.fftfreq(ngrid, d=1 / (2 * kmax))
        ky = kx
        kz = fft.rfftfreq(ngrid, d=1 / (2 * kmax))

        # we produce a 3-D grid of the Fourier coordinates
        kxx, kyy, kzz = meshgrid(kx, ky, kz, indexing="ij", sparse=True)
        kk = kxx * kxx + kyy * kyy + kzz * kzz + kmin ** 2

        random.seed(seed)

        # we sample the components of a vector potential, as we want
        # an incompresible velocity field
        xi1 = random.random(size=kk.shape)
        xi2 = random.random(size=kk.shape)
        c = kk ** ((npow - 2.0) / 4.0) * sqrt(-log(1 - xi1))
        phi = 2 * pi * xi2
        akx = c * exp(1j * phi)
        xi1 = random.random(size=kk.shape)
        xi2 = random.random(size=kk.shape)
        c = kk ** ((npow - 2.0) / 4.0) * sqrt(-log(1 - xi1))
        phi = 2 * pi * xi2
        aky = c * exp(1j * phi)
        xi1 = random.random(size=kk.shape)
        xi2 = random.random(size=kk.shape)
        c = kk ** ((npow - 2.0) / 4.0) * sqrt(-log(1 - xi1))
        phi = 2 * pi * xi2
        akz = c * exp(1j * phi)

        new_shape = akx.shape + (3,)
        kv = zeros(new_shape, dtype=akx.dtype)
        kv[:, :, :, 0] = 1j * kxx
        kv[:, :, :, 1] = 1j * kyy
        kv[:, :, :, 2] = 1j * kzz
        ak = zeros(new_shape, dtype=akx.dtype)
        ak[:, :, :, 0] = akx
        ak[:, :, :, 1] = aky
        ak[:, :, :, 2] = akz

        # the velocity vector in Fourier space is obtained by
        # taking the curl of A, which is
        vk = cross(kv, ak)

        self.ngrid = ngrid
        self.vx = fft.irfftn(vk[:, :, :, 0])
        self.vy = fft.irfftn(vk[:, :, :, 1])
        self.vz = fft.irfftn(vk[:, :, :, 2])

        print("\nInverse Fourier Transform took {:g}s.".format(time() - start))

    def coordinate_grid(self, xstart=0.0, xend=1.0):
        self.x = linspace(xstart, xend, self.ngrid)

    def add_turbulence(self, pos, vel):

        pos = array(pos).reshape(-1, 3)
        vel = array(vel).reshape(-1, 3)

        if not hasattr(self, "x"):
            print("WARNING: Interpolating with grid with default values.")
            print("         Please make sure this is what you want.")
            self.coordinate_grid()

        x = self.x
        vx = self.vx
        vy = self.vy
        vz = self.vz

        interp_func_x = RegularGridInterpolator((x, x, x), vx)
        interp_func_y = RegularGridInterpolator((x, x, x), vy)
        interp_func_z = RegularGridInterpolator((x, x, x), vz)
        vx = interp_func_x(pos)
        vy = interp_func_y(pos)
        vz = interp_func_z(pos)

        vel[:, 0] += vx
        vel[:, 1] += vy
        vel[:, 2] += vz

        return vel
