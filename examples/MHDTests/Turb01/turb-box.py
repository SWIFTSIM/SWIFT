from __future__ import print_function

import numpy as np
from libs.turbulence import VelocityGrid
from libs.uniform_cube import Cube
from libs.rotation import Rotation
from libs.const import G, msol, parsec
from libs.utils import save_particles
from libs.options_parser import OptionsParser


if __name__ == "__main__":

    op = OptionsParser()
    args = op.get_args()

    # mbox   = args.mass * 10E10 * msol
    # lside  = args.lside * parsec * 10E6
    mbox = args.mass  # * 10E10 * msol
    lside = args.lside  # * parsec * 10E6
    nside = args.nside
    print("We want {:f} gas cells per side in the box".format(nside))
    print("the Lenght of the box is {:f}".format(lside))

    # r_com = np.array([0.,0.,0.])

    # first, determine position of particles given total number of
    # desired cells
    box = Cube(n=nside, side=lside, mass=mbox)

    pos = box.pos
    dx = box.dx
    ngas = int(box.npart)
    mass = box.mass
    print("the h {:f}".format(dx[1]))
    print("the Npart {:d}".format(ngas))
    print("the mass {:f}".format(mass[1]))

    Gamma = 5.0 / 3.0
    P0 = 1.0
    Rho0 = 1.0

    ids = np.arange(1, ngas + 1)
    u = np.ones(ngas) * P0 / (Rho0 * (Gamma - 1.0))
    h = dx
    vel = np.zeros((ngas, 3))

    # produce the velocity grid for turbulent ICs
    vg = VelocityGrid(xmax=lside, dx=dx[1], npow=args.npow, ngrid=int(args.ngrid))
    vg.coordinate_grid(xstart=0.0, xend=lside)
    print("Adding turbulent velocity to particles.")
    vel = vg.add_turbulence(pos=pos, vel=vel)

    # now we need to normalize the velocity values
    # we do it according to the alpha value
    vtur = vel - np.mean(vel, axis=0)
    vt2 = np.linalg.norm(vtur, axis=1) ** 2
    etur = 0.5 * np.sum(mass * vt2)
    epot = 3.0 / 5.0 * G * np.sum(mass) ** 2 / lside  # we think the mass as a sphere
    kvel = np.sqrt(args.alpha * epot / etur)
    print("the Kvel {:f}".format(kvel))
    # kvel = args.alpha
    vel *= kvel

    # we manually add rotation if desired
    # rot = Rotation(beta=args.beta, alpha=args.alpha, epot=epot)
    # vel = rot.add_rotation(pos=pos, vel=vel, mass=mass)

    print("Writing output file {}...".format(args.outfile))
    save_particles(ids, pos, vel, mass, u, args.outfile, args.format, lside, h)

    print("done...bye!")
