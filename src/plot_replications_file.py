#!/bin/env python

import subprocess
import sys

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations


def draw_cube(x, y, z, r):
    """
    Draw a cube of side r centred on (x,y,z)
    """
    r = [-0.5*r, 0.5*r]
    for s, e in combinations(np.array(list(product(r, r, r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            pos = np.asarray(list(zip(s, e)))
            pos[0,:] += x
            pos[1,:] += y
            pos[2,:] += z
            ax.plot3D(pos[0,:], pos[1,:], pos[2,:], color="b")


def draw_sphere(x, y, z, r):
    """
    Draw a sphere of radius r centred on (x,y,z)
    """
    nres = 10
    u, v = np.mgrid[0:2*np.pi:2j*nres, 0:np.pi:1j*nres]
    x = x + np.cos(u)*np.sin(v) * r
    y = y + np.sin(u)*np.sin(v) * r
    z = z + np.cos(v) * r
    ax.plot_wireframe(x, y, z, color="r")


def read_replications(fname):
    """
    Read replications file written by swift
    """
    with open(fname, "r") as infile:
        infile.readline() # skip header
        obs_x, obs_y, obs_z = [float(field) for field in infile.readline().split(",")]
        infile.readline() # skip header
        boxsize, rmin, rmax = [float(field) for field in infile.readline().split(",")]
        infile.readline() # skip header
        data = np.loadtxt(infile, dtype=float, delimiter=",")
    rep_pos = data[:,0:3]
    rep_rmin2 = data[:,3]
    rep_rmax2 = data[:,4]
    obs_pos = np.asarray((obs_x, obs_y, obs_z), dtype=float)
    
    return obs_pos, boxsize, rmin, rmax, rep_pos, rep_rmin2, rep_rmax2


if __name__ == "__main__":

    # Read the input file
    (obs_pos, boxsize, rmin, rmax,
     rep_pos, rep_rmin2, rep_rmax2) = read_replications(sys.argv[1])
    nrep = rep_pos.shape[0]

    # Set up 3D plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect("equal")

    # Indicate inner and outer radii
    draw_sphere(obs_pos[0],
                obs_pos[1],
                obs_pos[2],
                rmin)
    draw_sphere(obs_pos[0],
                obs_pos[1],
                obs_pos[2],
                rmax)

    # Draw replications
    for i in range(nrep):
        draw_cube(rep_pos[i,0]+0.5*boxsize,
                  rep_pos[i,1]+0.5*boxsize,
                  rep_pos[i,2]+0.5*boxsize,
                  boxsize)

    plt.show()

