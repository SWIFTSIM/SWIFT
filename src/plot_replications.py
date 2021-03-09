#!/bin/env python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations

import subprocess


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


if __name__ == "__main__":

    boxsize = 50
    observer_position = (25, 25, 25)
    lc_rmin = 240
    lc_rmax = 250
    view_vector = (1.0, 0.0, 0.0)
    view_radius = np.radians(20.0)

    # Run the code    
    command = (["./test_replications", str(boxsize)] +
               [str(obs_pos) for obs_pos in observer_position] +
               [str(lc_rmin), str(lc_rmax)] +
               [str(v) for v in view_vector] + [str(view_radius),])

    process = subprocess.run(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)

    # Extract stdout to numpy arrays
    data = np.loadtxt((line for line in process.stdout.decode().split("\n")), delimiter=",")
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    rmin = data[:,3]
    rmax = data[:,4]
    nrep = data.shape[0]
    print(nrep)

    # Plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect("equal")

    # Indicate maximum radius
    draw_sphere(observer_position[0],
                observer_position[1],
                observer_position[2],
                lc_rmax)

    # Draw replications
    for i in range(nrep):
        draw_cube(x[i]*boxsize+0.5*boxsize,
                  y[i]*boxsize+0.5*boxsize,
                  z[i]*boxsize+0.5*boxsize,
                  boxsize)

    plt.show()

