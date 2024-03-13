"""
Create a close packed glass.
"""

import numpy as np
import h5py

from unyt import cm, g, s, erg
from unyt.unit_systems import cgs_unit_system

import matplotlib.pyplot as plt

from swiftsimio import Writer


def generate_close_packed_cell(num_on_side_x, num_on_side_y, num_on_side_z, side_length_x=1.0):
    """
    Generates a cube
    """

    slx = side_length_x
    sly = np.sqrt(3) * slx * num_on_side_y / (2 * num_on_side_x)
    slz = np.sqrt(2/3) * slx * num_on_side_z / num_on_side_x

    print("The shape of the unit cell is : ", slx, sly, slz)

    valx = np.linspace(0.0, slx, num_on_side_x + 1)[:-1]
    valy = np.linspace(0.0, sly, num_on_side_y + 1)[:-1]
    valz = np.linspace(0.0, slz, num_on_side_z + 1)[:-1]

    positions = np.zeros((num_on_side_x * num_on_side_y * num_on_side_z, 3), dtype=float)

    print("Positions array has shape : ", positions.shape) 

    slxi = slx / num_on_side_x
    tol = 1.0 #0.9999999

    for x in range(num_on_side_x):
        for y in range(num_on_side_y):
            for z in range(num_on_side_z):                               
                 
                    index = z*num_on_side_x*num_on_side_y + y * num_on_side_x + x

                    positions[index, 0] = (valx[x] + (y % 2) * 0.5 * slxi + (z % 3) * 0.5 * slxi) % (tol * slx)
                    positions[index, 1] = (valy[y] + (z % 3) * slxi / (2 * np.sqrt(3))) % (tol * sly) 
                    positions[index, 2] = valz[z] % (tol * slz)

    positions += 0.5 * np.array([0.5 * slxi, 0.5 * np.sqrt(3) * slxi, np.sqrt(2/3) * slxi])                    
    
    print("The maximum x is : ", max(positions[:,0]))             

    '''
    unx = np.unique(positions[:,0])
    uny = np.unique(positions[:,1])
    unz = np.unique(positions[:,2])

    c = ['k', 'r', 'b', 'g']
    
    masks_x = [positions[:,0] == unx[0], positions[:,0] == unx[1], positions[:,0] == unx[2]]
    masks_y = [positions[:,1] == uny[0], positions[:,1] == uny[1], positions[:,1] == uny[2]]
    masks_z = [positions[:,2] == unz[0], positions[:,2] == unz[1], positions[:,2] == unz[2]]
    
    masks = [masks_x, masks_y, masks_z]
    
    masks = [positions[:,2] == unz[0], positions[:,2] == unz[1], positions[:,2] == unz[2], positions[:,2] == unz[3]]

    fig, ax = plt.subplots(3,1)

    for ii in range(3):
        for ind, ci in enumerate(c):
            ax[ii].scatter(positions[:,ii][masks[ind]], positions[:,(ii+1)%3][masks[ind]], c=ci)
            ax[ii].set_aspect('equal')
    '''
    
    fig, ax = plt.subplots(3,1)

    for ii in range(3):
        ax[ii].scatter(positions[:,ii], positions[:,(ii+1)%3])
        ax[ii].set_aspect('equal')         
    
    plt.tight_layout()        
    plt.savefig("Lattice.png")
    
    return positions

def write_out_glass(filename, cube, num_on_side_x, num_on_side_y, num_on_side_z, side_length_x=1.0):
    
    slx = side_length_x
    sly = np.sqrt(3) * slx * num_on_side_y / (2 * num_on_side_x)
    slz = np.sqrt(2/3) * slx * num_on_side_z / num_on_side_x

    boxSize = [slx, sly, slz] * cm

    x = Writer(cgs_unit_system, boxSize)

    x.gas.coordinates = cube * cm

    x.gas.velocities = np.zeros_like(cube) * cm / s

    x.gas.masses = np.ones(cube.shape[0], dtype=float) * g

    x.gas.internal_energy = np.ones(cube.shape[0], dtype=float) * erg / g

    x.gas.generate_smoothing_lengths(boxsize=boxSize, dimension=3)

    x.write(filename)

    return


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(description="Generate a CP lattice")

    parser.add_argument(
        "-nx",
        "--numpartsx",
        help="Number on base unit cells replicated along the x axis. Default: 2",
        default=2,
        type=int,
    )

    parser.add_argument(
        "-ny",
        "--numpartsy",
        help="Number on base unit cells replicated along the y axis. Default: 2",
        default=2,
        type=int,
    )

    parser.add_argument(
        "-nz",
        "--numpartsz",
        help="Number on base unit cells replicated along the z axis. Default: 2",
        default=2,
        type=int,
    )

    args = parser.parse_args()

    output = "CPLglassCube_{}x_{}y_{}z.hdf5".format(args.numpartsx, args.numpartsy, args.numpartsz)

    glass_cube = generate_close_packed_cell(args.numpartsx, args.numpartsy, args.numpartsz)
    write_out_glass(output, glass_cube, args.numpartsx, args.numpartsy, args.numpartsz)
