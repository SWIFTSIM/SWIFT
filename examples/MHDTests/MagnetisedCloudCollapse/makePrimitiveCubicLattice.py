"""
Create a cubic lattice.
"""

import numpy as np
import h5py

from unyt import cm, g, s, erg
from unyt.unit_systems import cgs_unit_system

from swiftsimio import Writer


def generate_cube(num_on_side, side_length=1.0):
    """
    Generates a cube
    """

    values, step = np.linspace(0.0, side_length, num_on_side, endpoint=False, retstep=True)
    values += 0.5*step

    positions = np.empty((num_on_side ** 3, 3), dtype=float)

    for x in range(num_on_side):
        for y in range(num_on_side):
            for z in range(num_on_side):
                index = x * num_on_side + y * num_on_side ** 2 + z

                positions[index, 0] = values[x]
                positions[index, 1] = values[y]
                positions[index, 2] = values[z]

    return positions


def write_out_glass(filename, cube, side_length=1.0):
    x = Writer(cgs_unit_system, side_length * cm)

    x.gas.coordinates = cube * cm

    x.gas.velocities = np.zeros_like(cube) * cm / s

    x.gas.masses = np.ones(cube.shape[0], dtype=float) * g

    x.gas.internal_energy = np.ones(cube.shape[0], dtype=float) * erg / g

    x.gas.generate_smoothing_lengths(boxsize=side_length * cm, dimension=3)

    x.write(filename)

    return


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(description="Generate a BCC lattice")

    parser.add_argument(
        "-n",
        "--numparts",
        help="Number of particles on a side. Default: 32",
        default=32,
        type=int,
    )

    args = parser.parse_args()

    output = "CubicPrimitive_{}.hdf5".format(args.numparts)

    glass_cube = generate_cube(args.numparts)
    write_out_glass(output, glass_cube)
