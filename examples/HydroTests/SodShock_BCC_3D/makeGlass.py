"""
Create a fcc glass.
"""

import numpy as np
import h5py

from unyt import cm, g, s, erg
from unyt.unit_systems import cgs_unit_system

import swiftsimio as sw


def generate_cube(num_on_side, side_length=1.0):
    """
    Generates a cube
    """

    values = np.linspace(0.0, side_length, num_on_side + 1)[:-1]

    positions = np.empty((num_on_side**3, 3), dtype=float)

    for x in range(num_on_side):
        for y in range(num_on_side):
            for z in range(num_on_side):
                index = x * num_on_side + y * num_on_side**2 + z

                positions[index, 0] = values[x]
                positions[index, 1] = values[y]
                positions[index, 2] = values[z]

    return positions


def generate_two_cube(num_on_side, side_length=1.0):
    cube = generate_cube(num_on_side // 2, side_length)

    mips = side_length / num_on_side

    positions = np.concatenate([cube, cube + mips])

    return positions


def write_out_glass(filename, cube, side_length=1.0):
    boxsize_cosmo = sw.cosmo_array(
        [side_length, side_length, side_length],
        cm,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=1,
    )
    x = sw.Writer(unit_system=cgs_unit_system, boxsize=boxsize_cosmo)

    x.gas.coordinates = sw.cosmo_array(
        cube, cm, comoving=True, scale_factor=1.0, scale_exponent=1
    )
    x.gas.velocities = sw.cosmo_array(
        np.zeros_like(cube), cm / s, comoving=True, scale_factor=1.0, scale_exponent=0
    )
    x.gas.masses = sw.cosmo_array(
        np.ones(cube.shape[0], dtype=float),
        g,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=0,
    )
    x.gas.internal_energy = sw.cosmo_array(
        np.ones(cube.shape[0], dtype=float),
        erg / g,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=-2,
    )
    x.gas.generate_smoothing_lengths()

    x.write(filename)

    return


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(description="Generate a glass file with a FCC lattice")

    parser.add_argument(
        "-n",
        "--numparts",
        help="Number of particles on a side. Default: 64",
        default=1,
        type=int,
    )

    parser.add_argument(
        "-o", "--output", help="Output filename.", type=str, required=False
    )

    args = parser.parse_args()

    glass_cube = generate_two_cube(args.numparts)
    write_out_glass(args.output, glass_cube)
