"""
Creates the BCC ICs for the blob test.
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

    values = np.linspace(0.0, side_length, num_on_side + 1)[:-1]

    positions = np.empty((num_on_side ** 3, 3), dtype=float)

    for x in range(num_on_side):
        for y in range(num_on_side):
            for z in range(num_on_side):
                index = x * num_on_side + y * num_on_side ** 2 + z

                positions[index, 0] = values[x]
                positions[index, 1] = values[y]
                positions[index, 2] = values[z]

    return positions


def generate_bcc_lattice(num_on_side, side_length=1.0):
    num_per_cube = num_on_side // 2

    cube = generate_cube(num_per_cube, side_length)

    mips = side_length / num_per_cube

    positions = np.concatenate([cube, cube + mips * 0.5])

    return positions


def generate_outside_of_blob(
    num_on_side, num_copies_x, side_length, blob_centre_x, blob_radius, initial_velocity
):
    """
    Generates the area outside of the blob, including a
    cut out for where the blob will fit in.
    """

    bcc_lattice = generate_bcc_lattice(num_on_side)

    # We now need to duplicate.
    bcc_lattice_full = np.concatenate(
        [bcc_lattice + x * np.array([1.0, 0.0, 0.0]) for x in range(num_copies_x)]
    )

    # Now size it appropriately
    bcc_lattice_full *= side_length

    # Now we need to chop out the region that our blob will live in
    dx = bcc_lattice_full - np.array(
        [blob_centre_x, 0.5 * side_length, 0.5 * side_length]
    )
    squared_radius = np.sum(dx * dx, axis=1)

    cut_out_mask = squared_radius > blob_radius ** 2

    positions = bcc_lattice_full[cut_out_mask]

    # Now we can generate the velocities
    velocities = np.zeros_like(positions)
    velocities[:, 0] = initial_velocity

    return positions, velocities


def generate_blob(num_on_side, side_length, blob_centre_x, blob_radius):
    """
    Generate the positions and velocities for the blob.
    """

    bcc_lattice = generate_bcc_lattice(num_on_side)

    # Update to respect side length
    bcc_lattice *= side_length

    # Find the radii
    squared_radius = np.sum((bcc_lattice - 0.5 * side_length) ** 2, axis=1)
    inside_sphere = squared_radius < blob_radius * blob_radius

    # Now select out particles
    bcc_lattice = bcc_lattice[inside_sphere]

    # Move to the correct x_position
    bcc_lattice[:, 0] = bcc_lattice[:, 0] + blob_centre_x - 0.5 * side_length
    positions = bcc_lattice

    # Generate velocities
    velocities = np.zeros_like(positions)

    return positions, velocities


def write_out_ics(
    filename,
    num_on_side,
    side_length=1.0,
    duplications=4,
    blob_radius=0.1,
    blob_location_x=0.5,
    inside_density=10,
    velocity_outside=2.7,  # Actually the mach number
):
    x = Writer(
        cgs_unit_system, [side_length * duplications, side_length, side_length] * cm
    )

    positions_outside, velocities_outside = generate_outside_of_blob(
        num_on_side,
        duplications,
        side_length,
        blob_location_x,
        blob_radius,
        initial_velocity=1.0,
    )
    positions_inside, velocities_inside = generate_blob(
        int(num_on_side * np.cbrt(inside_density)),
        side_length,
        blob_location_x,
        blob_radius,
    )

    coordinates = np.concatenate([positions_inside, positions_outside])

    velocities = np.concatenate([velocities_inside, velocities_outside])

    x.gas.coordinates = coordinates * cm

    x.gas.velocities = velocities * cm / s

    x.gas.masses = (
        np.ones(coordinates.shape[0], dtype=float) * (side_length / num_on_side) * g
    )

    # Set to be (1 / n) c_s
    x.gas.internal_energy = (
        np.ones(coordinates.shape[0], dtype=float)
        * ((1.0 / velocity_outside) * x.gas.velocities[-1][0]) ** 2
        / (5.0 / 3.0 - 1)
    )

    # Need to correct the particles inside the sphere so that we are in pressure equilibrium everywhere
    x.gas.internal_energy[: len(positions_inside)] /= inside_density

    x.gas.generate_smoothing_lengths(
        boxsize=(duplications / 2) * side_length * cm, dimension=3
    )

    x.write(filename)

    return


if __name__ == "__main__":
    write_out_ics(filename="blob.hdf5", num_on_side=64)
