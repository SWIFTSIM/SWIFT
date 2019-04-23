"""
Creates a redshift-dependent cooling box such that it always has the same
_physical_ density at the given redshift.
"""

from swiftsimio import Writer
from swiftsimio.units import cosmo_units

from unyt import mh, cm, s, K, Mpc, kb
import numpy as np

import h5py

# Physics parameters.
boxsize = 1.0 * Mpc
physical_density = 0.1 * mh / cm ** 3
mu_hydrogen = 0.5  # Fully ionised
temperature = 1e7 * K
gamma = 5.0 / 3.0


def get_coordinates(glass_filename: str) -> np.array:
    """
    Gets the coordinates from the glass file.
    """

    with h5py.File(glass_filename, "r") as handle:
        coordinates = handle["PartType1/Coordinates"][...]

    return coordinates


def generate_ics(redshift: float, filename: str, glass_filename: str) -> None:
    """
    Generate initial conditions for the CoolingRedshiftDependence example.
    """

    scale_factor = 1 / (1 + redshift)
    comoving_boxsize = boxsize / scale_factor

    glass_coordinates = get_coordinates(glass_filename)
    number_of_particles = len(glass_coordinates)

    gas_particle_mass = physical_density * (boxsize ** 3) / number_of_particles

    writer = Writer(cosmo_units, comoving_boxsize)

    writer.gas.coordinates = glass_coordinates * comoving_boxsize

    writer.gas.velocities = np.zeros_like(glass_coordinates) * cm / s

    writer.gas.masses = np.ones(number_of_particles, dtype=float) * gas_particle_mass

    # Leave in physical units; handled by boxsize change.
    writer.gas.internal_energy = (
        np.ones(number_of_particles, dtype=float)
        * 3.0 / 2.0
        * (temperature * kb)
        / (mu_hydrogen * mh)
    )

    writer.gas.generate_smoothing_lengths(boxsize=comoving_boxsize, dimension=3)

    writer.write(filename)

    return


if __name__ == "__main__":
    """
    Sets up the initial parameters.
    """

    import argparse as ap

    parser = ap.ArgumentParser(
        description="""
            Sets up the initial conditions for the cooling test. Takes two
            redshifts, and produces two files: ics_high_z.hdf5 and
            ics_low_z.hdf5.
            """
    )

    parser.add_argument(
        "-a",
        "--high",
        help="The high redshift to generate initial conditions for. Default: 1.0",
        default=1,
        type=float,
    )

    parser.add_argument(
        "-b",
        "--low",
        help="The low redshift to generate initial conditions for. Default: 0.01",
        default=0.01,
        type=float,
    )

    parser.add_argument(
        "-n",
        "--nocosmo",
        help="Generate a non-cosmological box? Default: Truthy",
        default=1,
        type=bool,
    )

    parser.add_argument(
        "-g",
        "--glass",
        help="Glass filename. Default: gravity_glassCube_32.hdf5",
        type=str,
        default="gravity_glassCube_32.hdf5",
    )

    args = parser.parse_args()

    generate_ics(args.low, filename="ics_low_z.hdf5", glass_filename=args.glass)
    generate_ics(args.high, filename="ics_high_z.hdf5", glass_filename=args.glass)

    if args.nocosmo:
        generate_ics(0.0, filename="ics_no_z.hdf5", glass_filename=args.glass)

    exit(0)
