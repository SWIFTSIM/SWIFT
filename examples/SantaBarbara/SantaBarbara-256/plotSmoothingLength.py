"""
Plots the smoothing length (compared to the softening).
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py

from collections import namedtuple

SnapshotData = namedtuple(
    "SnapshotData",
    [
        "smoothing_lengths",
        "particle_ids",
        "softening",
        "internal_length",
        "snapshot_length",
    ],
)

HaloCatalogueData = namedtuple(
    "HaloCatalogueData", ["largest_halo", "particle_ids_in_largest_halo"]
)


def load_data(filename: str) -> SnapshotData:
    """
    Loads the data that we need, i.e. the smoothing lengths and the
    softening length, from the snapshot.
    """

    with h5py.File(filename, "r") as handle:
        data = SnapshotData(
            smoothing_lengths=handle["PartType0/SmoothingLength"][...],
            particle_ids=handle["PartType0/ParticleIDs"][...],
            softening=handle["GravityScheme"].attrs[
                "Comoving softening length [internal units]"
            ][0],
            internal_length=handle["InternalCodeUnits"].attrs[
                "Unit length in cgs (U_L)"
            ][0],
            snapshot_length=handle["Units"].attrs["Unit length in cgs (U_L)"][0],
        )

    return data


def load_halo_data(halo_filename: str) -> HaloCatalogueData:
    """
    Loads the halo data and finds the particle IDs that belong to
    the largest halo. The halo filename should be given without
    any extension as we need a couple of files to complete this.
    """

    catalogue_filename = f"{halo_filename}.properties"
    groups_filename = f"{halo_filename}.catalog_groups"
    particles_filename = f"{halo_filename}.catalog_particles"

    with h5py.File(catalogue_filename, "r") as handle:
        largest_halo = np.where(
            handle["Mass_200crit"][...] == handle["Mass_200crit"][...].max()
        )[0][0]

    with h5py.File(groups_filename, "r") as handle:
        offset_begin = handle["Offset"][largest_halo]
        offset_end = handle["Offset"][largest_halo + 1]

    with h5py.File(particles_filename, "r") as handle:
        particle_ids = handle["Particle_IDs"][offset_begin:offset_end]

    return HaloCatalogueData(
        largest_halo=largest_halo, particle_ids_in_largest_halo=particle_ids
    )


def make_plot(
    snapshot_filename: str,
    halo_filename: str,
    output_filename="smoothing_length_variation.png",
) -> None:
    """
    Makes the plot and saves it in output_filename.

    The halo filename should be provided without extension.
    """

    data = load_data(filename)
    halo_data = load_halo_data(halo_filename)

    smoothing_lengths_in_halo = data.smoothing_lengths[
        np.in1d(data.particle_ids, halo_data.particle_ids_in_largest_halo)
    ]

    softening = data.softening * (data.snapshot_length / data.internal_length)

    fig, ax = plt.subplots(1)

    ax.semilogy()

    ax.hist(data.smoothing_lengths, bins="auto", label="All particles")
    ax.hist(
        smoothing_lengths_in_halo,
        bins="auto",
        label=f"Particles in largest halo (ID={halo_data.largest_halo})",
    )
    ax.axvline(x=softening, label="Softening", ls="--", color="grey")

    ax.legend()

    ax.set_xlabel("Smoothing length")
    ax.set_ylabel("Number of particles")

    ax.set_xlim(0, ax.get_xlim()[1])

    fig.tight_layout()

    fig.savefig(output_filename, dpi=300)

    return


if __name__ == "__main__":
    import argparse as ap

    PARSER = ap.ArgumentParser(
        description="""
            Makes a plot of the smoothing lengths in the box, compared
            to the gravitational softening. Also splits out the particles
            that are contained in the largest halo according to the
            velociraptor outputs.
            """
    )

    PARSER.add_argument(
        "-s",
        "--snapshot",
        help="""
            Filename and path for the snapshot (without the .hdf5),
            Default: ./santabarbara_0153
            """,
        required=False,
        default="./santabarbara_0153",
    )

    PARSER.add_argument(
        "-v",
        "--velociraptor",
        help="""
            The filename and path of the velociraptor files, excluding
            the descriptors (i.e. without .catalog_particles).
            Default: ./halo/santabarbara
            """,
        required=False,
        default="./halo/santabarbara",
    )

    ARGS = vars(PARSER.parse_args())

    filename = f"{ARGS['snapshot']}.hdf5"

    make_plot(filename, ARGS["velociraptor"])
