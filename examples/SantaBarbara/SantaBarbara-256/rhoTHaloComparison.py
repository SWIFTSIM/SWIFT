"""
This script finds the temperatures inside all of the halos and
compares it against the virial temperature. This uses velociraptor
and the SWIFT snapshot.

Folkert Nobels (2018) nobels@strw.leidenuniv.nl
Josh Borrow (2019) joshua.borrow@durham.ac.uk
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

mH = 1.6733e-24  # g
kB = 1.38e-16  # erg/K


def virial_temp(mu, M, h=0.703, a=1.0):
    """
    Calculates the virial temperature according to

    https://arxiv.org/pdf/1105.5701.pdf

    Equation 1.
    """
    return  4e4 * (mu / 1.2) * (M * h / 1e8) ** (2 / 3) / (10 * a)


def calculate_group_sizes_array(offsets: np.array, total_size: int) -> np.array:
    """
    Calculates the group sizes array from the offsets and total size, i.e. it
    calculates the diff between all of the offsets.
    """

    # Does not include the LAST one
    group_sizes = [x - y for x, y in zip(offsets[1:], offsets[:-1])]
    group_sizes += [total_size - offsets[-1]]
    group_sizes = np.array(group_sizes, dtype=type(offsets[0]))

    return group_sizes


def create_group_array(group_sizes: np.array) -> np.array:
    """
    Creates an array that looks like:
    [GroupID0, GroupID0, ..., GroupIDN, GroupIDN]
    i.e. for each group create the correct number of group ids.
    This is used to be sorted alongside the particle IDs to track
    the placement of group IDs.
    """

    slices = []
    running_total_of_particles = type(offsets[0])(0)

    for group in group_sizes:
        slices.append([running_total_of_particles, group + running_total_of_particles])
        running_total_of_particles += group

    groups = np.empty(group_sizes.sum(), dtype=int)

    for group_id, group in enumerate(slices):
        groups[group[0] : group[1]] = group_id

    return groups


if __name__ == "__main__":
    import argparse as ap

    PARSER = ap.ArgumentParser(
        description="""
        Makes many plots comparing the virial temperature and the
        temperature of halos. Requires the velociraptor files and
        the SWIFT snapshot.
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

    # Grab some metadata before we begin.
    with h5py.File("%s.hdf5" % ARGS["snapshot"], "r") as handle:
        # Cosmology
        a = handle["Cosmology"].attrs["Scale-factor"][0]
        h = handle["Cosmology"].attrs["h"][0]

        # Gas
        hydro = handle["HydroScheme"].attrs
        X = hydro["Hydrogen mass fraction"][0]
        gamma = hydro["Adiabatic index"][0]
        mu = 1 / (X + (1 - X) / 4)

    # First we must construct a group array so we know which particles belong
    # to which group.
    with h5py.File("%s.catalog_groups" % ARGS["velociraptor"], "r") as handle:
        offsets = handle["Offset"][...]

    # Then, extract the particles that belong to the halos. For that, we need
    # the particle IDs:
    with h5py.File("%s.catalog_particles" % ARGS["velociraptor"], "r") as handle:
        ids_in_halos = handle["Particle_IDs"][...]

    number_of_groups = len(offsets)
    group_sizes = calculate_group_sizes_array(offsets, ids_in_halos.size)
    group_array = create_group_array(group_sizes)

    # We can now load the particle data from the snapshot.
    with h5py.File("%s.hdf5" % ARGS["snapshot"], "r") as handle:
        gas_particles = handle["PartType0"]

        particle_ids = gas_particles["ParticleIDs"][...]

        # Requires numpy 1.15 or greater.
        _, particles_in_halos_mask, group_array_mask = np.intersect1d(
            particle_ids,
            ids_in_halos,
            assume_unique=True,
            return_indices=True,
        )

        # We also need to re-index the group array to cut out DM particles
        group_array = group_array[group_array_mask]
        
        # Kill the spare
        del particle_ids

        # Now we can only read the properties that we require from the snapshot!
        temperatures = np.take(gas_particles["InternalEnergy"], particles_in_halos_mask)
        # This 1e10 should probably be explained somewhere...
        temperatures *= 1e10 * (gamma - 1) * mu * mH / kB

        densities = np.take(gas_particles["Density"], particles_in_halos_mask)

    # Just a quick check to make sure nothing's gone wrong.
    assert len(group_array) == len(temperatures)

    # Now we can loop through all the particles and find out the mean temperature and
    # density in each halo.

    particles_in_group = np.zeros(number_of_groups, dtype=int)
    temp_in_group = np.zeros(number_of_groups, dtype=float)
    dens_in_group = np.zeros(number_of_groups, dtype=float)

    for group, T, rho in zip(group_array, temperatures, densities):
        particles_in_group[group] += 1
        temp_in_group[group] += T
        dens_in_group[group] += rho

    # First get a mask to ensure no runtime warnings
    mask = particles_in_group != 0
    
    # Normalize
    temp_in_group[mask] /= particles_in_group[mask]
    dens_in_group[mask] /= particles_in_group[mask]

    # Now we can load the data according to the halo finder to compare with.
    with h5py.File("%s.properties" % ARGS["velociraptor"], "r") as handle:
        halo_masses = handle["Mass_200crit"][...]

    halo_temperatures = virial_temp(mu, halo_masses * 1e10, h=h, a=a)

    # Finally, the plotting!

    fig, ax = plt.subplots()
    ax.loglog()

    mask = np.logical_and.reduce([
         halo_temperatures != 0.0,
         temp_in_group != 0.0,
    ])

    temp_in_group = temp_in_group[mask]
    halo_temperatures = halo_temperatures[mask]

    mean_range = [temp_in_group.min(), temp_in_group.max()]
    halo_range = [halo_temperatures.min(), halo_temperatures.max()]

    bottom = min([halo_range[0], mean_range[0]])
    top = max([halo_range[1], mean_range[1]])

    plt.plot(
        [bottom, top],
        [bottom, top],
        lw=2, linestyle="--", color="grey", label="1:1"
    )
    
    ax.scatter(halo_temperatures, temp_in_group, s=2, edgecolor="none", label="Halos")

    ax.set_ylabel("Mean Group Temperature [K]")
    ax.set_xlabel("Halo Virial Temperature [K]")

    ax.set_ylim(mean_range)
    ax.set_xlim(halo_range)

    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig("temperature_comparison.png", dpi=300)
