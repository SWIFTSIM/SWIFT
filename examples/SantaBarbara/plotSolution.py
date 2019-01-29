"""
Plots the "solution" (i.e. some profiles) for the Santa Barbara cluster.

Invoke as follows:

python3 plotSolution.py <snapshot number> <catalogue directory> <number of bins (optional)>
"""

import matplotlib.pyplot as plt
import numpy as np

import h5py

from collections import namedtuple
from typing import Tuple

try:
    import makeImage

    create_images = True
except:
    create_images = False

# Simulation data
SimulationParticleData = namedtuple(
    "SimulationData", ["gas", "dark_matter", "metadata"]
)
ParticleData = namedtuple(
    "ParticleData", ["coordinates", "radii", "masses", "densities", "energies"]
)
MetaData = namedtuple("MetaData", ["header", "code", "hydroscheme"])
HaloData = namedtuple("HaloData", ["c", "Rvir", "Mvir", "center"])


def get_energies(handle: h5py.File):
    """
    Gets the energies with the correct units.
    """
    u = handle["PartType0/InternalEnergy"][:]
    unit_length_in_cgs = handle["/Units"].attrs["Unit length in cgs (U_L)"]
    unit_mass_in_cgs = handle["/Units"].attrs["Unit mass in cgs (U_M)"]
    unit_time_in_cgs = handle["/Units"].attrs["Unit time in cgs (U_t)"]
    gas_gamma = handle["/HydroScheme"].attrs["Adiabatic index"][0]
    a = handle["/Cosmology"].attrs["Scale-factor"][0]

    unit_length_in_si = 0.01 * unit_length_in_cgs
    unit_mass_in_si = 0.001 * unit_mass_in_cgs
    unit_time_in_si = unit_time_in_cgs

    u *= unit_length_in_si ** 2 / unit_time_in_si ** 2
    u /= a ** (3 * (gas_gamma - 1.))

    return u


def load_data(filename: str, center: np.array) -> SimulationParticleData:
    """
    Loads the relevant data for making the profiles, as well as some metadata
    for the plot.

    Center is the center of the SB cluster and is used to calculate the radial
    distances to the particles.
    """

    with h5py.File(filename, "r") as file:
        gas_handle = file["PartType0"]
        dm_handle = file["PartType1"]

        gas_data = ParticleData(
            coordinates=gas_handle["Coordinates"][...],
            radii=get_radial_distances(gas_handle["Coordinates"][...], center),
            masses=gas_handle["Masses"][...],
            energies=get_energies(file),
            densities=gas_handle["Density"][...],
        )

        dm_data = ParticleData(
            coordinates=dm_handle["Coordinates"][...],
            radii=get_radial_distances(dm_handle["Coordinates"][...], center),
            masses=dm_handle["Masses"][...],
            energies=None,
            densities=None,
        )

        metadata = MetaData(
            header=dict(file["Header"].attrs),
            code=dict(file["Code"].attrs),
            hydroscheme=dict(file["HydroScheme"].attrs),
        )

        simulation_data = SimulationParticleData(
            gas=gas_data, dark_matter=dm_data, metadata=metadata
        )

    return simulation_data


def get_halo_data(catalogue_filename: str) -> HaloData:
    """
    Gets the halo center of the largest halo (i.e. the SB cluster).

    You will want the .properties file, probably

    halo/santabarbara.properties

    that is given by VELOCIraptor.
    """


    with h5py.File(catalogue_filename, "r") as file:
        largest_halo = np.where(
            file["Mass_200crit"][...] == file["Mass_200crit"][...].max()
        )

        x = float(np.take(file["Xc"], largest_halo))
        y = float(np.take(file["Yc"], largest_halo))
        z = float(np.take(file["Zc"], largest_halo))
        Mvir = float(np.take(file["Mass_200crit"], largest_halo))
        Rvir = float(np.take(file["R_200crit"], largest_halo))
        c = float(np.take(file["cNFW"], largest_halo))

    return HaloData(c=c, Rvir=Rvir, Mvir=Mvir, center=np.array([x, y, z]))


def get_radial_distances(coordinates: np.ndarray, center: np.array) -> np.array:
    """
    Gets the radial distances for all particles.
    """
    dx = coordinates - center

    return np.sqrt(np.sum(dx * dx, axis=1))


def get_radial_density_profile(radii, masses, bins: int) -> Tuple[np.ndarray]:
    """
    Gets the radial gas density profile, after generating similar bins to those
    used in similar works.
    """

    bins = np.logspace(-2, 1, bins)

    histogram, bin_edges = np.histogram(a=radii, weights=masses, bins=bins)

    volumes = np.array(
        [
            (4. * np.pi / 3.) * (r_outer ** 3 - r_inner ** 3)
            for r_outer, r_inner in zip(bin_edges[1:], bin_edges[:-1])
        ]
    )

    return histogram / volumes, bin_edges  # densities


def mu(T, H_frac, T_trans):
    """
    Get the molecular weight as a function of temperature.
    """
    if T > T_trans:
        return 4. / (8. - 5. * (1. - H_frac))
    else:
        return 4. / (1. + 3. * H_frac)


def T(u, metadata: MetaData):
    """
    Temperature of primordial gas.
    """

    gas_gamma = metadata.hydroscheme["Adiabatic index"][0]
    H_frac = metadata.hydroscheme["Hydrogen mass fraction"][0]
    T_trans = metadata.hydroscheme["Hydrogen ionization transition temperature"][0]

    k_in_J_K = 1.38064852e-23
    mH_in_kg = 1.6737236e-27

    T_over_mu = (gas_gamma - 1.) * u * mH_in_kg / k_in_J_K
    ret = np.ones(np.size(u)) * T_trans

    # Enough energy to be ionized?
    mask_ionized = T_over_mu > (T_trans + 1) / mu(T_trans + 1, H_frac, T_trans)
    if np.sum(mask_ionized) > 0:
        ret[mask_ionized] = T_over_mu[mask_ionized] * mu(T_trans * 10, H_frac, T_trans)

    # Neutral gas?
    mask_neutral = T_over_mu < (T_trans - 1) / mu((T_trans - 1), H_frac, T_trans)
    if np.sum(mask_neutral) > 0:
        ret[mask_neutral] = T_over_mu[mask_neutral] * mu(0, H_frac, T_trans)

    return ret


def get_radial_temperature_profile(
    data: SimulationParticleData, bins: int
) -> np.ndarray:
    """
    Gets the radial gas density profile, after generating similar bins to those
    used in similar works.
    """

    temperatures = T(data.gas.energies, data.metadata)
    radii = data.gas.radii

    bins = np.logspace(-2, 1, bins)

    histogram, _ = np.histogram(a=radii, weights=temperatures, bins=bins)

    counts, _ = np.histogram(a=radii, weights=np.ones_like(radii), bins=bins)

    return histogram / counts  # need to get mean value in bin


def get_radial_entropy_profile(data: SimulationParticleData, bins: int) -> np.ndarray:
    """
    Gets the radial gas density profile, after generating similar bins to those
    used in similar works.
    """

    gas_gamma = data.metadata.hydroscheme["Adiabatic index"][0]
    gamma_minus_one = gas_gamma - 1.0

    entropies = (
        data.gas.energies * (gamma_minus_one) / data.gas.densities ** gamma_minus_one
    )
    print("Warning: Current entropy profile assumes all gas is ionised")
    radii = data.gas.radii

    bins = np.logspace(-2, 1, bins)

    histogram, _ = np.histogram(a=radii, weights=entropies, bins=bins)

    counts, _ = np.histogram(a=radii, weights=np.ones_like(radii), bins=bins)

    return histogram / counts  # need to get mean value in bin


def nfw(R, halo_data: HaloData):
    """
    NFW profile at radius R.
    """

    R_s = halo_data.Rvir / halo_data.c
    rho_0 = (4 * np.pi * R_s ** 3) / (halo_data.Mvir)
    rho_0 *= np.log(1 + halo_data.c) - halo_data.c / (halo_data.c + 1)
    rho_0 = 1.0 / rho_0

    ratio = R / R_s

    return rho_0 / (ratio * (1 + ratio) ** 2)


def create_plot(
    data: SimulationParticleData,
    halo_data: HaloData,
    bins: int,
    create_images: bool,
    image_data: np.ndarray,
):
    """
    Creates the figure and axes objects and plots the data on them.
    """

    fig, axes = plt.subplots(2, 3, figsize=(12, 8))

    gas_density, bin_edges = get_radial_density_profile(
        data.gas.radii, data.gas.masses, bins=bins
    )
    dm_density, _ = get_radial_density_profile(
        data.dark_matter.radii, data.dark_matter.masses, bins=bins
    )
    temperature = get_radial_temperature_profile(data, bins=bins)
    entropy = get_radial_entropy_profile(data, bins=bins)

    bin_centers = [0.5 * (x + y) for x, y in zip(bin_edges[:-1], bin_edges[1:])]
    nfw_R = np.logspace(-2, 1, bins * 100)
    nfw_rho = nfw(nfw_R, halo_data)

    axes[0][0].loglog()
    axes[0][0].plot(nfw_R, 0.1 * nfw_rho, ls="dashed", color="grey")
    axes[0][0].scatter(bin_centers, gas_density)
    axes[0][0].set_ylabel(r"$\rho_{\rm gas} (R)$ [$10^{10}$ M$_\odot$ Mpc$^{-3}$]")
    axes[0][0].set_xlabel(r"R [Mpc]")
    axes[0][0].set_xlim(0.01, 10)

    axes[0][1].semilogx()
    axes[0][1].scatter(bin_centers, np.log(entropy))
    axes[0][1].set_ylabel(
        r"Entropy $\log(A$ [K ($10^{10}$ M$_\odot$)$^{2/3}$ Mpc$^{-2}$])"
    )
    axes[0][1].set_xlabel(r"R [Mpc]")
    axes[0][1].set_xlim(0.01, 10)

    if create_images:
        axes[0][2].imshow(np.log10(image_data))

    axes[0][2].set_xticks([])
    axes[0][2].set_yticks([])

    axes[1][0].loglog()
    axes[1][0].scatter(bin_centers, temperature)
    axes[1][0].set_ylabel(r"$T_{\rm gas} (R)$ [K]")
    axes[1][0].set_xlabel(r"R [Mpc]")
    axes[1][0].set_xlim(0.01, 10)

    axes[1][1].loglog()
    axes[1][1].scatter(bin_centers, dm_density)
    axes[1][1].plot(nfw_R, 0.9 * nfw_rho, ls="dashed", color="grey")
    axes[1][1].set_ylabel(r"$\rho_{\rm DM} (R)$ [$10^{10}$ M$_\odot$ Mpc$^{-3}$]")
    axes[1][1].set_xlabel(r"R [Mpc]")
    axes[1][1].set_xlim(0.01, 10)
    axes[1][1].text(
        0.02,
        5,
        "$c_{{vir}} = {:2.2f}$\n$R_{{vir}} = {:2.2f}$ Mpc\n$M_{{vir}} = {:2.2f}$ $10^{{10}}$ M$_\odot$".format(
            halo_data.c, halo_data.Rvir, halo_data.Mvir
        ),
        va="bottom",
        ha="left",
    )

    axes[1][2].text(
        -0.49,
        0.7,
        "Santa Barbara with $\\gamma={:2.2f}$ in 3D".format(
            data.metadata.hydroscheme["Adiabatic index"][0]
        ),
    )

    scheme_list = data.metadata.hydroscheme["Scheme"].decode("utf-8").split(" ")
    i = 4
    while i < len(scheme_list):
        scheme_list.insert(i, "\n")
        i += 4 + 1
    wrapped_scheme = " ".join(scheme_list)
    wrapped_scheme.replace("\n ", "\n")

    axes[1][2].text(-0.49, 0.8, wrapped_scheme)

    axes[1][2].plot([-0.49, 0.1], [0.62, 0.62], "k-", lw=1)

    axes[1][2].text(
        -0.49, 0.5, f"SWIFT {data.metadata.code['Git Revision'].decode('utf-8')}"
    )

    axes[1][2].text(
        -0.49,
        0.3,
        data.metadata.hydroscheme["Kernel function"].decode("utf-8"),
        fontsize=10,
    )
    axes[1][2].text(
        -0.49,
        0.2,
        "{:2.3f} neighbours ($\\eta={:3.3f}$)".format(
            data.metadata.hydroscheme["Kernel target N_ngb"][0],
            data.metadata.hydroscheme["Kernel eta"][0],
        ),
    )
    axes[1][2].set_xlim(-0.5, 0.5)
    axes[1][2].set_ylim(0, 1)
    axes[1][2].axis("off")

    fig.tight_layout()

    return fig, axes


if __name__ == "__main__":
    import sys

    filename = "santabarbara_{:04d}.hdf5".format(int(sys.argv[1]))
    catalogue_filename = f"{sys.argv[2]}/santabarbara.properties"

    try:
        bins = int(sys.argv[3])
    except:
        bins = 25

    halo_data = get_halo_data(catalogue_filename)
    simulation_data = load_data(filename, halo_data.center)

    if create_images:
        data = makeImage.read_data_from_file(filename, part_type=0)
        _, image_data = makeImage.generate_views(data)
        del data
    else:
        image_data = None

    fig, _ = create_plot(
        data=simulation_data,
        halo_data=halo_data,
        bins=bins,
        create_images=create_images,
        image_data=image_data,
    )

    fig.savefig("santabarbara.png", dpi=300)
