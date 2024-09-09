from swiftsimio import load
from swiftsimio.visualisation.volume_render import render_gas
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from unyt import s, kb, mp
import glob


def analyse_snapshot(i, chi, rho_cl0, m_cl0, resolution):
    print("analysing snapshot", i)
    try:
        data = load("./{}/blob_{:0=4d}.hdf5".format(resolution, i))
    except OSError:
        print("snapshot {} is not found".format(i))
        exit()

    cloud = data.gas.densities.value > rho_cl0 / 3
    m_cl = np.sum(data.gas.masses.value[cloud])
    dense_mass = m_cl / m_cl0

    # Resample onto a 16^3 grid
    # Naive histogram
    x_cl = np.median(data.gas.coordinates[:, 0].value[cloud])
    bins_y = np.linspace(0, 1, 17)
    bins_z = np.linspace(0, 1, 17)
    bins_x = np.linspace(0, 4, 4 * 16 + 1)
    # Shift x_bins such that x_cl is in the centre of a bin
    bin_space_x = 1 / 16
    bins_x += x_cl - bin_space_x / 2
    bins_x = np.sort(bins_x % 4)
    print(bins_x)
    binned_masses, _ = np.histogramdd(
        data.gas.coordinates.v, weights=data.gas.masses.v, bins=[bins_x, bins_y, bins_z]
    )
    binned_masses = binned_masses * data.gas.masses.units
    densities_resampled = binned_masses / (data.metadata.boxsize[1] / 16) ** 3
    densities_resampled = densities_resampled.to(data.gas.densities.units)

    cloud_resampled = densities_resampled.value > rho_cl0 / 3
    m_cl_resampled = np.sum(binned_masses.value[cloud_resampled])
    # For volume rendering we compute the mass from the density and the volume
    dense_mass_resampled_16 = m_cl_resampled / m_cl0

    # Resample onto a 32^3 grid
    # Naive histogram
    x_cl = np.median(data.gas.coordinates[:, 0].value[cloud])
    bins_y = np.linspace(0, 1, 33)
    bins_z = np.linspace(0, 1, 33)
    bins_x = np.linspace(0, 4, 4 * 32 + 1)
    # Shift x_bins such that x_cl is in the centre of a bin
    bin_space_x = 1 / 32
    bins_x += x_cl - bin_space_x / 2
    bins_x = np.sort(bins_x % 4)
    print(bins_x)
    binned_masses, _ = np.histogramdd(
        data.gas.coordinates.v, weights=data.gas.masses.v, bins=[bins_x, bins_y, bins_z]
    )
    binned_masses = binned_masses * data.gas.masses.units
    densities_resampled = binned_masses / (data.metadata.boxsize[1] / 32) ** 3
    densities_resampled = densities_resampled.to(data.gas.densities.units)

    cloud_resampled = densities_resampled.value > rho_cl0 / 3
    m_cl_resampled = np.sum(binned_masses.value[cloud_resampled])
    # For volume rendering we compute the mass from the density and the volume
    dense_mass_resampled_32 = m_cl_resampled / m_cl0

    t = data.metadata.time.to(s).value

    # Convert internal energies to temperatures
    temperatures = (((5 / 3 - 1) * data.gas.internal_energies / kb) * mp).value
    data.gas.temperatures = temperatures

    if chi == 100:
        m_t_intermediate = (
            np.sum(
                data.gas.masses.value[
                    (data.gas.temperatures > 10 ** 4.5)
                    & (data.gas.temperatures < 10 ** 5.5)
                ]
            )
            / m_cl0
        )
    elif chi == 10:
        m_t_intermediate = (
            np.sum(
                data.gas.masses.value[
                    (data.gas.temperatures > 10 ** 4.25)
                    & (data.gas.temperatures < 10 ** 4.75)
                ]
            )
            / m_cl0
        )

    """
    resample the intermediate temperature gas to 32^3 resolution
    """
    x_cl = np.median(data.gas.coordinates[:, 0].value[cloud])
    bins_y = np.linspace(0, 1, 33)
    bins_z = np.linspace(0, 1, 33)
    bins_x = np.linspace(0, 4, 4 * 32 + 1)
    # Shift x_bins such that x_cl is in the centre of a bin
    bin_space_x = 1 / 32
    bins_x += x_cl - bin_space_x / 2
    bins_x = np.sort(bins_x % 4)
    print(bins_x)
    binned_masses, _ = np.histogramdd(
        data.gas.coordinates.v, weights=data.gas.masses.v, bins=[bins_x, bins_y, bins_z]
    )
    binned_masstemp, _ = np.histogramdd(
        data.gas.coordinates.v,
        weights=data.gas.masses.v * data.gas.temperatures,
        bins=[bins_x, bins_y, bins_z],
    )
    binned_temperatures = binned_masstemp / binned_masses
    binned_masses = binned_masses * data.gas.masses.units

    if chi == 100:
        mask = (binned_temperatures > 10 ** 4.5) & (binned_temperatures < 10 ** 5.5)
        m_t_intermediate_resampled_32 = np.sum(binned_masses.v[mask]) / m_cl0
    elif chi == 10:
        mask = (binned_temperatures > 10 ** 4.25) & (binned_temperatures < 10 ** 4.75)
        m_t_intermediate_resampled_32 = np.sum(binned_masses.v[mask]) / m_cl0

    """
    resample the intermediate temperature gas to 16^3 resolution 
    """
    x_cl = np.median(data.gas.coordinates[:, 0].value[cloud])
    bins_y = np.linspace(0, 1, 17)
    bins_z = np.linspace(0, 1, 17)
    bins_x = np.linspace(0, 4, 4 * 16 + 1)
    # Shift x_bins such that x_cl is in the centre of a bin
    bin_space_x = 1 / 16
    bins_x += x_cl - bin_space_x / 2
    bins_x = np.sort(bins_x % 4)
    print(bins_x)
    binned_masses, _ = np.histogramdd(
        data.gas.coordinates.v, weights=data.gas.masses.v, bins=[bins_x, bins_y, bins_z]
    )
    binned_masstemp, _ = np.histogramdd(
        data.gas.coordinates.v,
        weights=data.gas.masses.v * data.gas.temperatures,
        bins=[bins_x, bins_y, bins_z],
    )
    binned_temperatures = binned_masstemp / binned_masses
    binned_masses = binned_masses * data.gas.masses.units

    if chi == 100:
        mask = (binned_temperatures > 10 ** 4.5) & (binned_temperatures < 10 ** 5.5)
        m_t_intermediate_resampled_16 = np.sum(binned_masses.v[mask]) / m_cl0
    elif chi == 10:
        mask = (binned_temperatures > 10 ** 4.25) & (binned_temperatures < 10 ** 4.75)
        m_t_intermediate_resampled_16 = np.sum(binned_masses.v[mask]) / m_cl0

    return (
        t,
        dense_mass,
        m_t_intermediate,
        dense_mass_resampled_16,
        m_t_intermediate_resampled_16,
        dense_mass_resampled_32,
        m_t_intermediate_resampled_32,
    )


def load_first_snapshot(chi, resolution):
    # Load snapshot 0
    try:
        data_0 = load("./{}/blob_0000.hdf5".format(resolution))
    except OSError:
        print("snapshot 0 is not found")

    v_wind0 = data_0.gas.velocities[
        0, 0
    ]  # All wind particles have the same initial velocity
    r_cl0 = data_0.metadata.boxsize[1] * 0.1  # The cloud is 0.1 * boxsize

    rho_cl0 = chi * np.percentile(data_0.gas.densities.value, 50)
    cloud_0 = (
        data_0.gas.densities.value > rho_cl0 / 3
    )  # Create a mask for the dense gas
    m_cl0 = np.sum(data_0.gas.masses.value[cloud_0])

    t_cc = (np.sqrt(chi) * r_cl0 / v_wind0).to(s)  # Cloud crushing time

    # This will give the number of snapshots + 1 IC file
    hdf5counter = len(glob.glob1("./{}".format(resolution), "*.hdf5"))
    # Subtract the IC file from the count
    steps = hdf5counter - 1

    return rho_cl0, m_cl0, steps, t_cc


# What is the number of (particles)^1/3?
resolution = 128
# What is the density contrast?
chi = 10

rho_cl0, m_cl0, steps, t_cc = load_first_snapshot(chi, resolution)

time = np.zeros(steps)
dense_mass = np.zeros(steps)
intermediate_temp_mass = np.zeros(steps)

dense_mass_16 = np.zeros(steps)
intermediate_temp_mass_16 = np.zeros(steps)

dense_mass_32 = np.zeros(steps)
intermediate_temp_mass_32 = np.zeros(steps)

for i in range(steps):
    time[i], dense_mass[i], intermediate_temp_mass[i], dense_mass_16[
        i
    ], intermediate_temp_mass_16[i], dense_mass_32[i], intermediate_temp_mass_32[
        i
    ] = analyse_snapshot(
        i, chi, rho_cl0, m_cl0, resolution
    )

fig, ax = plt.subplots(3, 2, figsize=(10, 15))

ax[0, 0].plot(time / t_cc, dense_mass)
ax[0, 0].set_title("Dense Mass")

ax[0, 1].plot(time / t_cc, intermediate_temp_mass)
ax[0, 1].set_title("Intermediate-temperature Mass")

ax[1, 0].plot(time / t_cc, dense_mass_16)
ax[1, 0].set_title("Dense Mass resampled 16^3")

ax[1, 1].plot(time / t_cc, intermediate_temp_mass_16)
ax[1, 1].set_title("Intermediate-temperature Mass resampled 16^3")

ax[2, 0].plot(time / t_cc, dense_mass_32)
ax[2, 0].set_title("Dense Mass resampled 32^3")

ax[2, 1].plot(time / t_cc, intermediate_temp_mass_32)
ax[2, 1].set_title("Intermediate-temperature Mass resampled 32^3")

plt.savefig("cloud_evolution")
plt.close()
