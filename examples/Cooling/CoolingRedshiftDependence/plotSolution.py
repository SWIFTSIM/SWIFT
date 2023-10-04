"""
Plots the mean temperature.
"""

import matplotlib.pyplot as plt

from swiftsimio import load

from unyt import Myr, K, mh, cm, erg

import numpy as np


def get_data_dump(metadata):
    """
    Gets a big data dump from the SWIFT metadata
    """

    try:
        viscosity = metadata.viscosity_info
    except:
        viscosity = "No info"

    try:
        diffusion = metadata.diffusion_info
    except:
        diffusion = "No info"

    output = (
        "$\\bf{SWIFT}$\n"
        + metadata.code_info
        + "\n\n"
        + "$\\bf{Compiler}$\n"
        + metadata.compiler_info
        + "\n\n"
        + "$\\bf{Hydrodynamics}$\n"
        + metadata.hydro_info
        + "\n\n"
        + "$\\bf{Cooling}$\n"
        + metadata.subgrid_scheme["Cooling Model"].decode("utf-8")
        + "\n\n"
        + "$\\bf{Viscosity}$\n"
        + viscosity
        + "\n\n"
        + "$\\bf{Diffusion}$\n"
        + diffusion
    )

    return output


def setup_axes():
    """
    Sets up the axes and returns fig, ax
    """

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))

    ax = ax.flatten()

    ax[0].semilogy()
    ax[2].semilogy()

    ax[1].set_xlabel("Simulation time elapsed [Myr]")
    ax[2].set_xlabel("Simulation time elapsed [Myr]")

    ax[0].set_ylabel("Temperature of Universe [K]")
    ax[1].set_ylabel("Physical Density of Universe [$n_H$ cm$^{-3}$]")
    ax[2].set_ylabel("Physical Energy [erg]")

    for a in ax[:-1]:
        a.set_xlim(0, 100)

    return fig, ax


def get_data(handle: float, n_snaps: int):
    """
    Returns the elapsed simulation time, temperature, and density
    """

    t0 = 0.0

    times = []
    temps = []
    densities = []
    energies = []
    radiated_energies = []

    for snap in range(n_snaps):
        try:
            data = load(f"data/{handle}_{snap:04d}.hdf5")

            if snap == 0:
                t0 = data.metadata.t.to(Myr).value

            times.append(data.metadata.t.to(Myr).value - t0)
            temps.append(np.mean(data.gas.temperatures.to(K).value))
            densities.append(
                np.mean(data.gas.densities.to(mh / cm ** 3).value)
                / (data.metadata.scale_factor ** 3)
            )

            try:
                energies.append(
                    np.mean(
                        (data.gas.internal_energies * data.gas.masses)
                        .astype(np.float64)
                        .to(erg)
                        .value
                    )
                    * data.metadata.scale_factor ** (2)
                )
                radiated_energies.append(
                    np.mean(data.gas.radiated_energy.astype(np.float64).to(erg).value)
                )
            except AttributeError:
                continue
        except OSError:
            continue

    return times, temps, densities, energies, radiated_energies


def get_n_steps(timesteps_filename: str) -> int:
    """
    Reads the number of steps from the timesteps file.
    """

    data = np.genfromtxt(timesteps_filename)

    return int(data[-1][0])


def plot_single_data(
    handle: str, n_snaps: int, label: str, ax: plt.axes, n_steps: int = 0, run: int = 0
):
    """
    Takes the a single simulation and plots it on the axes.
    """

    data = get_data(handle, n_snaps)

    ax[0].plot(data[0], data[1], label=label, marker="o", ms=2, c=f"C{run}")

    ax[1].plot(
        data[0], data[2], label=f"Steps: {n_steps}", marker="o", ms=2, c=f"C{run}"
    )

    if run == 0:
        label_energy = "Particle Energy"
        label_radiated = "Radiated Energy"
        label_sum = "Sum"
    else:
        label_energy = ""
        label_radiated = ""
        label_sum = ""

    ax[2].plot(data[0], data[3], label=label_energy, ls="dotted", color=f"C{run}")

    # ax[2].plot(data[0], data[4], label=label_radiated, ls="dashed", C=f"C{run}")

    # ax[2].plot(
    #    data[0], [x + y for x, y in zip(*data[3:5])], label=label_sum, C=f"C{run}"
    # )

    return


def make_plot(handles, names, timestep_names, n_snaps=100):
    """
    Makes the whole plot and returns the fig, ax objects.
    """

    fig, ax = setup_axes()

    run = 0

    for handle, name, timestep_name in zip(handles, names, timestep_names):
        try:
            n_steps = get_n_steps(timestep_name)
        except:
            n_steps = "Unknown"

        plot_single_data(handle, n_snaps, name, ax, n_steps=n_steps, run=run)
        run += 1

    ax[0].legend()
    ax[1].legend()
    ax[2].legend()

    info_axis = ax[-1]

    try:
        info = get_data_dump(load(f"data/{handles[0]}_0000.hdf5").metadata)

        info_axis.text(
            0.5,
            0.45,
            info,
            ha="center",
            va="center",
            fontsize=7,
            transform=info_axis.transAxes,
        )
    except OSError:
        pass

    info_axis.axis("off")

    return fig, ax


if __name__ == "__main__":
    """
    Plot everything!
    """

    handles = [
        "redshift_dependence_no_z",
        "redshift_dependence_low_z",
        "redshift_dependence_high_z",
    ]
    names = ["No Cosmology", "Low Redshift ($z=0.01$)", "High Redshift ($z=1$)"]
    timestep_names = ["timesteps_no.txt", "timesteps_low.txt", "timesteps_high.txt"]

    fig, ax = make_plot(handles, names, timestep_names)

    fig.tight_layout(pad=0.5)
    fig.savefig("redshift_dependence.png", dpi=300)
