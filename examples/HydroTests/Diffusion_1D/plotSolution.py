"""
This file is part of SWIFT.
Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Plots the solution for the ContactDiscontinuity_1D test.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.integrate import solve_ivp
    solve_ode = True
except:
    solve_ode = False

from swiftsimio import load

matplotlib.use("Agg")

def solve_analytic(u_0, u_1, t_0, t_1, alpha=0.1):
    """
    Solves the analytic equation:

    $$
        \frac{d \Delta}{d t} = \kappa \alpha (
            \sqrt{u_0 + \Delta} + \sqrt{u_1 + \Delta}
        ) (
            u_1 - u_0 - 2 \Delta
        )
    $$

    from time t0 to t1.

    + alpha is the gradient term 
    + u_0 is the "low" state
    + u_1 is the "high" state.
    """

    if not solve_ode:
        return [0.0], [0.0]

    def gradient(t, u):
        """
        Returns du0/dt, du1/dt
        """

        common = alpha * (np.sqrt(u[0]) + np.sqrt(u[1])) * (u[0] - u[1])

        return np.array([-1.0 * common, 1.0 * common])

    ret = solve_ivp(gradient, t_span=[t_0.value, t_1.value], y0=[u_0.value, u_1.value], t_eval=np.linspace(t_0.value, t_1.value, 100))

    t = ret.t
    high = ret.y[1]
    low = ret.y[0]

    return t, (high - low) * u_0.units


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
        "SWIFT\n"
        + metadata.code_info
        + "\n\n"
        + "Compiler\n"
        + metadata.compiler_info
        + "\n\n"
        + "Hydrodynamics\n"
        + metadata.hydro_info
        + "\n\n"
        + "Viscosity\n"
        + viscosity
        + "\n\n"
        + "Diffusion\n"
        + diffusion
    )

    return output


def get_data_list(start: int, stop: int, handle: str):
    """
    Gets a list of swiftsimio objects that contains all of the data.
    """

    data = [load("{}_{:04d}.hdf5".format(handle, x)) for x in range(start, stop + 1)]

    return data


def setup_axes(size=[8, 8], dpi=300):
    """
    Sets up the axes with the correct labels, etc.
    """
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=size, dpi=dpi)

    ax = ax.flatten()

    ax[0].axis("off")
    ax[1].set_xlabel("Time [s]")
    ax[1].set_ylabel("Relative energy difference $u/\\left<u\\right>$")

    ax[2].set_xlabel("Time [s]")
    ax[2].set_ylabel("Deviation in position relative to MIPS [$\\times 10^{6}$]")

    ax[3].set_xlabel("Time [s]")

    return fig, ax

def mean_std_max_min(data):
    """
    Returns:
    mean, stdev, max, min
    for your data.
    """
    means = np.array([np.mean(x) for x in data])
    stdevs = np.array([np.std(x) for x in data])
    maxs = np.array([np.max(x) for x in data])
    mins = np.array([np.min(x) for x in data])

    return means, stdevs, maxs, mins


def extract_plottables_u(data_list):
    """
    Extracts the plottables for the internal energies. Returns:
    mean, stdev, max, min differences between adjacent internal energies
    """

    data = [
        np.diff(x.gas.internal_energy.value) / np.mean(x.gas.internal_energy.value)
        for x in data_list
    ]

    return mean_std_max_min(data)


def extract_plottables_x(data_list):
    """
    Extracts the plottables for positions. Returns:
    mean, stdev, max, min * 1e6 deviations from original position
    """

    n_part = data_list[0].metadata.n_gas
    boxsize = data_list[0].metadata.boxsize[0].value
    dx = boxsize / n_part

    original_x = np.arange(n_part, dtype=float) * dx + (0.5 * dx)
    
    deviations = [1e6 * abs(original_x - x.gas.coordinates.value[:, 0]) / dx for x in data_list]

    return mean_std_max_min(deviations)


def extract_plottables_rho(data_list):
    """
    Extracts the plottables for pressure. Returns:
    mean, stdev, max, min * 1e6 deviations from mean density 
    """

    P = [x.gas.density.value for x in data_list]
    mean_P = [np.mean(x) for x in P]
    deviations = [1e6 * (x - y) / x for x, y in zip(mean_P, P)]

    return mean_std_max_min(deviations)


def extract_plottables_diff(data_list):
    """
    Extracts the plottables for pressure. Returns:
    mean, stdev, max, min * 1e6 deviations from mean density 
    """

    P = [x.gas.diffusion.value for x in data_list]

    return mean_std_max_min(P)


def make_plot(start: int, stop: int, handle: str):
    """
    Makes the plot and returns the figure and axes objects.
    """
    fig, ax = setup_axes()
    data_list = get_data_list(start, stop, handle)
    data_dump = get_data_dump(data_list[0].metadata)
    t = [x.metadata.t for x in data_list]
    means, stdevs, maxs, mins = extract_plottables_u(data_list)
    x_means, x_stdevs, x_maxs, x_mins = extract_plottables_x(data_list)

    try:
        alpha = np.mean([np.mean(x.gas.diffusion) for x in data_list])
    except AttributeError:
        # Must be using a non-diffusive scheme.
        alpha = 0.0

    ax[0].text(
        0.5,
        0.5,
        data_dump,
        ha="center",
        va="center",
        fontsize=8,
        transform=ax[0].transAxes,
    )

    ax[1].fill_between(
        t, means - stdevs, means + stdevs, color="C0", alpha=0.5, edgecolor="none"
    )
    ax[1].plot(t, means, label="Mean", c="C0")
    ax[1].plot(t, maxs, label="Max", linestyle="dashed", c="C1")
    ax[1].plot(t, mins, label="Min", linestyle="dashed", c="C2")

    if solve_ode:
        times_ode, diff = solve_analytic(
            u_0=data_list[0].gas.internal_energy.min(),
            u_1=data_list[0].gas.internal_energy.max(),
            t_0=t[0],
            t_1=t[-1],
            alpha=(
                np.sqrt(5.0/3.0 * (5.0/3.0 - 1.0)) * 
                alpha / data_list[0].gas.smoothing_length[0].value
            )
        )

        ax[1].plot(
            times_ode,
            (diff) / np.mean(data_list[0].gas.internal_energy),
            label="Analytic",
            linestyle="dotted",
            c="C3"
        )

        #import pdb;pdb.set_trace()

    ax[2].fill_between(
        t, x_means - x_stdevs, x_means + x_stdevs, color="C0", alpha=0.5, edgecolor="none"
    )
    ax[2].plot(t, x_means, label="Mean", c="C0")
    ax[2].plot(t, x_maxs, label="Max", linestyle="dashed", c="C1")
    ax[2].plot(t, x_mins, label="Min", linestyle="dashed", c="C2")

    try:
        # Give diffusion info a go; this may not be present
        diff_means, diff_stdevs, diff_maxs, diff_mins = extract_plottables_diff(data_list)

        ax[3].set_ylabel(r"Diffusion parameter $\alpha_{diff}$")
        ax[3].fill_between(
            t, diff_means - diff_stdevs, diff_means + diff_stdevs, color="C0", alpha=0.5, edgecolor="none"
        )
        ax[3].plot(t, diff_means, label="Mean", c="C0")
        ax[3].plot(t, diff_maxs, label="Max", linestyle="dashed", c="C1")
        ax[3].plot(t, diff_mins, label="Min", linestyle="dashed", c="C2")

    except:
        # Diffusion info must not be present.
        rho_means, rho_stdevs, rho_maxs, rho_mins = extract_plottables_rho(data_list)

        ax[3].set_ylabel("Deviation from mean density $(\\rho_i - \\bar{\\rho}) / \\bar{\\rho}$ [$\\times 10^{6}$]")
        ax[3].fill_between(
            t, rho_means - rho_stdevs, rho_means + rho_stdevs, color="C0", alpha=0.5, edgecolor="none"
        )
        ax[3].plot(t, rho_means, label="Mean", c="C0")
        ax[3].plot(t, rho_maxs, label="Max", linestyle="dashed", c="C1")
        ax[3].plot(t, rho_mins, label="Min", linestyle="dashed", c="C2")

    ax[1].legend(loc=1, markerfirst=False)

    ax[1].set_xlim(t[0], t[-1])
    ax[2].set_xlim(t[0], t[-1])
    ax[3].set_xlim(t[0], t[-1])

    fig.tight_layout()

    return fig, ax


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(
        description="Makes a plot of the data from the ContactDiscontinuity_1D test."
    )

    parser.add_argument(
        "-i", "--initial", help="Initial snapshot. Default: 0", type=int, default=0
    )

    parser.add_argument(
        "-f", "--final", help="Final snapshot. Default: 50", type=int, default=50
    )

    parser.add_argument(
        "-s",
        "--snapshot",
        help="First part of the snapshot filename. Default: diffusion",
        type=str,
        default="diffusion",
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output filename. Default: diffusion.png",
        type=str,
        default="diffusion.png",
    )

    args = vars(parser.parse_args())

    fig, ax = make_plot(args["initial"], args["final"], args["snapshot"])

    fig.savefig(args["output"])
