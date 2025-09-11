###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
#                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################
import argparse
import sys

import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import stats
from scipy.special import gamma as Gamma
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import astropy.units as units


def calc_a(g, nu=3):
    """
    Calculate the exponents of the polynomials for the Sedov solution.

    Parameters
    ----------
    g : float
        The polytropic gas gamma.
    nu : int, optional
        The dimension of the explosion. The default is 3.

    Returns
    -------
    list
        A list of 8 coefficients for the Sedov solution.

    """
    a = [0] * 8

    a[0] = 2.0 / (nu + 2)
    a[2] = (1 - g) / (2 * (g - 1) + nu)
    a[3] = nu / (2 * (g - 1) + nu)
    a[5] = 2 / (g - 2)
    a[6] = g / (2 * (g - 1) + nu)

    a[1] = (((nu + 2) * g) / (2.0 + nu * (g - 1.0))) * (
        (2.0 * nu * (2.0 - g)) / (g * (nu + 2.0) ** 2) - a[2]
    )
    a[4] = a[1] * (nu + 2) / (2 - g)
    a[7] = (2 + nu * (g - 1)) * a[1] / (nu * (2 - g))
    return a


def calc_beta(v, g, nu=3):
    """
    Calculate the coefficients of the polynomials for the similarity variables.

    Parameters
    ----------
    v : numpy.ndarray
        The similarity variable array.
    g : float
        The polytropic gas gamma.
    nu : int, optional
        The dimension of the explosion. The default is 3.

    Returns
    -------
    numpy.ndarray
        A 4xN array of beta values for the Sedov solution.

    """

    beta = (
        (nu + 2)
        * (g + 1)
        * np.array(
            (
                0.25,
                (g / (g - 1)) * 0.5,
                -(2 + nu * (g - 1))
                / 2.0
                / ((nu + 2) * (g + 1) - 2 * (2 + nu * (g - 1))),
                -0.5 / (g - 1),
            ),
            dtype=np.float64,
        )
    )

    beta = np.outer(beta, v)

    beta += (g + 1) * np.array(
        (
            0.0,
            -1.0 / (g - 1),
            (nu + 2) / ((nu + 2) * (g + 1) - 2.0 * (2 + nu * (g - 1))),
            1.0 / (g - 1),
        ),
        dtype=np.float64,
    ).reshape((4, 1))

    return beta


def sedov(t, E0, rho0, g, n=1000, nu=3):
    """
    Solve the standard Sedov blast wave problem.

    This function computes the analytical solution for a blast wave in a uniform
    medium, returning the profiles for radius, pressure, density, and velocity.

    Parameters
    ----------
    t : float
        The time of the blast wave.
    E0 : float
        The initial explosion energy.
    rho0 : float
        The initial background density.
    g : float
        The polytropic gas gamma.
    n : int, optional
        The number of points in the similarity variable array. The default is 1000.
    nu : int, optional
        The dimension of the explosion. The default is 3.

    Returns
    -------
    tuple
        A tuple containing:
        - r (numpy.ndarray): Radial position profile.
        - p (numpy.ndarray): Pressure profile.
        - rho (numpy.ndarray): Density profile.
        - u (numpy.ndarray): Specific internal energy profile.
        - v (numpy.ndarray): Velocity profile.
        - r_s (float): Shock radius.
        - p_s (float): Post-shock pressure.
        - rho_s (float): Post-shock density.
        - u_s (float): Post-shock specific internal energy.
        - shock_speed (float): The speed of the shock front.
    """
    # the similarity variable
    v_min = 2.0 / ((nu + 2) * g)
    v_max = 4.0 / ((nu + 2) * (g + 1))

    v = v_min + np.arange(n) * (v_max - v_min) / (n - 1.0)

    a = calc_a(g, nu)
    beta = calc_beta(v, g=g, nu=nu)
    lbeta = np.log(beta)

    r = np.exp(-a[0] * lbeta[0] - a[2] * lbeta[1] - a[1] * lbeta[2])
    rho = ((g + 1.0) / (g - 1.0)) * np.exp(
        a[3] * lbeta[1] + a[5] * lbeta[3] + a[4] * lbeta[2]
    )
    p = np.exp(
        nu * a[0] * lbeta[0] + (a[5] + 1) * lbeta[3] + (a[4] - 2 * a[1]) * lbeta[2]
    )
    u = beta[0] * r * 4.0 / ((g + 1) * (nu + 2))
    p *= 8.0 / ((g + 1) * (nu + 2) * (nu + 2))

    # we have to take extra care at v=v_min, since this can be a special point.
    # It is not a singularity, however, the gradients of our variables (wrt v) are.
    # r -> 0, u -> 0, rho -> 0, p-> constant

    u[0] = 0.0
    rho[0] = 0.0
    r[0] = 0.0
    p[0] = p[1]

    # volume of an n-sphere
    vol = (np.pi ** (nu / 2.0) / Gamma(nu / 2.0 + 1)) * np.power(r, nu)

    # note we choose to evaluate the integral in this way because the
    # volumes of the first few elements (i.e near v=vmin) are shrinking
    # very slowly, so we dramatically improve the error convergence by
    # finding the volumes exactly. This is most important for the
    # pressure integral, as this is on the order of the volume.

    # (dimensionless) energy of the model solution
    de = rho * u * u * 0.5 + p / (g - 1)
    # integrate (trapezium rule)
    q = np.inner(de[1:] + de[:-1], np.diff(vol)) * 0.5

    # the factor to convert to this particular problem
    fac = (q * (t**nu) * rho0 / E0) ** (-1.0 / (nu + 2))

    # shock speed
    shock_speed = fac * (2.0 / (nu + 2))
    rho_s = ((g + 1) / (g - 1)) * rho0
    r_s = shock_speed * t * (nu + 2) / 2.0
    p_s = (2.0 * rho0 * shock_speed * shock_speed) / (g + 1)
    u_s = (2.0 * shock_speed) / (g + 1)

    r *= fac * t
    u *= fac
    p *= fac * fac * rho0
    rho *= rho0
    return r, p, rho, u, v, r_s, p_s, rho_s, u_s, shock_speed


def sedov_solution(t, E0, rho0, P0, g, n=1000, nu=3):
    """
    Solve the Sedov blast wave problem and extend the solution past the shock front.

    This function calls the standard Sedov solver and then appends ambient
    conditions to the solution arrays to represent the region outside the
    blast wave.

    Parameters
    ----------
    t : float
        The time of the blast wave.
    E0 : float
        The initial explosion energy.
    rho0 : float
        The initial background density.
    P0 : float
        The initial background pressure.
    g : float
        The polytropic gas gamma.
    n : int, optional
        The number of points for the standard Sedov solution. The default is 1000.
    nu : int, optional
        The dimension of the explosion. The default is 3.

    Returns
    -------
    tuple
        A tuple containing the extended analytical profiles and post-shock values:
        - r (numpy.ndarray): Extended radial position profile.
        - p (numpy.ndarray): Extended pressure profile.
        - rho (numpy.ndarray): Extended density profile.
        - u (numpy.ndarray): Extended specific internal energy profile.
        - r_s (float): Shock radius from the standard solution.
        - p_s (float): Post-shock pressure from the standard solution.
        - rho_s (float): Post-shock density from the standard solution.
        - u_s (float): Post-shock specific internal energy from the standard solution.
        - shock_speed (float): The speed of the shock front.

    """
    n_2 = n - 2

    r, p, rho, u, v, r_s, p_s, rho_s, u_s, shock_speed = sedov(t, E0, rho0, g, n_2, nu)

    # Append points for after the shock
    r_shock = np.asarray(r_s).flatten()
    r = np.insert(r, np.size(r), [r_shock[0], r_shock[0] * 1.5])
    rho = np.insert(rho, np.size(rho), [rho0, rho0])
    p = np.insert(p, np.size(p), [P0, P0])
    v = np.insert(v, np.size(v), [0, 0])

    u_0 = P0 / (rho0 * (g - 1.0))
    u = np.insert(u, np.size(u), [u_0, u_0])
    return r, p, rho, u, r_s, p_s, rho_s, u_s, shock_speed


def sedov_fit_model(r, E_0, rho_0, P_0, time, g=1.4, nu=3):
    """
    Return interpolated values for density, pressure, and velocity from a Sedov model.

    This function acts as a wrapper for the Sedov solution, making it suitable
    for use with `scipy.optimize.curve_fit`. It solves the Sedov problem for
    the given parameters and returns interpolated values at the input radii.

    Parameters
    ----------
    r_fit : numpy.ndarray
        The array of radii at which to evaluate the model.
    E_0 : float
        The fitted explosion energy.
    rho_0 : float
        The fitted background density.
    P_0 : float
        The fitted background pressure.
    time : float
        The time of the blast wave.
    g : float, optional
        The polytropic gas gamma. The default is 1.4.
    nu : int, optional
        The dimension of the explosion. The default is 3.

    Returns
    -------
    numpy.ndarray
        A flattened array containing the interpolated density, pressure, and
        velocity values.

    """
    # Solve Sedov blast wave problem
    r_s, P_s, rho_s, v_s, r_shock, _, _, _, _ = sedov_solution(
        t=time, E0=E_0, rho0=rho_0, P0=P_0, g=g, n=len(r), nu=nu
    )

    # Ensure monotonicity for interpolation (sort values if needed)
    sorted_indices = np.argsort(r_s)
    r_s = r_s[sorted_indices]
    P_s = P_s[sorted_indices]
    rho_s = rho_s[sorted_indices]
    v_s = v_s[sorted_indices]

    # Interpolation functions
    interp_rho = interp1d(r_s, rho_s, kind="linear", fill_value="extrapolate")
    interp_P = interp1d(r_s, P_s, kind="linear", fill_value="extrapolate")
    interp_v = interp1d(r_s, v_s, kind="linear", fill_value="extrapolate")

    # Interpolated values at given r
    rho_fit = interp_rho(r)
    P_fit = interp_P(r)
    v_fit = interp_v(r)

    # Assign ambient conditions for r > r_shock
    mask_post_shock = r > r_shock
    rho_fit[mask_post_shock] = rho_0
    P_fit[mask_post_shock] = P_0
    v_fit[mask_post_shock] = P_0 / (rho_0 * (g - 1.0))  # v_0

    # Flatten and return as a single vector
    return np.hstack((rho_fit, P_fit, v_fit))


def cost_function(log_params, r_bin, data, sigma_data, time, gas_gamma):
    """
    Calculate the sum of squared residuals for the Sedov fit in log-space.
    """
    # Unpack the log-transformed parameters
    log_E_0, log_rho_0, log_P_0 = log_params

    # Convert back to linear space
    E_0 = np.exp(log_E_0)
    rho_0 = np.exp(log_rho_0)
    P_0 = np.exp(log_P_0)

    # Calculate the model values
    model_data = sedov_fit_model(r_bin, E_0, rho_0, P_0, time, gas_gamma)

    # Calculate the weighted residuals, handling potential NaNs
    residuals = (data - model_data) / sigma_data
    residuals[np.isnan(residuals)] = 1e10

    return np.sum(residuals**2)


class ExplicitDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Formatter that prints default arguments if provided."""

    def _get_help_string(self, action):
        if isinstance(action.default, np.ndarray) or action.default not in (
            None,
            False,
        ):
            return super()._get_help_string(action)
        return action.help


class RawTextArgumentDefaultsHelpFormatter(
    ExplicitDefaultsHelpFormatter, argparse.RawTextHelpFormatter
):
    """Combines raw text formatting with default argument printing."""

    pass


def parse_options():
    """
    Parse command-line arguments for the Sedov blast wave analysis script.

    Returns
    -------
    argparse.Namespace
        An object containing the parsed command-line arguments.

    """
    parser = argparse.ArgumentParser(
        description="Compute and analyze the Sedov blast wave from a SWIFT simulation.",
        epilog="Example usage:\npython script.py output.hdf5 --rho_0 1.0 --E_0 1e51",
        formatter_class=RawTextArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("file", type=str, help="Path to the SWIFT output file.")

    parser.add_argument(
        "--rho_0",
        type=float,
        default=None,
        help="Background density (in simulation units).",
    )
    parser.add_argument(
        "--P_0",
        type=float,
        default=None,
        help="Background pressure (in simulation units).",
    )
    parser.add_argument(
        "--E_0",
        type=float,
        default=None,
        help="Explosion energy (in simulation units).",
    )

    parser.add_argument(
        "-o",
        type=str,
        dest="output_filename",
        default=None,
        help="Output filename for the plot. If not set, the plot is shown instead.",
    )

    parser.add_argument(
        "--image_format",
        type=str,
        default="png",
        help="Format of the output image (recognized by Matplotlib).",
    )
    parser.add_argument(
        "--dpi", type=float, default=300, help="Resolution of output images."
    )
    parser.add_argument("--step", type=float, default=0.01, help="Step for binning.")
    parser.add_argument(
        "--stylesheet", type=str, default=None, help="Matplotlib stylesheet."
    )

    return parser.parse_args()


# Parse command-line arguments
opt = parse_options()

# Apply Matplotlib stylesheet
if opt.stylesheet:
    plt.style.use(opt.stylesheet)

# Physical parameters
rho_0 = opt.rho_0  # Background Density
P_0 = opt.P_0  # Background Pressure
E_0 = opt.E_0  # Energy of the explosion
step = opt.step  # Binning step size

gas_gamma = 5.0 / 3.0  # Gas polytropic index
n_shock = 1.5  # Shock index

# Ensure the file exists
file = opt.file

output_filename = opt.output_filename
image_format = opt.image_format
dpi = opt.dpi

try:
    with h5py.File(file, "r") as sim:
        # Read unit conversion factors from the snapshot
        try:
            U_M = sim["/InternalCodeUnits"].attrs.get("Unit mass in cgs (U_M)", 1.0)[0]
            U_L = sim["/InternalCodeUnits"].attrs.get("Unit length in cgs (U_L)", 1.0)[
                0
            ]
            U_t = sim["/InternalCodeUnits"].attrs.get("Unit time in cgs (U_t)", 1.0)[0]

        except KeyError:
            print("Warning: /Units group not found. Assuming CGS units.")
            pass  # Continue with default CGS-like conversion factors

        # Read global simulation metadata
        box_size = sim["/Header"].attrs.get("BoxSize", [np.nan])[0]
        time = sim["/Header"].attrs.get("Time", np.nan)
        scheme = sim["/HydroScheme"].attrs.get("Scheme", "Unknown")
        kernel = sim["/HydroScheme"].attrs.get("Kernel function", "Unknown")
        n_ngb = sim["/HydroScheme"].attrs.get("Kernel target N_ngb", np.nan)
        eta = sim["/HydroScheme"].attrs.get("Kernel eta", np.nan)
        git_version = sim["Code"].attrs.get("Git Revision", "Unknown")

        # Read gas particle data
        pos = sim["/PartType0/Coordinates"][:]
        vel = sim["/PartType0/Velocities"][:]
        u = sim["/PartType0/InternalEnergies"][:]
        S = sim["/PartType0/Entropies"][:]
        P = sim["/PartType0/Pressures"][:]
        rho = sim["/PartType0/Densities"][:]
        mass = sim["/PartType0/Masses"][:]

        # Compute radial positions centered on explosion
        x, y, z = (
            pos[:, 0] - box_size / 2,
            pos[:, 1] - box_size / 2,
            pos[:, 2] - box_size / 2,
        )
        r = np.sqrt(x**2 + y**2 + z**2)

        # Compute radial velocity component
        v_r = (x * vel[:, 0] + y * vel[:, 1] + z * vel[:, 2]) / r

        # Convert all data to cgs in place
        box_size *= U_L
        time *= U_t
        pos *= U_L
        r *= U_L
        vel *= U_L / U_t
        v_r *= U_L / U_t
        u *= U_L**2 / U_t**2
        P *= U_M / (U_L * U_t**2)
        rho *= U_M / (U_L**3)
        mass *= U_M

        # We also need to convert the initial guess values from the command
        # line
        step = step * U_L
        if rho_0 is not None:
            rho_0 = rho_0 * U_M / (U_L**3)
        if P_0 is not None:
            P_0 = P_0 * U_M / (U_L * U_t**2)
        if E_0 is not None:
            E_0 = E_0 * U_M * U_L**2 / U_t**2

except FileNotFoundError:
    sys.exit(f"Error: File '{file}' not found.")
except KeyError as e:
    sys.exit(f"Error: Missing expected dataset {e} in file '{file}'.")
except Exception as e:
    sys.exit(f"Unexpected error: {e}")

# Bin the data
r_bin_edge = np.arange(0.0, r.max(), step)
r_bin = 0.5 * (r_bin_edge[1:] + r_bin_edge[:-1])
rho_bin, _, _ = stats.binned_statistic(r, rho, statistic="mean", bins=r_bin_edge)
v_bin, _, _ = stats.binned_statistic(r, v_r, statistic="mean", bins=r_bin_edge)
P_bin, _, _ = stats.binned_statistic(r, P, statistic="mean", bins=r_bin_edge)
S_bin, _, _ = stats.binned_statistic(r, S, statistic="mean", bins=r_bin_edge)
u_bin, _, _ = stats.binned_statistic(r, u, statistic="mean", bins=r_bin_edge)
rho2_bin, _, _ = stats.binned_statistic(r, rho**2, statistic="mean", bins=r_bin_edge)
v2_bin, _, _ = stats.binned_statistic(r, v_r**2, statistic="mean", bins=r_bin_edge)
P2_bin, _, _ = stats.binned_statistic(r, P**2, statistic="mean", bins=r_bin_edge)
S2_bin, _, _ = stats.binned_statistic(r, S**2, statistic="mean", bins=r_bin_edge)
u2_bin, _, _ = stats.binned_statistic(r, u**2, statistic="mean", bins=r_bin_edge)

# Remove non physical values
valid_mask = np.isfinite(rho_bin) & np.isfinite(P_bin) & np.isfinite(v_bin)
r_bin = r_bin[valid_mask]
rho_bin = rho_bin[valid_mask]
v_bin = v_bin[valid_mask]
P_bin = P_bin[valid_mask]
S_bin = S_bin[valid_mask]
u_bin = u_bin[valid_mask]
rho2_bin = rho2_bin[valid_mask]
v2_bin = v2_bin[valid_mask]
P2_bin = P2_bin[valid_mask]
S2_bin = S2_bin[valid_mask]
u2_bin = u2_bin[valid_mask]

# Compute the error bars
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin**2)
v_sigma_bin = np.sqrt(v2_bin - v_bin**2)
P_sigma_bin = np.sqrt(P2_bin - P_bin**2)
u_sigma_bin = np.sqrt(u2_bin - u_bin**2)

################################
# Now, work out the solution....
################################

if rho_0 is None:
    rho_0_guess = np.median(rho)
else:
    rho_0_guess = rho_0

if P_0 is None:
    P_0_guess = np.median(P)
else:
    P_0_guess = P_0

if E_0 is None:
    # For E_0, find the point of maximum change in density to estimate the shock
    # radius.
    d_rho_bin = np.gradient(rho_bin)
    r_shock_guess = r_bin[np.argmax(d_rho_bin)]

    # Guess E_0 using a known analytical relationship for a Sedov explosion.
    E_0_guess = rho_0_guess * (r_shock_guess**5) / (time**2)
    E_0_guess = E_0_guess[0]
else:
    E_0_guess = E_0

initial_guess = [E_0_guess, rho_0_guess, P_0_guess]
print(f"Initial Guess: (E_0, rho_0, P_0) = {initial_guess}")

# Gather the binned data
data = np.hstack((rho_bin, P_bin, v_bin))

# Use the standard deviation (sigma) to weight the fit.
sigma_data = np.hstack((rho_sigma_bin, P_sigma_bin, v_sigma_bin))

# We must ensure no zero or near-zero sigmas
sigma_data = np.maximum(sigma_data, 1e-15)

# Transform them into log-space
initial_guess_log = [np.log(E_0_guess), np.log(rho_0_guess), np.log(P_0_guess)]

bounds_log = [
    (np.log(1e45), np.log(1e55)),
    (np.log(1e-30), np.log(1e-15)),
    (np.log(1e-20), np.log(1e-10)),
]

result = minimize(
    cost_function,
    initial_guess_log,
    args=(r_bin, data, sigma_data, time, gas_gamma),
    method="L-BFGS-B",
    bounds=bounds_log,
)

# Check if the minimization was successful
if result.success:
    log_popt = result.x
    E_0_fit = np.exp(log_popt[0])
    rho_0_fit = np.exp(log_popt[1])
    P_0_fit = np.exp(log_popt[2])

    # E_0_fit, rho_0_fit, P_0_fit = popt
    print("Fitted values: ", E_0_fit, rho_0_fit, P_0_fit)
else:
    print("Optimization failed. Reason:", result.message)
    E_0_fit, rho_0_fit, P_0_fit = E_0_guess, rho_0_guess, P_0_guess

# The main properties of the solution
r_s, P_s, rho_s, v_s, _, r_shock, _, _, _, _ = sedov(
    time, E_0_fit, rho_0_fit, gas_gamma, 1000, 3
)

# Append points for after the shock
r_shock = np.asarray(r_shock).flatten()  # Ensures it's 1D
r_s = np.insert(r_s, np.size(r_s), [r_shock[0], r_shock[0] * 1.5])
rho_s = np.insert(rho_s, np.size(rho_s), [rho_0_fit, rho_0_fit])
P_s = np.insert(P_s, np.size(P_s), [P_0_fit, P_0_fit])
v_s = np.insert(v_s, np.size(v_s), [0, 0])

# Additional arrays
u_s = P_s / (rho_s * (gas_gamma - 1.0))  # internal energy

#################################
# Plot the interesting quantities
#################################

# --- Unit definitions for plotting ---
# Define the target units for the plots
# This makes it easy to change units in one place if needed
target_units = {
    "length": units.kpc,
    "velocity": units.km / units.s,
    "density": units.g / units.cm**3,
    "internal_energy": units.erg,
    "pressure": units.Pa,
}

# --- Conversion to plotting units ---
# Convert all simulation and analytical data from CGS to plotting units

# Convert length from cm to kpc
r_plot = (r * units.cm).to(target_units["length"]).value
r_bin_plot = (r_bin * units.cm).to(target_units["length"]).value
r_s_plot = (r_s * units.cm).to(target_units["length"]).value
r_shock_plot = (r_shock * units.cm).to(target_units["length"]).value

# Convert velocity from cm/s to km/s
v_r_plot = (v_r * units.cm / units.s).to(target_units["velocity"]).value
v_bin_plot = (v_bin * units.cm / units.s).to(target_units["velocity"]).value
v_s_plot = (v_s * units.cm / units.s).to(target_units["velocity"]).value
v_sigma_bin_plot = (v_sigma_bin * units.cm / units.s).to(target_units["velocity"]).value

# Convert density from g/cm^3 to g/cm^3 (no change, but explicit for clarity)
rho_plot = (rho * units.g / units.cm**3).to(target_units["density"]).value
rho_bin_plot = (rho_bin * units.g / units.cm**3).to(target_units["density"]).value
rho_s_plot = (rho_s * units.g / units.cm**3).to(target_units["density"]).value
rho_sigma_bin_plot = (
    (rho_sigma_bin * units.g / units.cm**3).to(target_units["density"]).value
)

# Convert specific internal energy to erg/s
u_plot = (u * (units.cm / units.s) ** 2).to(units.erg / units.g).value
u_bin_plot = (u_bin * (units.cm / units.s) ** 2).to(units.erg / units.g).value
u_s_plot = (u_s * (units.cm / units.s) ** 2).to(units.erg / units.g).value
u_sigma_bin_plot = (
    (u_sigma_bin * (units.cm / units.s) ** 2).to(units.erg / units.g).value
)

# Convert pressure from CGS to Pascals (Pa)
P_plot = (P * units.g / (units.cm * units.s**2)).to(target_units["pressure"]).value
P_bin_plot = (
    (P_bin * units.g / (units.cm * units.s**2)).to(target_units["pressure"]).value
)
P_s_plot = (P_s * units.g / (units.cm * units.s**2)).to(target_units["pressure"]).value
P_sigma_bin_plot = (
    (P_sigma_bin * units.g / (units.cm * units.s**2)).to(target_units["pressure"]).value
)

# Define colors and plot properties
line_color = "blueviolet"
binned_color = "tomato"
scatter_color = "royalblue"
binned_marker_size = 4

scatter_props = dict(
    marker=".",
    ms=1,
    markeredgecolor="none",
    alpha=0.8,
    zorder=-1,
    rasterized=True,
    linestyle="none",
    color=scatter_color,
)

errorbar_props = dict(color=binned_color, ms=binned_marker_size, fmt=".", lw=1.2)

# Create figure with subplots
fig, axes = plt.subplots(1, 4, figsize=(16, 4))

# Velocity profile
ax = axes[0]
ax.plot(r_plot, v_r_plot, **scatter_props)
ax.plot(r_s_plot, v_s_plot, "--", color=line_color, alpha=0.8, lw=1.2)
ax.errorbar(r_bin_plot, v_bin_plot, yerr=v_sigma_bin_plot, **errorbar_props)
ax.set_xlabel("$r$ [kpc]")
ax.set_ylabel("$v_r$ [km/s]")
ax.set_xlim(0, n_shock * r_shock_plot)

# Density profile
ax = axes[1]
ax.plot(r_plot, rho_plot, **scatter_props)
ax.plot(r_s_plot, rho_s_plot, "--", color=line_color, alpha=0.8, lw=1.2)
ax.errorbar(r_bin_plot, rho_bin_plot, yerr=rho_sigma_bin_plot, **errorbar_props)
ax.set_xlabel("$r$ [kpc]")
ax.set_ylabel(r"$\rho$ [g/cm$^3$]")
ax.set_xlim(0, n_shock * r_shock_plot)

# Internal energy profile
ax = axes[2]
ax.plot(r_plot, u_plot, **scatter_props)
ax.plot(r_s_plot, u_s_plot, "--", color=line_color, alpha=0.8, lw=1.2)
ax.errorbar(r_bin_plot, u_bin_plot, yerr=u_sigma_bin_plot, **errorbar_props)
ax.set_xlabel("$r$ [kpc]")
ax.set_ylabel(r"$u_{\mathrm{int}}$ [erg/g]")
ax.set_xlim(0, n_shock * r_shock_plot)

# Pressure profile
ax = axes[3]
ax.plot(r_plot, P_plot, **scatter_props)
ax.plot(r_s_plot, P_s_plot, "--", color=line_color, alpha=0.8, lw=1.2)
ax.errorbar(r_bin_plot, P_bin_plot, yerr=P_sigma_bin_plot, **errorbar_props)
ax.set_xlabel("$r$ [kpc]")
ax.set_ylabel("$P$ [Pa]")
ax.set_xlim(0, n_shock * r_shock_plot)

# Adjust layout and save figure
plt.tight_layout()

if output_filename is not None:
    plt.savefig(output_filename, format=image_format, dpi=dpi)
else:
    plt.show()
