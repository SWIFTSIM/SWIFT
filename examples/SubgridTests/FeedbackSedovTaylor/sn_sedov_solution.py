###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
#                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
import astropy.units as units
from scipy.optimize import leastsq

# %%


def calc_a(g, nu=3):
    """
    exponents of the polynomials of the sedov solution
    g - the polytropic gamma
    nu - the dimension
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
    beta values for the sedov solution (coefficients of the polynomials of the similarity variables)
    v - the similarity variable
    g - the polytropic gamma
    nu- the dimension
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
    solve the sedov problem
    t - the time
    E0 - the initial energy
    rho0 - the initial density
    n - number of points (10000)
    nu - the dimension
    g - the polytropic gas gamma
    """
    # the similarity variable
    v_min = 2.0 / ((nu + 2) * g)
    v_max = 4.0 / ((nu + 2) * (g + 1))

    v = v_min + np.arange(n) * (v_max - v_min) / (n - 1.0)

    a = calc_a(g, nu)
    beta = calc_beta(v, g=g, nu=nu)
    lbeta = np.log(beta)

    r = np. exp(-a[0] * lbeta[0] - a[2] * lbeta[1] - a[1] * lbeta[2])
    rho = ((g + 1.0) / (g - 1.0)) * np.exp(
        a[3] * lbeta[1] + a[5] * lbeta[3] + a[4] * lbeta[2]
    )
    p = np.exp(nu * a[0] * lbeta[0] + (a[5] + 1) *
               lbeta[3] + (a[4] - 2 * a[1]) * lbeta[2])
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
    fac = (q * (t ** nu) * rho0 / E0) ** (-1.0 / (nu + 2))

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
    return r, p, rho, u, r_s, p_s, rho_s, u_s, shock_speed

# %%


class ExplicitDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Formatter that prints default arguments if provided."""

    def _get_help_string(self, action):
        if isinstance(action.default, np.ndarray) or action.default not in (None, False):
            return super()._get_help_string(action)
        return action.help


class RawTextArgumentDefaultsHelpFormatter(ExplicitDefaultsHelpFormatter,
                                           argparse.RawTextHelpFormatter):
    """Combines raw text formatting with default argument printing."""
    pass


def parse_options():
    """
    Parses command-line arguments for computing the analytical solution
    of the 3D Sedov blast wave and comparing it to simulation data.
    """
    parser = argparse.ArgumentParser(
        description="Compute and analyze the Sedov blast wave from a SWIFT simulation.",
        epilog="Example usage:\npython script.py output.hdf5 --rho_0 1.0 --E_0 1e51",
        formatter_class=RawTextArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "file", type=str, help="Path to the SWIFT output file.")

    parser.add_argument("--rho_0", type=float, default=1.0,
                        help="Background density (in simulation units).")
    parser.add_argument("--P_0", type=float, default=1e-6,
                        help="Background pressure (in simulation units).")
    parser.add_argument("--E_0", type=float, default=1,
                        help="Explosion energy (in simulation units).")

    parser.add_argument("-o", type=str, dest="output_filename", default=None,
                        help="Output filename for the plot. If not set, the plot is shown instead.")

    parser.add_argument("--image_format", type=str, default="png",
                        help="Format of the output image (recognized by Matplotlib).")
    parser.add_argument("--dpi", type=float, default=300,
                        help="Resolution of output images.")
    parser.add_argument("--step", type=float, default=0.01,
                        help="Step for binning.")
    parser.add_argument("--stylesheet", type=str,
                        default=None, help="Matplotlib stylesheet.")

    return parser.parse_args()


# %%
# Parse command-line arguments
opt = parse_options()

# Apply Matplotlib stylesheet
if opt.stylesheet:
    plt.style.use(opt.stylesheet)

# Physical parameters
rho_0 = opt.rho_0  # Background Density
P_0 = opt.P_0      # Background Pressure
E_0 = opt.E_0      # Energy of the explosion
gas_gamma = 5.0 / 3.0  # Gas polytropic index
n_shock = 1.5  # Shock index
step = opt.step

# Ensure the file exists
file = opt.file
try:
    with h5py.File(file, "r") as sim:
        # Read unit conversion factors from the snapshot
        try:
            U_M = sim["/InternalCodeUnits"].attrs.get("Unit mass in cgs (U_M)", 1.0)[0]
            U_L = sim["/InternalCodeUnits"].attrs.get("Unit length in cgs (U_L)", 1.0)[0]
            U_t = sim["/InternalCodeUnits"].attrs.get("Unit time in cgs (U_t)", 1.0)[0]
            
        except KeyError:
            print("Warning: /Units group not found. Assuming CGS units.")
            pass # Continue with default CGS-like conversion factors

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
        x, y, z = pos[:, 0] - box_size / 2, pos[:, 1] - \
            box_size / 2, pos[:, 2] - box_size / 2
        r = np.sqrt(x**2 + y**2 + z**2)

        # Compute radial velocity component
        v_r = (x * vel[:, 0] + y * vel[:, 1] + z * vel[:, 2]) / r

        # --- CONVERT ALL DATA TO CGS IN PLACE ---
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
        
        # We also need to convert the initial guess values from the command line
        rho_0 = rho_0 * U_M / (U_L**3)
        P_0 = P_0 * U_M / (U_L * U_t**2)
        E_0 = E_0 * U_M * U_L**2 / U_t**2
        step = step * U_L

except FileNotFoundError:
    sys.exit(f"Error: File '{file}' not found.")
except KeyError as e:
    sys.exit(f"Error: Missing expected dataset {e} in file '{file}'.")
except Exception as e:
    sys.exit(f"Unexpected error: {e}")
    
# Bin the data
r_bin_edge = np.arange(0.0, r.max(), step)
r_bin = 0.5 * (r_bin_edge[1:] + r_bin_edge[:-1])
rho_bin, _, _ = stats.binned_statistic(
    r, rho, statistic="mean", bins=r_bin_edge)
v_bin, _, _ = stats.binned_statistic(r, v_r, statistic="mean", bins=r_bin_edge)
P_bin, _, _ = stats.binned_statistic(r, P, statistic="mean", bins=r_bin_edge)
S_bin, _, _ = stats.binned_statistic(r, S, statistic="mean", bins=r_bin_edge)
u_bin, _, _ = stats.binned_statistic(r, u, statistic="mean", bins=r_bin_edge)
rho2_bin, _, _ = stats.binned_statistic(
    r, rho ** 2, statistic="mean", bins=r_bin_edge)
v2_bin, _, _ = stats.binned_statistic(
    r, v_r ** 2, statistic="mean", bins=r_bin_edge)
P2_bin, _, _ = stats.binned_statistic(
    r, P ** 2, statistic="mean", bins=r_bin_edge)
S2_bin, _, _ = stats.binned_statistic(
    r, S ** 2, statistic="mean", bins=r_bin_edge)
u2_bin, _, _ = stats.binned_statistic(
    r, u ** 2, statistic="mean", bins=r_bin_edge)

# Remove unwanted values
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

# Compute sigma
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin ** 2)
v_sigma_bin = np.sqrt(v2_bin - v_bin ** 2)
P_sigma_bin = np.sqrt(P2_bin - P_bin ** 2)
u_sigma_bin = np.sqrt(u2_bin - u_bin ** 2)


################################
# Now, work our the solution....
################################

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

def sedov_2(t, E0, rho0, P0, g, n=1000, nu=3):
    """
    solve the sedov problem
    t - the time
    E0 - the initial energy
    rho0 - the initial density
    n - number of points (10000)
    nu - the dimension
    g - the polytropic gas gamma
    """
    n = n-2
    
    # the similarity variable
    v_min = 2.0 / ((nu + 2) * g)
    v_max = 4.0 / ((nu + 2) * (g + 1))

    v = v_min + np.arange(n) * (v_max - v_min) / (n - 1.0)

    a = calc_a(g, nu)
    beta = calc_beta(v, g=g, nu=nu)
    lbeta = np.log(beta)

    r = np. exp(-a[0] * lbeta[0] - a[2] * lbeta[1] - a[1] * lbeta[2])
    rho = ((g + 1.0) / (g - 1.0)) * np.exp(
        a[3] * lbeta[1] + a[5] * lbeta[3] + a[4] * lbeta[2]
    )
    p = np.exp(nu * a[0] * lbeta[0] + (a[5] + 1) *
               lbeta[3] + (a[4] - 2 * a[1]) * lbeta[2])
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
    fac = (q * (t ** nu) * rho0 / E0) ** (-1.0 / (nu + 2))

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
    A wrapper for the Sedov function to be used in curve fitting.
    Returns interpolated values of rho, P, and v at given radii r.
    """

    # Solve Sedov blast wave problem
    r_s, P_s, rho_s, v_s, r_shock, _, _, _, _ = sedov_2(
        t=time,  
        E0=E_0,
        rho0=rho_0,
        P0=P_0,
        g=g,
        n=len(r),  
        nu=nu
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

rho_0_guess = np.median(rho)
P_0_guess = np.median(P)

# Step 1: Create a more robust initial guess from the data itself.
# Find the point of maximum change in density to estimate the shock radius.
d_rho_bin = np.gradient(rho_bin)
r_shock_guess = r_bin[np.argmax(d_rho_bin)]

# Guess E_0 using a known analytical relationship for a Sedov explosion.
# nu is the dimension, which is 3 here.
nu = 3
# A dimensionless constant for E_0, approx 1.15 for nu=3, g=5/3
E_0_guess = rho_0_guess * (r_shock_guess ** 5) / (time ** 2)
E_0_guess = E_0_guess[0]

initial_guess = [E_0_guess, rho_0_guess, P_0_guess]

print("Improved Initial Guess:", initial_guess)

data = np.hstack((rho_bin, P_bin, v_bin))


# Check for NaNs and Infs in your data
print("Number of NaNs in rho_bin:", np.sum(np.isnan(rho_bin)))
print("Number of NaNs in P_bin:", np.sum(np.isnan(P_bin)))
print("Number of NaNs in v_bin:", np.sum(np.isnan(v_bin)))
print("Number of Infs in data:", np.sum(np.isinf(data)))

valid_mask = np.isfinite(rho_bin) & np.isfinite(P_bin) & np.isfinite(v_bin)
r_bin = r_bin[valid_mask]
rho_bin = rho_bin[valid_mask]
P_bin = P_bin[valid_mask]
v_bin = v_bin[valid_mask]

print("Number of NaNs in rho_bin:", np.sum(np.isnan(rho_bin)))
print("Number of NaNs in P_bin:", np.sum(np.isnan(P_bin)))
print("Number of NaNs in v_bin:", np.sum(np.isnan(v_bin)))
print("Number of Infs in data:", np.sum(np.isinf(data)))

data = np.hstack((rho_bin, P_bin, v_bin))

# Step 2: Use the standard deviation (sigma) to weight the fit.
# This guides the fit to prioritize regions with less scatter.
sigma_data = np.hstack((rho_sigma_bin, P_sigma_bin, v_sigma_bin))

# We must ensure no zero or near-zero sigmas
sigma_data = np.maximum(sigma_data, 1e-15)

popt, pcov = curve_fit(
        lambda r_fit, E_0, rho_0, P_0: sedov_fit_model(r_fit, E_0, rho_0, P_0, time, gas_gamma, nu),
        r_bin,
        data,
        p0=initial_guess,
        sigma=sigma_data,  # Pass the uncertainties here
        absolute_sigma=True, # Treat sigma as absolute errors
    )

# Extract the fitted parameters
E_0_fit, rho_0_fit, P_0_fit = popt

print("Initial guess: ", E_0, rho_0, P_0)
print("Fitted values: ", E_0_fit, rho_0_fit, P_0_fit)

E_0, rho_0, P_0 = E_0_fit, rho_0_fit, P_0_fit

# The main properties of the solution
r_s, P_s, rho_s, v_s, r_shock, _, _, _, _ = sedov(
    time, E_0, rho_0, gas_gamma, 1000, 3)

# Append points for after the shock
r_shock = np.asarray(r_shock).flatten()  # Ensures it's 1D
r_s = np.insert(r_s, np.size(r_s), [r_shock[0], r_shock[0] * 1.5])
rho_s = np.insert(rho_s, np.size(rho_s), [rho_0, rho_0])
P_s = np.insert(P_s, np.size(P_s), [P_0, P_0])
v_s = np.insert(v_s, np.size(v_s), [0, 0])

# Additional arrays
u_s = P_s / (rho_s * (gas_gamma - 1.0))  # internal energy

#################################
# Plot the interesting quantities
#################################
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

errorbar_props = dict(color=binned_color,
                      ms=binned_marker_size, fmt=".", lw=1.2)

# Function for unit conversion
def convert_units(quantity, unit):
    return (quantity * 1e10 * units.Msun / units.kpc**3).to(unit).value

# Convert density-related quantities
rho_s = convert_units(rho_s, units.g / units.cm**3)
rho_bin = convert_units(rho_bin, units.g / units.cm**3)
rho = convert_units(rho, units.g / units.cm**3)
rho_sigma_bin = convert_units(rho_sigma_bin, units.g / units.cm**3)

# Convert internal energy-related quantities
def convert_energy(quantity):
    return (quantity * 1e10 * units.Msun * units.kpc**2 / units.Gyr**2).to(units.erg).value

u = convert_energy(u)
u_s = convert_energy(u_s)
u_bin = convert_energy(u_bin)
u_sigma_bin = convert_energy(u_sigma_bin)

# Create figure with subplots
fig, axes = plt.subplots(1, 4, figsize=(16, 4))

# Velocity profile
ax = axes[0]
ax.plot(r, v_r, **scatter_props)
ax.plot(r_s, v_s, "--", color=line_color, alpha=0.8, lw=1.2)
ax.errorbar(r_bin, v_bin, yerr=v_sigma_bin, **errorbar_props)
ax.set_xlabel("$r$ [kpc]")
ax.set_ylabel("$v_r$ [km/s]")
ax.set_xlim(0, n_shock * r_shock)

# Density profile
ax = axes[1]
ax.plot(r, rho, **scatter_props)
ax.plot(r_s, rho_s, "--", color=line_color, alpha=0.8, lw=1.2)
ax.errorbar(r_bin, rho_bin, yerr=rho_sigma_bin, **errorbar_props)
ax.set_xlabel("$r$ [kpc]")
ax.set_ylabel(r"$\rho$ [g/cm$^3$]")
ax.set_xlim(0, n_shock * r_shock)

# Internal energy profile
ax = axes[2]
ax.plot(r, u, **scatter_props)
ax.plot(r_s, u_s, "--", color=line_color, alpha=0.8, lw=1.2)
ax.errorbar(r_bin, u_bin, yerr=u_sigma_bin, **errorbar_props)
ax.set_xlabel("$r$ [kpc]")
ax.set_ylabel("$u_{\mathrm{int}}$ [erg]")
ax.set_xlim(0, n_shock * r_shock)

# Pressure profile
ax = axes[3]
ax.plot(r, P, **scatter_props)
ax.plot(r_s, P_s, "--", color=line_color, alpha=0.8, lw=1.2)
ax.errorbar(r_bin, P_bin, yerr=P_sigma_bin, **errorbar_props)
ax.set_xlabel("$r$ [kpc]")
ax.set_ylabel("$P$")
ax.set_xlim(0, n_shock * r_shock)

# Adjust layout and save figure
plt.tight_layout()
plt.savefig("Sedov.pdf", format="pdf", dpi=300)
