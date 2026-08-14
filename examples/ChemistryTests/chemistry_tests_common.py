#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################################################################################
# This file is part of SWIFT.
# Copyright (c)  2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
################################################################################
"""Helpers shared by the ChemistryTests metal_profile.py analysis scripts."""
import numpy as np
from scipy.special import ive as scaled_bessel


def get_fe_metal_mass(data):
    """Get the Fe metal mass (in code units) from swiftsimio data."""
    if hasattr(data.gas.metal_mass_fractions, "fe"):
        m_fe = data.gas.metal_mass_fractions.fe * data.gas.masses
    else:  # This case happens when we run the simulation without --feedback
        m_fe = data.gas.metal_mass_fractions[:, 0] * data.gas.masses
    return m_fe


def lattice_bin_count(x0, L, window=None, cap=512, floor=16):
    """Bin count matched to the seed lattice spacing over `window` (defaults
    to the full box L; pass the actual plotted width if it's a sub-range),
    else bins that don't evenly divide the static IC grid alias into
    spurious spikes: a bin boundary landing mid-lattice-spacing pulls two
    planes' worth of mass into one bin while its neighbour gets none."""
    ux = np.unique(np.round(x0 % L, 8))
    if len(ux) < 2:
        return cap
    spacing = np.median(np.diff(ux))
    if spacing <= 0:
        return cap
    if window is None:
        window = L
    return int(np.clip(round(window / spacing), floor, cap))


def radial_profile(value, r, r_max, cross_section_area, n_bins=30):
    """Compute a 1D density profile (mass / bin volume) across -r_max to
    r_max. Summing mass and dividing by bin volume (rather than averaging
    per-particle values) makes the result resolution-independent -- a
    per-particle mean shrinks by construction as particle count grows,
    even when the underlying physical density is converged."""
    r_min = -r_max
    r_bins = np.linspace(r_min, r_max, n_bins + 1)
    bin_width = (r_max - r_min) / n_bins
    bin_volume = bin_width * cross_section_area
    bin_indices = np.digitize(r, bins=r_bins) - 1
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    densities = np.zeros(n_bins)
    for i in range(n_bins):
        in_bin = bin_indices == i
        if np.any(in_bin):
            densities[i] = np.sum(value[in_bin]) / bin_volume
    return r_centers, densities


def hyperbolic_diffusion_solution(x, t, q_0, x_0, tau, kappa):
    """Smooth interior (Bessel) part of the exact 1D point-source Green's
    function of the Cattaneo (telegrapher) equation, unit-normalised by q_0."""
    c = np.sqrt(kappa / tau)
    Delta_x = x - x_0
    within_causal_region = np.abs(Delta_x) <= c * t
    radial_term = np.zeros_like(Delta_x)
    radial_term[within_causal_region] = np.sqrt(
        c**2 * t**2 - Delta_x[within_causal_region] ** 2
    )
    # Scaled Bessels: iv(n,z)*exp(-at) overflows for t >> tau (z ~ at can
    # exceed 700); ive(n,z)*exp(z-at) is exact and safe since z <= at.
    z = (c / (2 * kappa)) * radial_term[within_causal_region]
    at = c**2 * t / (2 * kappa)
    scale = np.exp(z - at)
    I0 = scaled_bessel(0, z) * scale
    I1 = scaled_bessel(1, z) * scale
    u = np.zeros_like(Delta_x)
    if t > 0:
        u[within_causal_region] = (
            0.5
            * q_0
            * (
                (c / (2 * kappa)) * I0
                + (c**2 / (2 * kappa)) * t * I1 / radial_term[within_causal_region]
            )
        )
    return u


def hyperbolic_diffusion_solution_convolved(
    r,
    t,
    total_mass,
    epsilon,
    x_0,
    tau,
    kappa,
    n_source_points=401,
    with_front_term=False,
    source_shape="segment",
):
    """Exact response to a finite seed of total mass total_mass and radius
    epsilon centred at x_0 (not a point source): convolve
    hyperbolic_diffusion_solution() against the seed's x-marginal, and add
    the ballistic front delta (Masoliver, Porra & Weiss, Phys. Rev. E 48,
    939 (1993)) at x_0 +/- c*t with weight 1/2 * exp(-t/2tau).
    source_shape: 'segment' (1D top-hat; uniform weight 1/(2 eps)) or
    'ball' (3D uniform sphere; x-marginal (3/(4 eps))(1-(xi/eps)^2)). For a
    3D run profiled along x this 1D convolution is exact: the y,z-marginal
    of the 3D Green's function is the 1D Green's function. Self-normalising
    (uses total_mass directly, no amplitude fit needed)."""
    c = np.sqrt(kappa / tau)
    if source_shape == "ball":

        def w_src(xi):
            return 0.75 / epsilon * np.maximum(1 - (xi / epsilon) ** 2, 0.0)

    else:

        def w_src(xi):
            return np.where(np.abs(xi) < epsilon, 0.5 / epsilon, 0.0)

    x_prime = np.linspace(-epsilon, epsilon, n_source_points)
    w = w_src(x_prime)
    if t == 0:
        return total_mass * w_src(r - x_0)
    u = np.zeros_like(r)
    for i, ri in enumerate(r):
        g = hyperbolic_diffusion_solution(ri - x_0 - x_prime, t, 1.0, 0.0, tau, kappa)
        u[i] = total_mass * np.trapezoid(w * g, x_prime)
    if with_front_term:
        front_weight = 0.5 * np.exp(-t / (2 * tau)) * total_mass
        for front_x0 in (x_0 - c * t, x_0 + c * t):
            u += front_weight * w_src(r - front_x0)
    return u
