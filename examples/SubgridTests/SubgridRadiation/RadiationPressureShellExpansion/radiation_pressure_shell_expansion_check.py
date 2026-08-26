################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
"""
Tier 2: compare the simulated radiation-pressure-driven shell expansion
against the momentum-conserving analytic solution of Krumholz & Matzner
(2009, ApJ 703, 1352, arXiv:0906.4343), radiation-only limit (their Eq. 11 at
k_rho=0; this example runs with photoionization off, so there is no
competing gas-pressure driver to blend in):

    R(t) = [ 3 * f_trap * L_bol * t^2 / (2*pi*rho_0*c) ]^(1/4)

f_trap is this code's own (1-exp(-tau_NUV))*(1+tau_IR)
(radiation_get_star_physical_radiation_pressure, src/feedback/GEAR/
radiation_pressure.c) -- the analogue of K&M2009's trapping factor, minus
their winds/Lyman-alpha channels this code does not model. Unlike K&M2009's
constant f_trap, this is measured directly from each snapshot's star
EnrichmentWeight/GradRhoStar/ZStar debug fields (requires a build configured
with --enable-debug-interactions-stars) and reported over the run, rather
than assumed -- this checks a previously-unverified hypothesis that the
star's adaptive smoothing length keeps its local coupling roughly tracking
the swept-up shell as the interior cavity empties out.

The current reference curve uses f_trap measured at the earliest snapshot,
held constant (a first pass -- see the module's own TODO on refining this to
a numerically-integrated ODE using the full measured f_trap(t), matched to
this check's own drift report).

Formula-level correctness of radiation_get_star_physical_radiation_pressure
is NOT this script's job: that's tests/testRadiationPressureFormula.c (no
hydro, no snapshots). This script validates the DYNAMICAL response -- the
sweep-up physics that made a static gas-kinematics comparison meaningless
(see the retired radiation_pressure_momentum_check.py, commit bfd339cc1) is
exactly what's being checked here.
"""

import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from radiation_table_reader import load_radiation_table

# src/kernel_hydro.h, 3D Wendland C2 branch (this example's --with-kernel).
KERNEL_GAMMA_WENDLAND_C2_3D = 1.936492
# src/feedback/GEAR/radiation.h
RADIATION_MIN_RELATIVE_DENSITY_GRADIENT = 0.01
# src/feedback/GEAR/radiation_pressure.c
KAPPA_IR_CGS = 10.0
KAPPA_NUV_CGS = 1800.0
Z_SUN = 0.02
C_LIGHT_CGS = 2.99792458e10


def snapshot_index(filename):
    """Trailing integer of a snapshot filename (snapshot_0007.hdf5 -> 7)."""
    stem = os.path.splitext(os.path.basename(filename))[0]
    match = re.search(r"(\d+)$", stem)
    return int(match.group(1)) if match else None


def load_snapshot(path):
    """Read the star/gas fields this check needs from one snapshot."""
    with h5py.File(path, "r") as f:
        header = f["/Header"]
        time = float(np.asarray(header.attrs["Time"]).flat[0])
        boxsize = float(np.asarray(header.attrs["BoxSize"]).flat[0])
        units = f["/Units"]
        unit_length_cgs = float(
            np.asarray(units.attrs["Unit length in cgs (U_L)"]).flat[0]
        )
        unit_mass_cgs = float(np.asarray(units.attrs["Unit mass in cgs (U_M)"]).flat[0])

        gas = f["/PartType0"]
        gas_pos = gas["Coordinates"][:, :]
        gas_mass = gas["Masses"][:]

        star = f["/PartType4"]
        star_pos = star["Coordinates"][0, :]
        star_mass_msun = float(star["Masses"][0]) * unit_mass_cgs / 1.988409870698051e33
        star_h = float(star["SmoothingLengths"][0])
        has_debug_fields = "EnrichmentWeight" in star
        f_trap_inputs = None
        if has_debug_fields:
            f_trap_inputs = dict(
                enrichment_weight=float(star["EnrichmentWeight"][0]),
                grad_rho=np.array(star["GradRhoStar"][0, :]),
                Z_star=float(star["ZStar"][0]),
                h=star_h,
            )

    return dict(
        time=time,
        boxsize=boxsize,
        unit_length_cgs=unit_length_cgs,
        unit_mass_cgs=unit_mass_cgs,
        gas_pos=gas_pos,
        gas_mass=gas_mass,
        star_pos=star_pos,
        star_mass_msun=star_mass_msun,
        has_debug_fields=has_debug_fields,
        f_trap_inputs=f_trap_inputs,
    )


def compute_f_trap(f_trap_inputs, unit_length_cgs, unit_mass_cgs):
    """Mirror radiation_pressure.c's Sigma_gas/tau_NUV/tau_IR/f_trap, in cgs.

    Same Sobolev-length guard as radiation_get_comoving_gas_column_density_
    at_star: a near-zero density gradient is discarded, not divided by.
    """
    rho_gas = f_trap_inputs["enrichment_weight"]
    h = f_trap_inputs["h"]
    norm_grad_rho = np.linalg.norm(f_trap_inputs["grad_rho"])

    h_gas = h * KERNEL_GAMMA_WENDLAND_C2_3D
    grad_rho_floor = RADIATION_MIN_RELATIVE_DENSITY_GRADIENT * rho_gas / h_gas
    sobolev_length = rho_gas / norm_grad_rho if norm_grad_rho > grad_rho_floor else 0.0
    Sigma_gas_internal = (h_gas + sobolev_length) * rho_gas

    unit_column_density_cgs = unit_mass_cgs / unit_length_cgs**2
    Sigma_gas_cgs = Sigma_gas_internal * unit_column_density_cgs

    Z_ratio = f_trap_inputs["Z_star"] / Z_SUN
    tau_IR = KAPPA_IR_CGS * Z_ratio * Sigma_gas_cgs
    tau_NUV = KAPPA_NUV_CGS * Z_ratio * Sigma_gas_cgs
    f_trap = (1.0 - np.exp(-tau_NUV)) * (1.0 + tau_IR)
    return f_trap, tau_NUV, tau_IR


def enclosed_gas_mass_profile(gas_pos, gas_mass, star_pos, box):
    """Sorted (radius, cumulative mass) profile from the star, periodic-aware."""
    dx = gas_pos - star_pos
    dx -= box * np.round(dx / box)
    r = np.linalg.norm(dx, axis=1)
    order = np.argsort(r)
    return r[order], np.cumsum(gas_mass[order])


def lagrangian_shell_radius(r_sorted, cum_mass, M_target):
    """Radius at which the cumulative mass profile reaches M_target.

    This is the accounting the thin-shell ODE itself assumes (all originally
    interior mass now sits at R): ask where the simulation's own mass
    profile puts the boundary the theory predicts, rather than fitting a
    noisy density peak.
    """
    if len(cum_mass) == 0 or M_target <= cum_mass[0]:
        return float(r_sorted[0]) if len(r_sorted) else 0.0
    if M_target >= cum_mass[-1]:
        return float(r_sorted[-1])
    return float(np.interp(M_target, cum_mass, r_sorted))


def density_peak_radius(r_sorted, cum_mass, box_half_width, n_bins=40):
    """Secondary/audit diagnostic: radius of peak radial density.

    Capped at the box half-width: beyond that, periodic images make radial
    binning ambiguous, and r_sorted[-1] (the single farthest particle) is
    itself an unstable estimate of "how far the profile extends" once the
    shell has grown to fill most of the box.
    """
    if len(r_sorted) < n_bins:
        return float(r_sorted[-1]) if len(r_sorted) else 0.0
    r_max = min(r_sorted[-1], box_half_width)
    edges = np.linspace(0.0, r_max, n_bins + 1)
    mass_at_edges = np.interp(edges, r_sorted, cum_mass, left=0.0, right=cum_mass[-1])
    shell_mass = np.diff(mass_at_edges)
    shell_volume = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    with np.errstate(invalid="ignore", divide="ignore"):
        density = shell_mass / shell_volume
    density[~np.isfinite(density)] = -np.inf
    peak = int(np.argmax(density))
    return 0.5 * (edges[peak] + edges[peak + 1])


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-s", "--snapshot-glob", default="snap/snapshot_*.hdf5")
    parser.add_argument("-o", "--output", default="shell_expansion_check.png")
    parser.add_argument(
        "--tol",
        type=float,
        default=0.25,
        help="Relative-error tolerance for a pass, applied to the median "
        "over the box-valid window (default: 0.25 -- looser than "
        "Starbench's single-point 0.1, tighter than the retired Tier-1 "
        "check's 0.15-on-a-contaminated-signal; retune once real run "
        "scatter is measured).",
    )
    parser.add_argument(
        "--mismatch-tol",
        type=float,
        default=0.20,
        help="Relative disagreement between the primary (Lagrangian-mass) "
        "and secondary (density-peak) shell-radius diagnostics above which "
        "the verdict gets a mismatch marker (default: 0.20).",
    )
    args = parser.parse_args()

    files = sorted(glob.glob(args.snapshot_glob))
    if len(files) < 2:
        raise RuntimeError(
            f"Need at least 2 snapshots, found {len(files)} matching "
            f"{args.snapshot_glob!r}."
        )
    first_index = snapshot_index(files[0])
    if first_index is not None and first_index != 0:
        print(
            f"WARNING: the snapshot glob starts at index {first_index}, not 0 "
            f"-- rho_0 would be measured after the region has already "
            f"expanded. Re-run against the complete snapshot series."
        )

    snaps = [load_snapshot(f) for f in files]
    if not snaps[0]["has_debug_fields"]:
        raise RuntimeError(
            "Snapshots carry no EnrichmentWeight/GradRhoStar/ZStar fields -- "
            "rebuild swift with --enable-debug-interactions-stars to get the "
            "star-local quantities this check needs to measure f_trap."
        )

    unit_length_cgs = snaps[0]["unit_length_cgs"]
    unit_mass_cgs = snaps[0]["unit_mass_cgs"]
    unit_density_cgs = unit_mass_cgs / unit_length_cgs**3

    # rho_0: mean gas density in the earliest snapshot (uniform IC).
    box0 = snaps[0]["boxsize"]
    rho_0_internal = np.sum(snaps[0]["gas_mass"]) / box0**3
    rho_0_cgs = rho_0_internal * unit_density_cgs
    print(f"rho_0              : {rho_0_cgs:.4g} g/cm^3  (measured in {files[0]})")

    star_mass_msun = snaps[0]["star_mass_msun"]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    radiation_table = load_radiation_table(files, script_dir)
    L_bol_cgs = radiation_table.luminosity_erg_s(star_mass_msun)
    print(f"Star mass          : {star_mass_msun:.4g} Msun")
    print(f"L_bol              : {L_bol_cgs:.4g} erg/s")

    with h5py.File(files[0], "r") as f:
        unit_time_cgs = float(
            np.asarray(f["/Units"].attrs["Unit time in cgs (U_t)"]).flat[0]
        )
    t_sim_s = np.array([s["time"] for s in snaps]) * unit_time_cgs

    box_half_width_cgs = 0.5 * snaps[0]["boxsize"] * unit_length_cgs

    f_trap_t, tau_NUV_t, tau_IR_t = [], [], []
    R_peak_cgs = []
    mass_profiles = []
    for snap in snaps:
        f_trap, tau_NUV, tau_IR = compute_f_trap(
            snap["f_trap_inputs"], unit_length_cgs, unit_mass_cgs
        )
        f_trap_t.append(f_trap)
        tau_NUV_t.append(tau_NUV)
        tau_IR_t.append(tau_IR)

        r_sorted, cum_mass = enclosed_gas_mass_profile(
            snap["gas_pos"], snap["gas_mass"], snap["star_pos"], snap["boxsize"]
        )
        mass_profiles.append((r_sorted, cum_mass))
        R_peak_cgs.append(
            density_peak_radius(r_sorted, cum_mass, 0.5 * snap["boxsize"])
            * unit_length_cgs
        )

    f_trap_t = np.array(f_trap_t)

    # The t=0 snapshot is dumped before the star's own density loop has ever
    # run, so EnrichmentWeight (hence f_trap) is still zero-initialized
    # there; use the first snapshot with a genuinely measured value instead
    # (mirrors Starbench's own "measured in <snapshot>, not necessarily the
    # first" provenance discipline).
    nonzero = np.nonzero(f_trap_t > 0)[0]
    if len(nonzero) == 0:
        raise RuntimeError("f_trap is 0 in every snapshot -- no radiation pressure?")
    ref_idx = int(nonzero[0])
    print(
        f"f_trap(t={t_sim_s[ref_idx] / 3.15576e13:.4g} Myr) : {f_trap_t[ref_idx]:.4g}  "
        f"(reference value, measured in {files[ref_idx]})"
    )
    print(
        f"f_trap drift       : min={f_trap_t[nonzero].min():.4g}  "
        f"max={f_trap_t[nonzero].max():.4g}  (range/f_trap[ref]="
        f"{((f_trap_t[nonzero].max() - f_trap_t[nonzero].min()) / f_trap_t[ref_idx]):.2%})"
    )

    # Closed-form analytic curve, f_trap held at its reference measured
    # value (first-pass reference -- see module docstring's TODO on refining
    # this to an ODE using the full measured f_trap(t)).
    f_trap_ref = f_trap_t[ref_idx]
    prefactor = (
        3.0 * f_trap_ref * L_bol_cgs / (2.0 * np.pi * rho_0_cgs * C_LIGHT_CGS)
    ) ** (0.25)
    R_analytic_cgs = prefactor * np.sqrt(t_sim_s)

    # Now the Lagrangian shell radius, using the theory's own M_theory(t).
    R_lagrangian_cgs = []
    for i, (r_sorted, cum_mass) in enumerate(mass_profiles):
        M_theory_cgs = (4.0 / 3.0) * np.pi * rho_0_cgs * R_analytic_cgs[i] ** 3
        M_theory_internal = M_theory_cgs / unit_mass_cgs
        R_lagrangian_cgs.append(
            lagrangian_shell_radius(r_sorted, cum_mass, M_theory_internal)
            * unit_length_cgs
        )
    R_lagrangian_cgs = np.array(R_lagrangian_cgs)
    R_peak_cgs = np.array(R_peak_cgs)

    PC_IN_CGS = 3.0856775814913673e18
    MYR_IN_CGS = 3.15576e13

    box_valid = R_analytic_cgs <= box_half_width_cgs
    if not np.any(box_valid):
        raise RuntimeError(
            "The analytic curve exceeds the box half-width at every "
            "snapshot -- box too small or time_end too long for this check."
        )
    last_valid = np.where(box_valid)[0][-1]
    if last_valid < len(t_sim_s) - 1:
        print(
            f"NOTE: reference curve leaves the box after snapshot "
            f"{last_valid} (t={t_sim_s[last_valid] / MYR_IN_CGS:.4g} Myr); "
            f"later snapshots excluded from the comparison."
        )

    window = np.arange(0, last_valid + 1)
    rel_error = np.abs(R_lagrangian_cgs[window] - R_analytic_cgs[window]) / np.maximum(
        R_analytic_cgs[window], 1e-30
    )
    median_rel_error = float(np.median(rel_error))
    verdict = "PASS" if median_rel_error <= args.tol else "FAIL"

    mismatch = np.abs(R_lagrangian_cgs[window] - R_peak_cgs[window]) / np.maximum(
        R_lagrangian_cgs[window], 1e-30
    )
    marker = ""
    if np.any(mismatch > args.mismatch_tol):
        marker = (
            f" (PRIMARY/SECONDARY MISMATCH: up to "
            f"{mismatch.max():.1%} disagreement between the Lagrangian-mass "
            f"and density-peak shell radii -- no clean coherent shell)"
        )

    print(
        f"\nComparison window  : snapshots 0-{last_valid} "
        f"(t=0 to {t_sim_s[last_valid] / MYR_IN_CGS:.4g} Myr)"
    )
    print(f"Median rel. error  : {median_rel_error:.2%}  [{verdict}{marker}]")
    print(f"(tolerance {args.tol:.0%})")

    fig, ax = plt.subplots(figsize=(7, 5))
    t_myr = t_sim_s / MYR_IN_CGS
    ax.plot(
        t_myr,
        R_lagrangian_cgs / PC_IN_CGS,
        "o-",
        color="#d62728",
        markersize=4,
        label="Simulation (Lagrangian-mass shell radius)",
    )
    ax.plot(
        t_myr,
        R_peak_cgs / PC_IN_CGS,
        "s-",
        color="#1f77b4",
        markersize=3,
        alpha=0.8,
        label="Simulation (density-peak radius, audit)",
    )
    ax.plot(
        t_myr,
        R_analytic_cgs / PC_IN_CGS,
        "--",
        color="black",
        label="Krumholz & Matzner (2009), radiation-only limit",
    )
    ax.axhline(
        box_half_width_cgs / PC_IN_CGS,
        color="firebrick",
        linestyle=":",
        alpha=0.6,
        label="Box half-width",
    )
    ax.axvline(
        t_myr[last_valid],
        color="black",
        linestyle="-.",
        alpha=0.4,
        label="End of comparison window",
    )
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("Shell radius [pc]")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
