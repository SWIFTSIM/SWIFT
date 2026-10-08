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
"""Check that the total ISRF field decays away once the star dies.

Tracks the box-total mass-weighted PE and LW specific energies,
``sum(mass*PESpecificEnergies)``/``sum(mass*LWSpecificEnergies)``, across
every snapshot. The star's death time is not read off any snapshot field
(a dead star's last-computed luminosity is not necessarily reset to zero,
so it is not a reliable death marker on its own); it is instead computed
analytically from the star's own birth mass and metallicity, using GEAR's
Poirier main-sequence-lifetime fit read directly from the run's own yields
table (mirrors ``lifetime_get_log_lifetime_from_mass``,
``src/feedback/GEAR/lifetime.h``).

Reports three things; only (b) and (c) are pass/fail gates, (a) is
informational:

(a) INFORMATIONAL. The post-death e-folding time of each band, fitted
    from the decay itself, against the analytic absorption timescale
    ``1/(c_hyp*kappa)`` (``c_hyp`` and ``kappa`` reconstructed from the
    run's own smoothing length/time-step/density/metallicity, mirroring
    ``radiation_isrf.c``'s own formulas). The two need not match exactly,
    since the shipped configuration also runs the negativity-triggered
    artificial-dissipation term, which reshapes the field spatially but
    should not change the total-energy decay rate; a large mismatch is
    worth a look, not an automatic failure.
(b) GATED. The residual fraction at the final snapshot relative to the
    pre-death level (the median box-total energy over the last three
    pre-death snapshots), gated at 1% provided at least 5 e-folds of decay
    time have elapsed (5 e-folds of pure exponential decay already leaves
    under 1% residual, so requiring both jointly rules out a stalled-but-
    still-below-1%-by-coincidence read).
(c) GATED. Whether any particle's PE or LW specific energy goes negative
    after death, excluding snapshots within ``--negativity-settle-myr`` of
    death (the switch-off transient: the field is still relaxing towards
    its post-death form there, and is reported separately, not scored).
"""

import argparse
import glob
import sys

import h5py
import matplotlib
import numpy as np
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Same constants as the sibling ISRF*/checks (src/feedback/GEAR/radiation.h).
SIGMA_D_PE_CGS = 9e-22
SIGMA_D_LW_CGS = 1.5e-21
MU_H = 1.4
M_H_CGS = 1.6726219e-24
# Grackle's SolarMetalFractionByMass default, src/feedback/GEAR/radiation.h
GRACKLE_SOLAR_Z = 0.01295
SPEED_OF_LIGHT_KM_S = 2.99792458e5
# src/physical_constants_cgs.h: const_year_cgs
YEAR_CGS = 3.15569251e7
# src/physical_constants_cgs.h: const_solar_mass_cgs (non-GADGET2 default)
SOLAR_MASS_CGS = 1.98841e33


def parse_options() -> argparse.Namespace:
    """Parse command-line options.

    Returns
    -------
    argparse.Namespace
        The parsed options.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-s",
        "--snapshot",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for snapshots to consider (default: %(default)s)",
    )
    parser.add_argument(
        "--yields-table",
        default="PopII_parsec_spectral.hdf5",
        help="Yields table to read the lifetime fit from (default: %(default)s)",
    )
    parser.add_argument(
        "--timesteps-log",
        default="timesteps.txt",
        help="SWIFT's timesteps log, for the bulk gas time-step used to "
        "reconstruct c_hyp (default: %(default)s)",
    )
    parser.add_argument(
        "--c-hyp-margin",
        type=float,
        default=0.5,
        help="GEARFeedback:ISRF_c_hyp_margin used by the run (default: "
        "%(default)s).",
    )
    parser.add_argument(
        "--residual-tol",
        type=float,
        default=0.01,
        help="Max allowed residual fraction (band energy at the final "
        "snapshot over the pre-death level) for a pass (default: "
        "%(default)s).",
    )
    parser.add_argument(
        "--min-efolds",
        type=float,
        default=5.0,
        help="Minimum number of measured e-folds of post-death decay time "
        "required for the residual check to count as a pass (default: "
        "%(default)s).",
    )
    parser.add_argument(
        "--negativity-settle-myr",
        type=float,
        default=None,
        help="Post-death window (Myr) excluded from the negativity gate as "
        "a switch-off transient; snapshots inside it are still printed, "
        "tagged separately, and not scored (default: one fine output-list "
        "interval from --output-list, or 0 if that cannot be read).",
    )
    parser.add_argument(
        "--output-list",
        default="output_list_isrf_source_turn_off.txt",
        help="SWIFT Snapshots:output_list file, used only to size the "
        "default --negativity-settle-myr (default: %(default)s).",
    )
    parser.add_argument(
        "--used-parameters",
        default="used_parameters.yml",
        help="SWIFT's used_parameters.yml, a fallback source for dt_max if "
        "--timesteps-log is unavailable (default: %(default)s).",
    )
    parser.add_argument(
        "--output",
        default="isrf_source_turn_off_check.png",
        help="Output plot filename.",
    )
    return parser.parse_args()


def modal_bulk_dt(path: str, n_total: int) -> float:
    """Return the bulk gas time-step from SWIFT's timesteps log.

    The modal interval between successive whole-box updates
    (``Updates == n_total``); mirrors the sibling ISRFCausalReach/
    ISRFDissipation checks' own ``modal_bulk_dt``.

    Parameters
    ----------
    path : str
        Path to SWIFT's ``timesteps.txt``.
    n_total : int
        Number of gas particles (``Updates`` counts gas only).

    Returns
    -------
    float
        The modal (most common) gas time-step, in internal units.
    """
    times = []
    updates = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            step = int(parts[0])
            if step == 0:
                continue
            times.append(float(parts[1]))
            updates.append(int(parts[7]))
    times_arr, updates_arr = np.array(times), np.array(updates)
    full = times_arr[updates_arr == n_total]
    if len(full) < 2:
        vals, counts = np.unique(times_arr[1:] - times_arr[:-1], return_counts=True)
        return float(vals[np.argmax(counts)])
    diffs = np.round(np.diff(full), 12)
    vals, counts = np.unique(diffs, return_counts=True)
    return float(vals[np.argmax(counts)])


def get_dt_max(timesteps_log: str, used_parameters: str, n_total: int) -> float:
    """Return the run's dt_max, from the timesteps log or the parameters.

    Tries :func:`modal_bulk_dt` on ``timesteps_log`` first (the step
    actually taken, which need not equal ``TimeIntegration:dt_max`` if the
    limiter or another timestep constraint binds); falls back to the
    literal ``TimeIntegration:dt_max`` in ``used_parameters`` if the log is
    unavailable.

    Parameters
    ----------
    timesteps_log : str
        Path to SWIFT's ``timesteps.txt``.
    used_parameters : str
        Path to SWIFT's ``used_parameters.yml``.
    n_total : int
        Number of gas particles, passed to :func:`modal_bulk_dt`.

    Returns
    -------
    float
        dt_max, in internal units, or ``numpy.nan`` if neither source is
        readable.
    """
    try:
        return modal_bulk_dt(timesteps_log, n_total)
    except (OSError, ValueError, IndexError):
        pass
    try:
        with open(used_parameters) as f:
            used = yaml.safe_load(f)
        return float(used["TimeIntegration"]["dt_max"])
    except (OSError, KeyError, TypeError, ValueError):
        return float("nan")


def default_negativity_settle_myr(output_list: str, unit_time_myr: float) -> float:
    """Return the default switch-off transient window, in Myr.

    One fine output-list interval (the smallest spacing between two
    entries of ``output_list``); 0 if that file cannot be read.

    Parameters
    ----------
    output_list : str
        Path to SWIFT's ``Snapshots:output_list`` file.
    unit_time_myr : float
        Internal-time-to-Myr conversion factor for this run.

    Returns
    -------
    float
        The default settle window, in Myr.
    """
    try:
        with open(output_list) as f:
            vals = np.array(
                [float(line) for line in f if line.strip() and not line.startswith("#")]
            )
    except (OSError, ValueError):
        return 0.0
    if len(vals) < 2:
        return 0.0
    return float(np.min(np.diff(np.sort(vals)))) * unit_time_myr


def load_snapshot(path: str) -> dict:
    """Read the fields this check needs from one snapshot.

    Parameters
    ----------
    path : str
        Path to the snapshot file.

    Returns
    -------
    dict
        Time, unit conversion factors, and the gas/star fields used below.
    """
    with h5py.File(path, "r") as f:
        header = f["/Header"]
        time = float(np.asarray(header.attrs["Time"]).flat[0])
        units = f["/Units"]
        unit_length_cgs = float(
            np.asarray(units.attrs["Unit length in cgs (U_L)"]).flat[0]
        )
        unit_mass_cgs = float(np.asarray(units.attrs["Unit mass in cgs (U_M)"]).flat[0])
        unit_time_cgs = float(np.asarray(units.attrs["Unit time in cgs (U_t)"]).flat[0])

        gas = f["/PartType0"]
        h = gas["SmoothingLengths"][:].astype(np.float64)
        rho = gas["Densities"][:].astype(np.float64)
        mass = gas["Masses"][:].astype(np.float64)
        u_pe = gas["PESpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
        Z = gas["MetalMassFractions"][:, -1].astype(np.float64)

        star = f["/PartType4"]
        star_mass = float(star["Masses"][0])
        star_birth_time = float(star["BirthTimes"][0])

    return dict(
        time=time,
        unit_length_cgs=unit_length_cgs,
        unit_mass_cgs=unit_mass_cgs,
        unit_time_cgs=unit_time_cgs,
        h=h,
        rho=rho,
        mass=mass,
        u_pe=u_pe,
        u_lw=u_lw,
        Z=Z,
        star_mass=star_mass,
        star_birth_time=star_birth_time,
    )


def lifetime_coefficients(
    yields_table: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read GEAR's Poirier main-sequence lifetime fit coefficients.

    Mirrors ``lifetime_read_from_tables``/``lifetime_init``
    (``src/feedback/GEAR/lifetime.h``): the table stores the fit as a flat
    9-element array (quadratic, linear, constant terms in metallicity, in
    that order); the constant term's log10(lifetime) is stored in years and
    is shifted here to Myr, exactly as ``lifetime_init`` does.

    Parameters
    ----------
    yields_table : str
        Path to the yields table (e.g. ``PopII_parsec_spectral.hdf5``).

    Returns
    -------
    tuple of numpy.ndarray
        ``(quadratic, linear, constant)``, each a 3-element array of the
        metallicity-polynomial coefficients.
    """
    with h5py.File(yields_table, "r") as f:
        coeff = np.asarray(f["Data/LifeTimes/coeff_z"][:], dtype=np.float64)
    quadratic, linear, constant = coeff[0].copy(), coeff[1].copy(), coeff[2].copy()
    constant[-1] -= 6.0  # yr -> Myr
    return quadratic, linear, constant


def lifetime_myr(
    mass_msun: float, Z: float, coeffs: tuple[np.ndarray, np.ndarray, np.ndarray]
) -> float:
    """Evaluate the Poirier lifetime fit.

    Parameters
    ----------
    mass_msun : float
        Star mass in solar masses.
    Z : float
        Total metal mass fraction.
    coeffs : tuple of numpy.ndarray
        Output of :func:`lifetime_coefficients`.

    Returns
    -------
    float
        Main-sequence lifetime, in Myr.
    """
    quadratic, linear, constant = coeffs
    q = (quadratic[0] * Z + quadratic[1]) * Z + quadratic[2]
    l = (linear[0] * Z + linear[1]) * Z + linear[2]
    c = (constant[0] * Z + constant[1]) * Z + constant[2]
    log_m = np.log10(mass_msun)
    return 10.0 ** ((q * log_m + l) * log_m + c)


def kappa_eff_mass_opacity_cgs(Z: float, sigma_d_cgs: float) -> float:
    """Band dust mass opacity (area/mass, cgs).

    Mirrors ``radiation_get_dust_mass_opacity`` (``radiation_isrf.c``),
    with ``local_dust_to_gas_ratio`` left at Grackle's own default (this
    example's ``params.yml`` leaves it unset), so ``D_relative`` reduces to
    ``Z/GRACKLE_SOLAR_Z``.

    Parameters
    ----------
    Z : float
        Total metal mass fraction. Must be > 0: this dust-opacity model is
        undefined at Z=0 (see the caller's guard).
    sigma_d_cgs : float
        Band cross-section per hydrogen nucleon, cgs.

    Returns
    -------
    float
        Mass opacity, cgs (area/mass).
    """
    D_relative = Z / GRACKLE_SOLAR_Z
    return sigma_d_cgs * D_relative / (MU_H * M_H_CGS)


def fit_efold_time(t_rel: np.ndarray, energy: np.ndarray, plateau: float) -> float:
    """Fit an exponential decay's e-folding time.

    Fits ``ln(energy/plateau)`` linearly against ``t_rel`` over the points
    where the signal is still well above the numerical floor
    (``energy > 1e-6*plateau``), and returns ``-1/slope``.

    Parameters
    ----------
    t_rel : numpy.ndarray
        Time since death.
    energy : numpy.ndarray
        Band energy at each time.
    plateau : float
        Reference value the decay is measured relative to (e.g. the
        pre-death level).

    Returns
    -------
    float
        Fitted e-folding time, or ``numpy.nan`` if fewer than 3 points are
        above the floor.
    """
    ratio = energy / plateau
    mask = ratio > 1e-6
    if np.sum(mask) < 3:
        return float("nan")
    slope, _ = np.polyfit(t_rel[mask], np.log(ratio[mask]), 1)
    return -1.0 / slope


def main() -> None:
    """Run the source-turn-off check and print/plot its report."""
    opt = parse_options()
    files = sorted(glob.glob(opt.snapshot))
    if not files:
        raise RuntimeError(f"No snapshots match {opt.snapshot}")

    snaps = [load_snapshot(fn) for fn in files]
    order = np.argsort([s["time"] for s in snaps])
    snaps = [snaps[i] for i in order]

    times = np.array([s["time"] for s in snaps])

    # A NaN specific energy is neither < 0 nor summable, so it would otherwise
    # read as a false PASS/NaN in the two gated legs below rather than a FAIL.
    n_nonfinite_pe = np.array([int(np.sum(~np.isfinite(s["u_pe"]))) for s in snaps])
    n_nonfinite_lw = np.array([int(np.sum(~np.isfinite(s["u_lw"]))) for s in snaps])
    bad = (n_nonfinite_pe > 0) | (n_nonfinite_lw > 0)
    if np.any(bad):
        i = int(np.argmax(bad))
        raise RuntimeError(
            f"Non-finite specific energy at snapshot {i} (t={times[i]:.6e}): "
            f"{n_nonfinite_pe[i]} PE / {n_nonfinite_lw[i]} LW particles out of "
            f"{len(snaps[i]['u_pe'])} are NaN/inf. This is a corrupted radiation "
            "field (a simulation-side defect), not a check-script threshold "
            "issue: fix the run, not this gate."
        )

    E_pe = np.array([np.sum(s["mass"] * s["u_pe"]) for s in snaps])
    E_lw = np.array([np.sum(s["mass"] * s["u_lw"]) for s in snaps])
    min_u_pe = np.array([np.min(s["u_pe"]) for s in snaps])
    min_u_lw = np.array([np.min(s["u_lw"]) for s in snaps])

    unit_length_cgs = snaps[0]["unit_length_cgs"]
    unit_mass_cgs = snaps[0]["unit_mass_cgs"]
    unit_time_cgs = snaps[0]["unit_time_cgs"]
    unit_time_myr = unit_time_cgs / (1e6 * YEAR_CGS)

    # Birth mass: the star's mass at the first snapshot. Valid here because
    # stellar winds and supernovae are both off before death (params.yml),
    # so the star's mass cannot have changed yet.
    star_mass_msun = snaps[0]["star_mass"] * unit_mass_cgs / SOLAR_MASS_CGS
    star_birth_time = snaps[0]["star_birth_time"]
    Z = float(np.median(snaps[0]["Z"]))

    if Z <= 0.0:
        raise ValueError(
            "Z = 0 (GEARChemistry:initial_metallicity must be > 0 for this "
            "example): the dust mass opacity kappa_eff is "
            "exactly zero at Z = 0, so the analytic decay timescale "
            "1/(c_hyp*kappa) this check relies on (both the informational "
            "e-fold comparison and the gated residual check) is undefined."
        )

    coeffs = lifetime_coefficients(opt.yields_table)
    t_death_myr = lifetime_myr(star_mass_msun, Z, coeffs)
    t_death = star_birth_time + t_death_myr / unit_time_myr

    print(f"Star mass: {star_mass_msun:.3f} Msun, Z: {Z:.4e}")
    print(f"Main-sequence lifetime: {t_death_myr:.4f} Myr")
    print(
        f"Death time: {t_death:.6e} (internal) = " f"{t_death * unit_time_myr:.4f} Myr"
    )

    pre_mask = times < t_death
    post_mask = times > t_death
    if np.sum(pre_mask) < 3 or np.sum(post_mask) < 3:
        raise RuntimeError(
            "Need at least 3 snapshots on each side of the death time; got "
            f"{np.sum(pre_mask)} before and {np.sum(post_mask)} after."
        )

    pre_indices = np.where(pre_mask)[0]

    # Pre-death level (b)/(c) are each measured against: the box-total
    # energy for the residual gate, the per-particle specific energy for
    # the negativity report. Both are the median over the last three
    # pre-death snapshots, kept as two distinct quantities rather than one
    # rescaled by mass, so each name matches what it actually computes.
    pre_death_level_pe = float(np.median(E_pe[pre_mask][-3:]))
    pre_death_level_lw = float(np.median(E_lw[pre_mask][-3:]))
    pre_death_median_u_pe = float(
        np.median([np.median(snaps[i]["u_pe"]) for i in pre_indices[-3:]])
    )
    pre_death_median_u_lw = float(
        np.median([np.median(snaps[i]["u_lw"]) for i in pre_indices[-3:]])
    )

    t_rel_post = times[post_mask] - t_death
    tau_pe_measured_myr = (
        fit_efold_time(t_rel_post, E_pe[post_mask], pre_death_level_pe) * unit_time_myr
    )
    tau_lw_measured_myr = (
        fit_efold_time(t_rel_post, E_lw[post_mask], pre_death_level_lw) * unit_time_myr
    )

    # Analytic expected e-folding time, 1/(c_hyp*kappa), reconstructed from
    # the run's own bulk gas time-step/smoothing length/density (mirrors
    # radiation_isrf.c's own c_hyp/kappa formulas). INFORMATIONAL only (see
    # module docstring): dt_bulk is the step actually taken (timesteps.txt),
    # deliberately not the get_dt_max() value used below for the negativity
    # report, since the two need not agree (see get_dt_max's docstring).
    n_gas = len(snaps[0]["mass"])
    dt_bulk = modal_bulk_dt(opt.timesteps_log, n_gas)
    last_pre = pre_indices[-1]
    h_med = float(np.median(snaps[last_pre]["h"]))
    rho_med_internal = float(np.median(snaps[last_pre]["rho"]))
    c_hyp = min(opt.c_hyp_margin * h_med / dt_bulk, SPEED_OF_LIGHT_KM_S)
    c_hyp_cgs = c_hyp * 1e5
    rho_cgs = rho_med_internal * unit_mass_cgs / unit_length_cgs**3

    tau_pe_expected_s = 1.0 / (
        c_hyp_cgs * kappa_eff_mass_opacity_cgs(Z, SIGMA_D_PE_CGS) * rho_cgs
    )
    tau_lw_expected_s = 1.0 / (
        c_hyp_cgs * kappa_eff_mass_opacity_cgs(Z, SIGMA_D_LW_CGS) * rho_cgs
    )
    tau_pe_expected_myr = tau_pe_expected_s / (1e6 * YEAR_CGS)
    tau_lw_expected_myr = tau_lw_expected_s / (1e6 * YEAR_CGS)

    print()
    print("--- (a) e-folding time, measured vs. analytic (INFORMATIONAL) ---")
    print(f"c_hyp (reconstructed): {c_hyp:.4f} km/s")
    print(
        f"PE: measured e-fold={tau_pe_measured_myr:.4f} Myr, "
        f"expected={tau_pe_expected_myr:.4f} Myr, "
        f"ratio={tau_pe_measured_myr / tau_pe_expected_myr:.3f}"
    )
    print(
        f"LW:  measured e-fold={tau_lw_measured_myr:.4f} Myr, "
        f"expected={tau_lw_expected_myr:.4f} Myr, "
        f"ratio={tau_lw_measured_myr / tau_lw_expected_myr:.3f}"
    )

    residual_pe = E_pe[-1] / pre_death_level_pe
    residual_lw = E_lw[-1] / pre_death_level_lw
    total_decay_time_myr = (times[-1] - t_death) * unit_time_myr
    n_efolds_pe = total_decay_time_myr / tau_pe_measured_myr
    n_efolds_lw = total_decay_time_myr / tau_lw_measured_myr

    print()
    print("--- (b) residual-fraction check (GATED) ---")
    print(
        f"PE: residual at time_end = {residual_pe:.3e} "
        f"({n_efolds_pe:.2f} e-folds elapsed)"
    )
    print(
        f"LW:  residual at time_end = {residual_lw:.3e} "
        f"({n_efolds_lw:.2f} e-folds elapsed)"
    )

    pass_residual = (
        residual_pe < opt.residual_tol
        and residual_lw < opt.residual_tol
        and n_efolds_pe >= opt.min_efolds
        and n_efolds_lw >= opt.min_efolds
    )

    # (c) Negativity check (GATED), excluding the switch-off transient.
    settle_myr = opt.negativity_settle_myr
    if settle_myr is None:
        settle_myr = default_negativity_settle_myr(opt.output_list, unit_time_myr)
    settle_internal = settle_myr / unit_time_myr

    dt_max_report = get_dt_max(opt.timesteps_log, opt.used_parameters, n_gas)

    print()
    print("--- (c) negativity check (GATED) ---")
    print(
        f"Switch-off transient window: {settle_myr:.6f} Myr after death; "
        "snapshots inside it are excluded from the gate below and reported "
        "separately."
    )

    post_indices = np.where(post_mask)[0]
    t_rel_internal_post = times[post_indices] - t_death
    is_transient = t_rel_internal_post <= settle_internal
    is_offender = (min_u_pe[post_indices] < 0) | (min_u_lw[post_indices] < 0)

    if not np.any(is_offender):
        print("No offending post-death snapshot (all specific energies stayed >= 0).")
    for k, i in enumerate(post_indices):
        if not is_offender[k]:
            continue
        t_rel_myr = t_rel_internal_post[k] * unit_time_myr
        if np.isfinite(dt_max_report) and dt_max_report > 0:
            dt_units_str = f", {t_rel_internal_post[k] / dt_max_report:.2f} dt_max"
        else:
            dt_units_str = " (dt_max unavailable)"
        n_neg_pe_i = int(np.sum(snaps[i]["u_pe"] < 0))
        n_neg_lw_i = int(np.sum(snaps[i]["u_lw"] < 0))
        tag = " [switch-off transient, not scored]" if is_transient[k] else ""
        print(
            f"  snapshot {i}: t-t_death = {t_rel_myr:.5f} Myr{dt_units_str}, "
            f"negative PE particles = {n_neg_pe_i}, "
            f"negative LW particles = {n_neg_lw_i}, "
            f"min(u_pe)/pre-death median = "
            f"{min_u_pe[i] / pre_death_median_u_pe:.3e}, "
            f"min(u_lw)/pre-death median = "
            f"{min_u_lw[i] / pre_death_median_u_lw:.3e}{tag}"
        )

    scored = ~is_transient
    n_negative_scored_pe = int(np.sum((min_u_pe[post_indices] < 0) & scored))
    n_negative_scored_lw = int(np.sum((min_u_lw[post_indices] < 0) & scored))
    pass_negativity = n_negative_scored_pe == 0 and n_negative_scored_lw == 0

    print(
        f"Scored snapshots with a negative PE specific energy: "
        f"{n_negative_scored_pe}/{int(np.sum(scored))}"
    )
    print(
        f"Scored snapshots with a negative LW specific energy: "
        f"{n_negative_scored_lw}/{int(np.sum(scored))}"
    )

    print()
    print(f"Residual-decay check (gated): {'PASS' if pass_residual else 'FAIL'}")
    print(f"Negativity check (gated): {'PASS' if pass_negativity else 'FAIL'}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(times, E_pe, "o-", color="C0", label="PE (total)")
    ax.semilogy(times, E_lw, "s-", color="C1", label="LW (total)")
    ax.axvline(t_death, color="k", linestyle="--", label="star death")
    ax.set_xlabel("time (internal units)")
    ax.set_ylabel(r"$\sum_i m_i u_i$ (internal units)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(opt.output, dpi=150)
    print(f"Plot saved to {opt.output}")

    if not (pass_residual and pass_negativity):
        sys.exit(1)


if __name__ == "__main__":
    main()
