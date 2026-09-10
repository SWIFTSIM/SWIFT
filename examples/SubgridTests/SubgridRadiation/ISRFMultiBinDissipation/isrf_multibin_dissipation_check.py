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
Per-run, report-only metrics for ISRFMultiBinDissipation. Reconstructs the
force loop's `dissipation_u_FUV`/`dissipation_u_LW` accumulators from
snapshot fields alone (the accumulator itself is not a snapshot field), the
same technique validated in the root-cause log's `diss_dipole.py` and the
Phase-2 adjudication's `adjudicate_wake.py`. Never exits nonzero; writes
`multibin_metrics.json` for `multibin_compare.py` to gate on.
"""

import argparse
import glob
import json
import os
import sys

import h5py
import numpy as np
import yaml
from scipy.spatial import cKDTree

# src/feedback/GEAR/radiation.h
SIGMA_D_FUV_CGS = 9e-22
SIGMA_D_LW_CGS = 1.5e-21
MU_H = 1.4
M_H_CGS = 1.6726219e-24
GRACKLE_SOLAR_Z = 0.01295  # Grackle's SolarMetalFractionByMass default
GAMMA_3D = 1.936492  # Wendland C2, 3D: src/kernel_hydro.h kernel_gamma
SPEED_OF_LIGHT_KM_S = 2.99792458e5

NEG_WEIGHT_VOID_THRESHOLD = 0.10
DRHO_VOID_THRESHOLD = 0.10
DX_GAS_VOID_THRESHOLD_H = 0.10
RHO_MATCH_THRESHOLD = 0.15
BIN_DT_DISAGREEMENT_THRESHOLD = 0.20

# Wendland C2, 3D coefficients (0 < x < 1 branch); src/kernel_hydro.h.
WENDLAND_C2_3D_COEFFS = [4.0, -15.0, 20.0, -10.0, 0.0, 1.0]
KERNEL_CONSTANT_3D = 21.0 * np.pi**-1 / 2.0


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-s", "--snapshot", default="snap/snapshot_*.hdf5")
    parser.add_argument("--timesteps-log", default="timesteps.txt")
    parser.add_argument("--ic-json", default="multibin_ic.json")
    parser.add_argument("--used-parameters", default="used_parameters.yml")
    parser.add_argument("--c-hyp-margin", type=float, default=0.5)
    parser.add_argument(
        "--c-hyp-pin",
        type=float,
        default=0.0,
        help="LW_FUV_c_hyp_pin_for_debugging used by the run; > 0 makes "
        "c_hyp a global constant, matching radiation_isrf.c's own closure.",
    )
    parser.add_argument("--r-cut-h", type=float, default=6.0)
    parser.add_argument("--json-out", default="multibin_metrics.json")
    parser.add_argument("--output", default="isrf_multibin_dissipation_check.png")
    return parser.parse_args()


# -----------------------------------------------------------------------------
# Family idioms, copied verbatim per this family's own convention (no shared
# module): modal_bulk_dt / load_snapshot / radial_distance.
# -----------------------------------------------------------------------------
def modal_bulk_dt(path, n_total):
    """The coarse phase's own real period: the modal spacing between
    successive whole-box (Updates == n_total) rows of timesteps.txt."""
    times, updates = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            step = int(parts[0])
            if step == 0:
                continue
            times.append(float(parts[1]))
            updates.append(int(parts[7]))  # Updates (gas); col 5,6 are min/max Time-bins
    times, updates = np.array(times), np.array(updates)
    full = times[updates == n_total]
    if len(full) < 2:
        vals, counts = np.unique(times[1:] - times[:-1], return_counts=True)
        return float(vals[np.argmax(counts)])
    diffs = np.round(np.diff(full), 12)
    vals, counts = np.unique(diffs, return_counts=True)
    return float(vals[np.argmax(counts)])


def modal_time_step_column(path):
    """Mode of the raw `Time-step` column itself (excluding step 0): the
    finest active bin's own step size, since it fires every row."""
    steps = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if int(parts[0]) == 0:
                continue
            steps.append(float(parts[4]))  # Time-step (col 4; col 3 is Redshift)
    steps = np.round(np.array(steps), 12)
    vals, counts = np.unique(steps, return_counts=True)
    return float(vals[np.argmax(counts)]), vals, counts


def load_snapshot(path):
    with h5py.File(path, "r") as f:
        header = f["/Header"]
        time = float(np.asarray(header.attrs["Time"]).flat[0])
        boxsize = np.asarray(header.attrs["BoxSize"], dtype=float).flatten()[0]
        units = f["/Units"]
        unit_length_cgs = float(
            np.asarray(units.attrs["Unit length in cgs (U_L)"]).flat[0]
        )
        unit_mass_cgs = float(np.asarray(units.attrs["Unit mass in cgs (U_M)"]).flat[0])
        gas = f["/PartType0"]
        pos = gas["Coordinates"][:, :]
        h = gas["SmoothingLengths"][:].astype(np.float64)
        mass = gas["Masses"][:].astype(np.float64)
        rho = gas["Densities"][:].astype(np.float64)
        ids = gas["ParticleIDs"][:]
        u_fuv = gas["FUVSpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
        alpha_fuv = gas["FUVArtificialDissipationCoefficients"][:].astype(np.float64)
        alpha_lw = gas["LWArtificialDissipationCoefficients"][:].astype(np.float64)
        Z = gas["MetalMassFractions"][:, -1].astype(np.float64)
        star = f["/PartType4"]
        star_pos = star["Coordinates"][:, :]
    return dict(
        time=time,
        boxsize=boxsize,
        unit_length_cgs=unit_length_cgs,
        unit_mass_cgs=unit_mass_cgs,
        pos=pos,
        h=h,
        mass=mass,
        rho=rho,
        ids=ids,
        u_fuv=u_fuv,
        u_lw=u_lw,
        alpha_fuv=alpha_fuv,
        alpha_lw=alpha_lw,
        Z=Z,
        star_pos=star_pos,
    )


def radial_distance(pos, ref, boxsize):
    dx = pos - ref
    dx -= boxsize * np.round(dx / boxsize)
    return dx, np.sqrt(np.sum(dx**2, axis=1))


def classify_hot(mass, m_cold, m_hot):
    """Nearest-match phase label from Masses (constraint C10): np.isclose's
    default atol (1e-8) swamps these internal-unit masses (~1e-11 to
    1e-13), so an absolute-tolerance comparison is unusable here."""
    return np.abs(mass - m_hot) < np.abs(mass - m_cold)


# -----------------------------------------------------------------------------
# Kernel and dissipation reconstruction.
# -----------------------------------------------------------------------------
def kernel_deval_np(u):
    """Vectorized port of src/kernel_hydro.h `kernel_deval`, Wendland C2 3D.
    `u` is r/h (each particle's own h), matching the C call site."""
    kernel_gamma_inv = 1.0 / GAMMA_3D
    x = np.asarray(u, dtype=np.float64) * kernel_gamma_inv
    c = WENDLAND_C2_3D_COEFFS
    w = c[0] * x + c[1]
    dw = np.full_like(x, c[0])
    for k in range(2, 6):
        dw = dw * x + w
        w = x * w + c[k]
    w = np.maximum(w, 0.0)
    dw = np.minimum(dw, 0.0)
    mask = x < 1.0
    w = np.where(mask, w, 0.0)
    dw = np.where(mask, dw, 0.0)
    kgid = kernel_gamma_inv**3
    kgidp1 = kernel_gamma_inv**4
    W = w * KERNEL_CONSTANT_3D * kgid
    dWdx = dw * KERNEL_CONSTANT_3D * kgidp1
    return W, dWdx


def w_dr(r, h):
    """`wi_dr`/`wj_dr` in radiation_propagation_iact.h: dW/dx(r/h) * h^-4
    (3D pow_dimension_plus_one)."""
    _, dWdx = kernel_deval_np(r / h)
    return dWdx * (1.0 / h) ** 4


def kappa_eff_mass_opacity_cgs(Z, sigma_d_cgs):
    """Mirrors radiation_get_dust_mass_opacity's cgs stage; local_dust_to_gas_ratio
    left unset (params.yml), so D_relative = Z/Z_grackle_sun."""
    D_relative = Z / GRACKLE_SOLAR_Z
    return sigma_d_cgs * D_relative / (MU_H * M_H_CGS)


def kappa_internal(Z, rho_internal, unit_length_cgs, unit_mass_cgs, sigma_d_cgs):
    """radiation_get_part_linear_absorption_rate: kappa_eff(Z) [internal
    area/mass] * rho_p [internal]; internal units of 1/length."""
    kappa_eff_cgs = kappa_eff_mass_opacity_cgs(Z, sigma_d_cgs)
    kappa_eff_internal = kappa_eff_cgs * unit_mass_cgs / unit_length_cgs**2
    return kappa_eff_internal * rho_internal


def reconstruct_pairs(pos_w, h_w, boxsize):
    """All (i, j) index pairs (local to the window) with
    r < GAMMA_3D*max(h_i, h_j) -- the SPH neighbour criterion both the
    symmetric and non-symmetric dissipation interactions use."""
    tree = cKDTree(pos_w, boxsize=boxsize)
    h_max = float(np.max(h_w)) if len(h_w) else 0.0
    pairs = np.array(
        sorted(tree.query_pairs(r=GAMMA_3D * h_max)), dtype=np.int64
    )
    if pairs.size == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64), None, None
    i_idx, j_idx = pairs[:, 0], pairs[:, 1]
    dx = pos_w[i_idx] - pos_w[j_idx]
    dx -= boxsize * np.round(dx / boxsize)
    r = np.sqrt(np.sum(dx**2, axis=1))
    keep = r < GAMMA_3D * np.maximum(h_w[i_idx], h_w[j_idx])
    return i_idx[keep], j_idx[keep], dx[keep], r[keep]


def dissipation_accumulator(i_idx, j_idx, r, h_w, mass_w, rho_w, c_hyp_w, alpha_w, u_w):
    """`radiation_dissipation_force_accumulate_band`, Stage 1 (no
    reconstruction branch: RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION is
    compiled out by default)."""
    n = len(mass_w)
    diss_u = np.zeros(n)
    if len(i_idx) == 0:
        return diss_u
    wi = w_dr(r, h_w[i_idx])
    wj = w_dr(r, h_w[j_idx])
    Wbar = 0.5 * (wi + wj)
    d_ij = rho_w[i_idx] * u_w[i_idx] - rho_w[j_idx] * u_w[j_idx]
    v_sig = np.maximum(alpha_w[i_idx], alpha_w[j_idx]) * np.minimum(
        c_hyp_w[i_idx], c_hyp_w[j_idx]
    )
    Psi = v_sig * d_ij * Wbar / (rho_w[i_idx] * rho_w[j_idx])
    np.add.at(diss_u, i_idx, mass_w[j_idx] * Psi)
    np.add.at(diss_u, j_idx, -mass_w[i_idx] * Psi)
    return diss_u


def wendland_c2_eval(q):
    w = np.zeros_like(q)
    sel = q < 1.0
    t = 1.0 - q[sel]
    w[sel] = (21.0 / (2.0 * np.pi)) * t**4 * (4.0 * q[sel] + 1.0)
    return w


def kernel_mean_at_point(point, pos, h, field, boxsize):
    """Kernel-weighted mean of `field` at `point`, same idiom as
    shear_compare.py's sph_interpolate but for one point."""
    dx, r = radial_distance(pos, point, boxsize)
    H = GAMMA_3D * h
    sel = r < H
    if not np.any(sel):
        return np.nan
    q = r[sel] / H[sel]
    w = wendland_c2_eval(q) / H[sel] ** 3
    wsum = w.sum()
    if wsum <= 0:
        return np.nan
    return float(np.sum(w * field[sel]) / wsum)


# -----------------------------------------------------------------------------
def main():
    opt = parse_options()
    files = sorted(glob.glob(opt.snapshot))
    if not files:
        raise RuntimeError(f"No snapshots match {opt.snapshot}")

    with open(opt.ic_json) as f:
        ic = json.load(f)
    L = ic["L"]
    h_median = ic["h_median"]
    m_cold = ic["m_cold"]
    m_hot = ic["m_hot"]
    bin_delta = ic["bin_delta"]
    variant = ic["variant"]
    star_positions = ic["star_positions"]  # name -> [x, y, z]
    star_names = ic["star_names"]

    used_params = {}
    if os.path.exists(opt.used_parameters):
        with open(opt.used_parameters) as f:
            used_params = yaml.safe_load(f)

    with h5py.File(files[0], "r") as f:
        n_gas = f["/PartType0/Coordinates"].shape[0]

    # --- bin reconstruction from timesteps.txt ---
    # dt_cold_realized: the coarse phase's own real period, from the modal
    # spacing between full-box (Updates == n_gas) syncs -- empirically
    # robust, cross-checked against the analytic floor below.
    #
    # dt_hot_realized is NOT read from the raw `Time-step` column's mode:
    # that column reports the current GLOBAL finest active tick, which is
    # permanently contaminated once any interface-layer particle lands on a
    # bin finer than the hot bulk's own (confirmed in S0: the column reads
    # 1.5625e-6 from step 28 on, while step 1's clean, uncontaminated
    # Updates == n_hot row shows the hot bulk's own dt is 3.125e-6, exactly
    # the floored-analytic prediction). Instead dt_hot_realized is derived
    # from dt_cold_realized and the intended bin separation: SWIFT floors
    # dt to the same power-of-two timeline on both phases, and that floor
    # logic is independently confirmed exact on the cold side (see the
    # dt_cold_realized vs dt_cold_analytic check below), so applying it via
    # the ratio rather than re-detecting it empirically is the more robust
    # route. bin_delta is dN, a bin-LEVEL count: each level is a factor of
    # 2 in dt, so the ratio is 2**bin_delta -- NOT 4**bin_delta (that would
    # conflate the temperature ratio with the dt ratio; the design doc's
    # own Sec 3.3.1 states the dt ratio "is exactly 4 by construction" at
    # bin_delta=2, i.e. 2**bin_delta, contradicting its Sec 3.4 cross-check
    # formula, which is corrected here).
    dt_cold_realized = modal_bulk_dt(opt.timesteps_log, n_gas)
    ratio_expected = 2.0**bin_delta
    dt_hot_realized = dt_cold_realized / ratio_expected if variant == "twophase" else dt_cold_realized
    ratio_realized = ratio_expected if variant == "twophase" else 1.0
    bin_span = float(np.log2(ratio_realized)) if ratio_realized > 0 else 0.0

    # Diagnostic-only global-finest-tick histogram (for bin_report): NOT
    # used for the dt_hot_realized derivation above.
    _global_finest_mode, dt_vals, dt_counts = modal_time_step_column(opt.timesteps_log)
    order = np.argsort(-dt_counts)
    bin_report = [
        {"dt": float(dt_vals[k]), "count": int(dt_counts[k])} for k in order[:8]
    ]

    # The one genuine empirical cross-check: dt_cold_realized must be the
    # nearest-power-of-two-floor of dt_cold_analytic. A naive relative
    # disagreement bar is unusable here -- power-of-two flooring alone
    # spreads the realized value uniformly in log space across a factor of
    # 2, so up to ~50% disagreement is ordinary quantization, not a defect
    # (S0 measured 22.5%, passing this floor check cleanly). A dt off by an
    # integer bin (a real defect) fails this test instead.
    dt_cold_analytic = ic["dt_cold_analytic"]
    dt_analytic_disagreement = (
        abs(dt_cold_realized - dt_cold_analytic) / dt_cold_analytic
    )
    dt_cold_is_floor_of_analytic = (
        dt_cold_realized <= dt_cold_analytic < 2.0 * dt_cold_realized
    )
    invalid_bin_dt = (variant == "twophase") and (not dt_cold_is_floor_of_analytic)

    v_over_c_hyp_realized = ic["v_star"] / (
        opt.c_hyp_margin * h_median / dt_cold_realized
    )

    print(f"n_gas={n_gas}  variant={variant}  bin_delta={bin_delta}")
    print(f"dt_hot (derived: dt_cold/2**bin_delta) = {dt_hot_realized:.6e}")
    print(f"dt_cold (modal_bulk_dt, empirical)     = {dt_cold_realized:.6e}")
    print(f"ratio (by construction)                = {ratio_realized:.4f}  (expected {ratio_expected:.1f})")
    print(f"bin span (log2 ratio)                  = {bin_span:.3f}")
    print(f"dt_cold vs analytic disagreement       = {dt_analytic_disagreement*100:.2f}%  "
          f"(floor check: {dt_cold_realized:.4e} <= {dt_cold_analytic:.4e} < {2*dt_cold_realized:.4e} "
          f"-> {'PASS' if dt_cold_is_floor_of_analytic else 'FAIL'})")
    print(f"global finest-tick mode (diagnostic only, NOT dt_hot) = {_global_finest_mode:.6e}")
    print(f"realized v_star/c_hyp_cold              = {v_over_c_hyp_realized:.4f}")
    if invalid_bin_dt:
        print("INVALID: dt_cold_realized is not the nearest power-of-two floor of dt_cold_analytic")

    # --- per-snapshot bookkeeping ---
    snaps = [load_snapshot(fn) for fn in files]
    first, last = snaps[0], snaps[-1]

    # SWIFT's particle array order is NOT stable across snapshots (confirmed:
    # cell-based spatial sorting reorders the array, especially with two
    # very different particle masses in the same box); any first-vs-last
    # comparison must align by ParticleIDs, not by raw array index, or it
    # silently compares unrelated particles. ids are a dense 0..N-1 range
    # (makeIC.py), so a direct index map suffices: `align[k]` is last's
    # array position of the same physical particle as first["..."][k].
    last_slot_by_id = np.empty(n_gas, dtype=np.int64)
    last_slot_by_id[last["ids"]] = np.arange(n_gas)
    align = last_slot_by_id[first["ids"]]

    # M5 drho, dx_gas
    is_hot0 = classify_hot(first["mass"], m_cold, m_hot) if variant == "twophase" else np.zeros(n_gas, bool)
    drho_by_phase = {}
    last_rho_aligned = last["rho"][align]
    for label, sel in (("cold", ~is_hot0), ("hot", is_hot0)):
        if not np.any(sel):
            continue
        rho0 = np.median(first["rho"][sel])
        rho1 = np.median(last_rho_aligned[sel])
        drho_by_phase[label] = abs(rho1 - rho0) / rho0 if rho0 > 0 else np.nan
    drho_max = max(drho_by_phase.values()) if drho_by_phase else 0.0

    last_pos_aligned = last["pos"][align]
    dx, dr = radial_distance(last_pos_aligned, first["pos"], L)
    dx_gas_h = float(np.median(dr)) / h_median

    void_drho = drho_max > DRHO_VOID_THRESHOLD
    void_dx_gas = dx_gas_h > DX_GAS_VOID_THRESHOLD_H
    print(f"drho (max over phases)  = {drho_max*100:.3f}%  -> {'VOID' if void_drho else 'ok'}")
    print(f"dx_gas (median, /h)     = {dx_gas_h:.4f}  -> {'VOID' if void_dx_gas else 'ok'}")

    # rho_match, A vs B (twophase only)
    rho_match = None
    ab_matched = True
    if variant == "twophase" and "A" in star_positions and "B" in star_positions:
        rho_A = kernel_mean_at_point(
            np.array(star_positions["A"]), last["pos"], last["h"], last["rho"], L
        )
        rho_B = kernel_mean_at_point(
            np.array(star_positions["B"]), last["pos"], last["h"], last["rho"], L
        )
        rho_match = abs(rho_A / rho_B - 1.0) if rho_B else np.nan
        ab_matched = rho_match <= RHO_MATCH_THRESHOLD
        print(f"rho_match |rho_A/rho_B-1| = {rho_match*100:.3f}%  -> {'matched' if ab_matched else 'A/B not matched'}")

    # --- M1/M2/M3 per star, per band, per snapshot ---
    R_cut = opt.r_cut_h
    C_hyp = opt.c_hyp_margin
    c_pin = opt.c_hyp_pin

    per_star_series = {name: {"FUV": [], "LW": []} for name in star_names}

    # M4 conservation trace (global sum(m*u)), per band.
    conservation = {"FUV": [], "LW": []}

    for snap in snaps:
        pos, h, mass, rho = snap["pos"], snap["h"], snap["mass"], snap["rho"]
        is_hot = classify_hot(mass, m_cold, m_hot) if variant == "twophase" else np.zeros(n_gas, bool)
        dt_i = np.where(is_hot, dt_hot_realized, dt_cold_realized)
        if c_pin > 0.0:
            c_hyp = np.full(n_gas, c_pin)
        else:
            c_hyp = np.minimum(C_hyp * h / dt_i, SPEED_OF_LIGHT_KM_S)

        kappa = {}
        for band, sigma in (("FUV", SIGMA_D_FUV_CGS), ("LW", SIGMA_D_LW_CGS)):
            kappa[band] = kappa_internal(
                snap["Z"], rho, snap["unit_length_cgs"], snap["unit_mass_cgs"], sigma
            )

        for band, u_field, alpha_field in (
            ("FUV", "u_fuv", "alpha_fuv"),
            ("LW", "u_lw", "alpha_lw"),
        ):
            u = snap[u_field]
            alpha = snap[alpha_field]
            conservation[band].append((snap["time"], float(np.sum(mass * u))))

            for name in star_names:
                star_pos = np.array(star_positions[name])
                dxs, rs = radial_distance(pos, star_pos, L)
                window = rs < (R_cut + 2.0) * h_median
                if window.sum() < 10:
                    continue
                w_idx = np.where(window)[0]
                pos_w = pos[w_idx]
                h_w = h[w_idx]
                mass_w = mass[w_idx]
                rho_w = rho[w_idx]
                c_hyp_w = c_hyp[w_idx]
                alpha_w = alpha[w_idx]
                u_w = u[w_idx]

                i_loc, j_loc, dx_pair, r_pair = reconstruct_pairs(pos_w, h_w, L)
                diss_u_w = dissipation_accumulator(
                    i_loc, j_loc, r_pair, h_w, mass_w, rho_w, c_hyp_w, alpha_w, u_w
                )

                # Tight subset (R_cut*h) for the actual moment.
                dx_tight, r_tight = radial_distance(pos_w, star_pos, L)
                tight = r_tight < R_cut * h_median
                if tight.sum() < 5:
                    continue
                xi = dx_tight[tight]
                mt = mass_w[tight]
                ut = u_w[tight]
                dt_u = diss_u_w[tight]
                c_hyp_t = c_hyp_w[tight]
                kappa_t = kappa[band][w_idx][tight]

                neg_share = (
                    np.sum(mt * np.maximum(-ut, 0.0)) / np.sum(mt * np.abs(ut))
                    if np.sum(mt * np.abs(ut)) > 0
                    else 0.0
                )
                is_void = neg_share > NEG_WEIGHT_VOID_THRESHOLD

                M0_pos = np.sum(mt * np.maximum(ut, 0.0))
                M0_signed = np.sum(mt * ut)
                c_hyp_S = float(np.median(c_hyp_t))
                kappa_S = float(np.median(kappa_t))
                tau = 1.0 / (c_hyp_S * kappa_S) if (c_hyp_S > 0 and kappa_S > 0) else np.nan

                d_vec = np.full(3, np.nan)
                if M0_pos > 0 and np.isfinite(tau):
                    for k in range(3):
                        d_vec[k] = -tau * np.sum(mt * xi[:, k] * dt_u) / M0_pos

                per_star_series[name][band].append(
                    dict(
                        time=snap["time"],
                        d_x=float(d_vec[0]) if np.isfinite(d_vec[0]) else None,
                        d_y=float(d_vec[1]) if np.isfinite(d_vec[1]) else None,
                        d_z=float(d_vec[2]) if np.isfinite(d_vec[2]) else None,
                        d_x_h=float(d_vec[0] / h_median) if np.isfinite(d_vec[0]) else None,
                        d_y_h=float(d_vec[1] / h_median) if np.isfinite(d_vec[1]) else None,
                        d_z_h=float(d_vec[2] / h_median) if np.isfinite(d_vec[2]) else None,
                        d_total_h=(
                            float(np.linalg.norm(d_vec) / h_median)
                            if np.all(np.isfinite(d_vec))
                            else None
                        ),
                        M0_pos=float(M0_pos),
                        M0_signed=float(M0_signed),
                        neg_weight_share=float(neg_share),
                        void=bool(is_void),
                        tau=float(tau) if np.isfinite(tau) else None,
                        c_hyp_S=c_hyp_S,
                        kappa_S=kappa_S,
                        n_particles=int(tight.sum()),
                    )
                )

    # --- M3 mechanism diagnostic at star A, last snapshot ---
    m3 = {}
    if "A" in star_positions:
        snap = last
        pos, h, mass, rho = snap["pos"], snap["h"], snap["mass"], snap["rho"]
        is_hot = classify_hot(mass, m_cold, m_hot) if variant == "twophase" else np.zeros(n_gas, bool)
        dt_i = np.where(is_hot, dt_hot_realized, dt_cold_realized)
        c_hyp = (
            np.full(n_gas, c_pin)
            if c_pin > 0.0
            else np.minimum(C_hyp * h / dt_i, SPEED_OF_LIGHT_KM_S)
        )
        star_pos = np.array(star_positions["A"])
        _, rs = radial_distance(pos, star_pos, L)
        window = rs < (R_cut + 2.0) * h_median
        w_idx = np.where(window)[0]
        pos_w, h_w, mass_w = pos[w_idx], h[w_idx], mass[w_idx]
        for band, alpha_field in (("FUV", "alpha_fuv"), ("LW", "alpha_lw")):
            alpha_w = snap[alpha_field][w_idx]
            c_hyp_w = c_hyp[w_idx]
            i_loc, j_loc, dx_pair, r_pair = reconstruct_pairs(pos_w, h_w, L)
            if len(i_loc) == 0:
                m3[band] = None
                continue
            _, r_tight = radial_distance(pos_w, star_pos, L)
            both_tight = (r_tight[i_loc] < R_cut * h_median) & (
                r_tight[j_loc] < R_cut * h_median
            )
            i_t, j_t = i_loc[both_tight], j_loc[both_tight]
            is_hot_w = is_hot[w_idx]
            cross = is_hot_w[i_t] != is_hot_w[j_t]
            E = np.abs(alpha_w[i_t] - alpha_w[j_t])
            Dbar = np.maximum(alpha_w[i_t], alpha_w[j_t]) * np.minimum(
                c_hyp_w[i_t], c_hyp_w[j_t]
            )
            E_cross = float(np.mean(E[cross])) if np.any(cross) else None
            E_same = float(np.mean(E[~cross])) if np.any(~cross) else None
            Dbar_cross = float(np.mean(Dbar[cross])) if np.any(cross) else None
            Dbar_same = float(np.mean(Dbar[~cross])) if np.any(~cross) else None
            m3[band] = dict(
                E_cross=E_cross,
                E_same=E_same,
                E_cross_over_same=(E_cross / E_same if E_cross and E_same else None),
                Dbar_cross=Dbar_cross,
                Dbar_same=Dbar_same,
                Dbar_cross_over_same=(
                    Dbar_cross / Dbar_same if Dbar_cross and Dbar_same else None
                ),
                n_cross_pairs=int(np.sum(cross)),
                n_same_pairs=int(np.sum(~cross)),
            )

    # --- M4 conservation drift ---
    # The field is not seeded in the IC (Sec 4.1), so the very first
    # snapshot is always exactly zero for every run; using it as the
    # fractional-change baseline would divide by zero. Use the first
    # NONZERO snapshot instead (the first snapshot after injection begins).
    m4 = {}
    for band in ("FUV", "LW"):
        series = conservation[band]
        nonzero = [(t, v) for t, v in series if v != 0.0]
        t_last, c_last = series[-1]
        if nonzero:
            t0, c0 = nonzero[0]
            frac = (c_last - c0) / c0 if c0 != 0 else np.nan
        else:
            t0, c0, frac = None, 0.0, np.nan
        m4[band] = dict(
            first=c0,
            first_time=t0,
            last=c_last,
            last_time=t_last,
            fractional_change=frac,
            full_trace=series,
        )

    print("\n--- M1/M2 summary (last snapshot) ---")
    for name in star_names:
        for band in ("FUV", "LW"):
            series = per_star_series[name][band]
            if not series:
                continue
            last_pt = series[-1]
            print(
                f"star {name} {band}: d_x/h={last_pt['d_x_h']}  d_y/h={last_pt['d_y_h']}  "
                f"|d|/h={last_pt['d_total_h']}  neg_weight={last_pt['neg_weight_share']*100:.2f}%  "
                f"void={last_pt['void']}  n={last_pt['n_particles']}"
            )

    print("\n--- M3 (star A, last snapshot) ---")
    for band, res in m3.items():
        print(f"{band}: {res}")

    print("\n--- M4 conservation ---")
    for band, res in m4.items():
        print(f"{band}: first={res['first']:.6e}  last={res['last']:.6e}  "
              f"frac_change={res['fractional_change']*100:.4f}%")

    status = "INVALID" if invalid_bin_dt else ("VOID" if (void_drho or void_dx_gas) else "ok")
    print(f"\nRun status (report-only): {status}")

    out = dict(
        variant=variant,
        bin_delta=bin_delta,
        n_gas=n_gas,
        dt_hot_realized=dt_hot_realized,
        dt_cold_realized=dt_cold_realized,
        ratio_realized=ratio_realized,
        ratio_expected=ratio_expected,
        bin_span=bin_span,
        dt_analytic_disagreement=dt_analytic_disagreement,
        invalid_bin_dt=bool(invalid_bin_dt),
        v_over_c_hyp_realized=v_over_c_hyp_realized,
        bin_report=bin_report,
        drho_by_phase=drho_by_phase,
        drho_max=drho_max,
        void_drho=bool(void_drho),
        dx_gas_h=dx_gas_h,
        void_dx_gas=bool(void_dx_gas),
        rho_match=rho_match,
        ab_matched=ab_matched,
        per_star_series=per_star_series,
        m3_mechanism_diagnostic=m3,
        m4_conservation=m4,
        used_parameters=used_params,
        ic=ic,
        c_hyp_margin=opt.c_hyp_margin,
        c_hyp_pin=opt.c_hyp_pin,
        r_cut_h=R_cut,
        status=status,
    )
    with open(opt.json_out, "w") as f:
        json.dump(out, f, indent=2, default=lambda o: None)
    print(f"\n{opt.json_out} saved.")

    # Never exits nonzero: report-only, per the family's run.sh convention.
    sys.exit(0)


if __name__ == "__main__":
    main()
