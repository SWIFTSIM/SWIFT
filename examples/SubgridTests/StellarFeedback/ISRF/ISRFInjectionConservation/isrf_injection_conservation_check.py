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
Checks one ISRF star-feedback injection pass against its own exact identities.

Zero-metallicity leg (default)
------------------------------
At Z=0 the receiver-side dust extinction factor is exp(0)=1 exactly (no
dust), so summing the injected energy over every gas particle the star's
kernel reaches must reproduce Delta_t * L_band exactly, up to
floating-point precision:

    sum_j(u_PE_j * m_j) == Delta_t * L_PE_star
    sum_j(u_LW_j  * m_j) == Delta_t * L_LW_star

This is an exact conservation identity (the injection formula's own
`sum_j weight_j = 1` normalization), not an approximate physics comparison
-- unlike Tier 1's propagation-profile fit, no loose tolerance is
expected here.

Nonzero-metallicity leg (--dusty)
---------------------------------
With dust the extinction factor is no longer 1 and the kernel-weighted
mean of exp(-tau) does not factor out of the sum, so the identity above
has no direct analogue. Two weaker but still exact statements replace it.

G2, the primary gate. The two bands differ only through their dust
cross-section per hydrogen nucleon, so on the same particle the kernel
weight, Delta_t, the particle mass, Z_j, the dust-to-gas ratio, mu_H, m_H
and every unit conversion cancel in the band ratio:

    u_PE_j / u_LW_j = (L_PE / L_LW) * exp(+(1 - sigma_PE/sigma_LW) tau_LW_j)

    R_j = ln(u_PE_j / u_LW_j) - ln(L_PE / L_LW) - 0.4 * tau_LW_j,   target 0

with sigma_PE / sigma_LW = 9.0e-22 / 1.5e-21 = 0.6 exactly (radiation.h).
The gate is max_j |R_j| over the illuminated set.

G1, the secondary gate. G2 cancels everything common to the two bands and
therefore cannot see a common-mode error (a wrong Delta_t, a broken
`sum_j weight_j = 1`, a uniform factor on both bands). Since the weights
are non-negative and sum to 1, the global sum is rigorously bracketed:

    min_j exp(-tau_b_j) <= sum_j m_j u_b_j / (Delta_t L_b) <= max_j exp(-tau_b_j)

What this tests, and what it does not
-------------------------------------
It tests that the extinction is WIRED per receiving particle with the
correct band-dependent cross-section ratio: that the call is made at all,
that the two bands are not swapped, that the column uses the receiver's
own h and rho with the configured path length and kernel_gamma, that the
comoving-to-physical scaling is applied once, and that the metallicity is
read from the array the code reads.

It does NOT test whether the column model, the cross-sections or the
dust-to-gas convention are physically right: the reconstruction reuses
all three. The formula itself is unit-tested in
tests/testRadiationISRFFormula.c.

It is also weak on one wiring error: applying the extinction source-side,
once for the whole star, instead of receiver-side per particle. In a
uniform glass box the column barely varies between particles, so the two
give nearly the same answer here. A density-gradient fixture would be
needed to separate them.

Delta_t (the star's own feedback-timestep at the checked snapshot) is read
from swift's own step-table log, not assumed to equal TimeIntegration:dt_max.
"""

import argparse
import glob
import os
import sys
from typing import Optional, Tuple

import h5py
import numpy as np

# Mirrors of src/feedback/GEAR/radiation.h. The reconstruction reuses the
# code's own constants on purpose: this check gates the wiring, not the
# physical values (see the module docstring).
SIGMA_D_PE_CGS = 9e-22
SIGMA_D_LW_CGS = 1.5e-21
MU_H = 1.4
HYDROGEN_MASS_CGS = 1.6726219e-24
GRACKLE_SOLAR_METAL_FRACTION = 0.01295
GRACKLE_DEFAULT_DUST_TO_GAS_RATIO = 0.009387

# src/kernel_hydro.h, Wendland C2 in 3D.
KERNEL_GAMMA_DEFAULT = 1.936492

# Half an ulp of float32, the unit of the reported error budget.
FLOAT32_ULP = 2.0**-24


def parse_options() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--snapshot",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for snapshots to consider (default: %(default)s). "
        "The last one (highest time) is checked.",
    )
    parser.add_argument(
        "--log",
        default="output.log",
        help="Run log to read the star's own feedback Delta_t from "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-5,
        help="Max allowed relative error on the summed energy (default: "
        "%(default)s). This is an exact identity, so the only expected "
        "discrepancy is float32 rounding in the kernel-weight normalization, "
        "~1e-7. The default leaves two decades of margin above that floor. "
        "With --dusty the same value is the relative slack allowed on the "
        "G1 bracket edges.",
    )
    parser.add_argument(
        "--dusty",
        action="store_true",
        help="Check the nonzero-metallicity extinction identities (G1 and G2) "
        "instead of the Z=0 summed-energy identity. Without it a nonzero "
        "metallicity is an error, not a mode.",
    )
    parser.add_argument(
        "--dust-tol",
        type=float,
        default=1e-4,
        help="--dusty only: max allowed max_j |R_j| on the band-ratio gate G2 "
        "(default: %(default)s). The residual is pinned to 0 by an identity "
        "with no free parameters, so the only expected discrepancy is float32 "
        "rounding, measured at ~3e-7 at this example's defaults. The default "
        "leaves two decades of margin above that floor, matching how --tol was "
        "set. Pass a negative value to report the residual without gating on "
        "it, which is how the floor was measured.",
    )
    parser.add_argument(
        "--kernel-gamma",
        type=float,
        default=KERNEL_GAMMA_DEFAULT,
        help="--dusty only: the binary's kernel_gamma, which sets the "
        "extinction column together with the path length (default: "
        "%(default)s, Wendland C2 in 3D). A run built with another kernel "
        "needs the matching value from src/kernel_hydro.h.",
    )
    parser.add_argument(
        "--extinction-path",
        type=float,
        default=None,
        help="--dusty only: GEARFeedback:ISRF_extinction_path in kernel "
        "support radii (constant_kernel_path = its own float, "
        "kernel_diameter = 2, kernel_radius = 1). Read from the run's "
        "used_parameters.yml when omitted.",
    )
    parser.add_argument(
        "--dust-to-gas-ratio",
        type=float,
        default=GRACKLE_DEFAULT_DUST_TO_GAS_RATIO,
        help="--dusty only: the resolved Grackle "
        "chemistry_data.local_dust_to_gas_ratio of the run (default: "
        "%(default)s, Grackle's own compiled default). The run log prints it "
        "to three digits only, so it is an argument here and the log value is "
        "used to cross-check it.",
    )
    return parser.parse_args()


def read_delta_t(log_path: str, time: float) -> float:
    """Read the star's feedback Delta_t from swift's own step-table log, by
    matching the row whose Time column equals the checked snapshot's time.
    Never assume Delta_t == TimeIntegration:dt_max: individual
    time-stepping can settle on a smaller step (see README)."""
    with open(log_path, "r") as f:
        for line in f:
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                t = float(fields[1])
                dt = float(fields[4])
            except ValueError:
                continue
            if abs(t - time) < 1e-4 * abs(time) or (time == 0.0 and t == 0.0):
                return dt
    raise RuntimeError(
        f"Could not find a step in {log_path!r} matching snapshot time {time:.6e}."
    )


def read_log_dust_to_gas_ratio(log_path: str) -> Optional[float]:
    """Return the local_dust_to_gas_ratio the run log reports, or None.

    Grackle's value is printed with three significant digits, which is too
    coarse to build the optical depth from; it is read only to confirm that
    the value the reconstruction uses is the one the run actually had.
    """
    key = "grackle_chemistry_data.local_dust_to_gas_ratio"
    try:
        with open(log_path, "r") as f:
            for line in f:
                if key in line:
                    return float(line.rsplit("=", 1)[1])
    except (OSError, ValueError):
        return None
    return None


def read_extinction_path(pattern: str, given: Optional[float]) -> float:
    """Return the extinction path in kernel radii, from the argument or the
    run's used_parameters.yml."""
    if given is not None:
        return given
    import yaml

    directory = os.path.dirname(os.path.dirname(sorted(glob.glob(pattern))[0]))
    with open(os.path.join(directory, "used_parameters.yml")) as handle:
        used = yaml.safe_load(handle)["GEARFeedback"]
    if "ISRF_extinction_path" not in used:
        # A run archived before the key existed recorded no value at all, and
        # the path then in force was two kernel support radii.
        return 2.0
    name = used["ISRF_extinction_path"]
    if name == "constant_kernel_path":
        return float(used["ISRF_extinction_path_in_kernel_radii"])
    paths = {"kernel_diameter": 2.0, "kernel_radius": 1.0}
    if name in paths:
        return paths[name]
    if name in ("pair_separation", "temperature_capped_jeans"):
        raise RuntimeError(
            f"GEARFeedback:ISRF_extinction_path {name!r} does not build the "
            "column from a multiple of the kernel support radius, so the "
            "closed form this check mirrors does not apply to it. Rerun the "
            "fixture with a constant_kernel_path mechanism, or extend this "
            "check to that mechanism's own length."
        )
    raise RuntimeError(f"Unknown GEARFeedback:ISRF_extinction_path {name!r}")


def physical(dataset: h5py.Dataset, a: float, unit_cgs: float) -> np.ndarray:
    """Return a snapshot dataset in physical CGS."""
    exponent = float(np.atleast_1d(dataset.attrs["a-scale exponent"])[0])
    return dataset[:].astype(np.float64) * a**exponent * unit_cgs


def optical_depths(
    h_cgs: np.ndarray,
    rho_cgs: np.ndarray,
    Z: np.ndarray,
    path_in_kernel_radii: float,
    kernel_gamma: float,
    local_dust_to_gas_ratio: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the PE and LW dust optical depths of each gas particle.

    Mirrors radiation_get_part_ISRF_extinction_factors and the chain below
    it in src/feedback/GEAR/radiation_isrf.c, entirely in physical CGS: the
    comoving-to-physical a^-2 of the code is exactly what converting the
    comoving column to a physical one does here.

    Parameters
    ----------
    h_cgs
        Physical smoothing lengths, cm.
    rho_cgs
        Physical mass densities, g cm^-3.
    Z
        Metal mass fractions, as the code reads them (the smoothed array,
        cast to float32).
    path_in_kernel_radii
        GEARFeedback:ISRF_extinction_path, in kernel support radii.
    kernel_gamma
        The binary's kernel_gamma.
    local_dust_to_gas_ratio
        The run's resolved Grackle chemistry_data.local_dust_to_gas_ratio.

    Returns
    -------
    tau_pe, tau_lw
        Dimensionless optical depths of each particle, in the two bands.
    """
    sigma_gas = path_in_kernel_radii * kernel_gamma * h_cgs * rho_cgs
    d_relative = (
        np.maximum(Z, 0.0)
        / GRACKLE_SOLAR_METAL_FRACTION
        * (local_dust_to_gas_ratio / GRACKLE_DEFAULT_DUST_TO_GAS_RATIO)
    )
    prefactor = d_relative / (MU_H * HYDROGEN_MASS_CGS) * sigma_gas
    return SIGMA_D_PE_CGS * prefactor, SIGMA_D_LW_CGS * prefactor


def gate(label: str, worst: float, bar: float) -> bool:
    """Print a pass/fail line, failing closed on a non-finite value."""
    ok = bool(np.isfinite(worst)) and worst <= bar
    print(f"  {'PASS' if ok else 'FAIL'}: {label}: {worst:.3e} vs bar {bar:.3e}")
    return ok


def report_sum(band: str, lhs: float, rhs: float, tol: float) -> bool:
    """Print and gate one band's summed-energy identity at Z=0."""
    rel_err = abs(lhs - rhs) / rhs
    status = "PASS" if np.isfinite(rel_err) and rel_err < tol else "FAIL"
    print(
        f"{band}: sum(u*m)={lhs:.10e}, Delta_t*L={rhs:.10e}, "
        f"rel_err={rel_err:.3e} -> {status}"
    )
    return status == "PASS"


def check_dusty(
    opt: argparse.Namespace,
    mass: np.ndarray,
    u_pe: np.ndarray,
    u_lw: np.ndarray,
    tau_pe: np.ndarray,
    tau_lw: np.ndarray,
    L_PE: float,
    L_LW: float,
    Delta_t: float,
) -> bool:
    """Gate the two nonzero-metallicity extinction identities, G2 then G1."""
    if not (np.isfinite(L_PE) and np.isfinite(L_LW)) or L_PE <= 0.0 or L_LW <= 0.0:
        print(
            f"  FAIL: star luminosities must be finite and positive, got {L_PE}, {L_LW}"
        )
        return False

    # Explicitly, before anything selects on `u > 0`: a NaN compares false
    # against every bound, so it would otherwise drop out of the illuminated
    # set unnoticed rather than fail the gate.
    for name, array in (
        ("masses", mass),
        ("u_PE", u_pe),
        ("u_LW", u_lw),
        ("tau_PE", tau_pe),
        ("tau_LW", tau_lw),
    ):
        bad = int(np.sum(~np.isfinite(array)))
        if bad:
            print(f"  FAIL: {bad} non-finite {name} value(s) in the snapshot")
            return False

    lit_pe = u_pe > 0.0
    lit_lw = u_lw > 0.0
    if int(np.sum(lit_pe)) != int(np.sum(lit_lw)):
        print(
            f"  FAIL: the two bands illuminate different particle counts, "
            f"{int(np.sum(lit_pe))} PE against {int(np.sum(lit_lw))} LW"
        )
        return False
    lit = lit_pe
    if not np.any(lit):
        print("  FAIL: no illuminated gas particle")
        return False

    ratio = SIGMA_D_PE_CGS / SIGMA_D_LW_CGS
    residual = (
        np.log(u_pe[lit] / u_lw[lit])
        - np.log(L_PE / L_LW)
        - (1.0 - ratio) * tau_lw[lit]
    )
    worst = float(np.max(np.abs(residual)))

    print(f"  tau_LW: min {np.min(tau_lw[lit]):.4f}, max {np.max(tau_lw[lit]):.4f}")
    print(f"  tau_PE: min {np.min(tau_pe[lit]):.4f}, max {np.max(tau_pe[lit]):.4f}")
    print(
        f"  band-ratio signal (1 - sigma_PE/sigma_LW) tau_LW: "
        f"{(1.0 - ratio) * float(np.max(tau_lw[lit])):.4f}"
    )

    # Float32 error budget of the reconstruction, each term about half an
    # ulp: the two specific-energy stores, the two expf calls, and the
    # float32 chain that builds tau (h, rho, and about four roundings in the
    # opacity), the last three scaling with the signal itself.
    signal = abs((1.0 - ratio) * float(np.max(tau_lw[lit])))
    budget = 2.0 * FLOAT32_ULP + 2.0 * FLOAT32_ULP + 6.0 * FLOAT32_ULP * signal
    print(f"  float32 budget of |R_j|: {budget:.2e}")

    if opt.dust_tol < 0.0:
        print(
            f"  REPORT (no bar given): band-ratio residual (G2), "
            f"max_j |R_j| = {worst:.3e}"
        )
        ok = bool(np.isfinite(worst))
        if not ok:
            print("  FAIL: the band-ratio residual is not finite")
    else:
        ok = gate("band-ratio residual (G2), max_j |R_j|", worst, opt.dust_tol)

    for band, u, tau, L in (("PE", u_pe, tau_pe, L_PE), ("LW", u_lw, tau_lw, L_LW)):
        measured = float(np.sum(mass * u)) / (Delta_t * L)
        low = float(np.min(np.exp(-tau[lit])))
        high = float(np.max(np.exp(-tau[lit])))
        inside = (
            bool(np.isfinite(measured))
            and np.isfinite(low)
            and np.isfinite(high)
            and high > 0.0
            and low * (1.0 - opt.tol) <= measured <= high * (1.0 + opt.tol)
        )
        print(
            f"  {'PASS' if inside else 'FAIL'}: weighted-mean extinction bracket "
            f"(G1) {band}: {measured:.8f} in [{low:.8f}, {high:.8f}]"
        )
        ok &= inside
    return ok


def main() -> int:
    opt = parse_options()
    files = sorted(glob.glob(opt.snapshot))
    if not files:
        raise RuntimeError(f"No snapshots match {opt.snapshot}")

    # Use the last snapshot. ISRF_propagation is off, so the field is
    # reset and fully re-injected on every star-feedback pass (see README):
    # any snapshot after at least one pass is a clean, self-contained check
    # of that pass's own injection, independent of how many earlier passes
    # ran.
    snap_path = files[-1]
    with h5py.File(snap_path, "r") as f:
        header = f["/Header"].attrs
        time = float(np.asarray(header["Time"]).flat[0])
        a = (
            float(np.atleast_1d(header["Scale-factor"])[0])
            if "Scale-factor" in header
            else 1.0
        )
        units = f["/Units"].attrs
        length_cgs = float(np.atleast_1d(units["Unit length in cgs (U_L)"])[0])
        mass_cgs = float(np.atleast_1d(units["Unit mass in cgs (U_M)"])[0])
        gas = f["/PartType0"]
        mass = gas["Masses"][:].astype(np.float64)
        u_pe = gas["PESpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
        # The unsmoothed array, for the Z=0 leg's own precondition only.
        Z = gas["MetalMassFractions"][:, -1]
        if opt.dusty:
            # The code reads the SMOOTHED metal mass fraction
            # (chemistry_get_total_metal_mass_fraction_for_cooling), as a
            # float32. The two arrays coincide only at Z=0.
            Z_used = gas["SmoothedMetalMassFractions"][:, -1].astype(np.float32)
            h_cgs = physical(gas["SmoothingLengths"], a, length_cgs)
            rho_cgs = physical(gas["Densities"], a, mass_cgs / length_cgs**3)

        star = f["/PartType4"]
        L_PE = float(star["PELuminosities"][0])
        L_LW = float(star["LWLuminosities"][0])

    if not opt.dusty and np.any(Z != 0.0):
        raise RuntimeError(
            "This check requires GEARChemistry:initial_metallicity=0 (exact "
            "extinction=1); found nonzero metallicity in the snapshot. Pass "
            "--dusty to check the extinction identities instead."
        )

    Delta_t = read_delta_t(opt.log, time)

    n_illuminated = int(np.sum(u_pe > 0))
    print(f"Snapshot: {snap_path} (t={time:.6e})")
    print(f"Delta_t (from {opt.log}): {Delta_t:.6e}")
    print(f"Star: L_PE={L_PE:.10e}, L_LW={L_LW:.10e}")
    print(f"Illuminated gas particles: {n_illuminated}")

    if not opt.dusty:
        ok_pe = report_sum("PE", float(np.sum(u_pe * mass)), Delta_t * L_PE, opt.tol)
        ok_lw = report_sum("LW", float(np.sum(u_lw * mass)), Delta_t * L_LW, opt.tol)
        return 0 if (ok_pe and ok_lw) else 1

    if not np.any(Z_used > 0.0):
        raise RuntimeError(
            "--dusty needs a nonzero metallicity; the snapshot's smoothed "
            "metal mass fraction is everywhere zero."
        )

    path = read_extinction_path(opt.snapshot, opt.extinction_path)
    logged = read_log_dust_to_gas_ratio(opt.log)
    if logged is not None and abs(logged - opt.dust_to_gas_ratio) > 5e-3 * logged:
        raise RuntimeError(
            f"--dust-to-gas-ratio {opt.dust_to_gas_ratio} disagrees with the "
            f"{logged} the run log reports."
        )
    print(
        f"Extinction column: path {path:g} kernel radii, kernel_gamma "
        f"{opt.kernel_gamma:g}, local_dust_to_gas_ratio {opt.dust_to_gas_ratio:g}"
        f"{'' if logged is None else f' (log reports {logged:g})'}"
    )
    print(f"Smoothed Z: min {np.min(Z_used):.6e}, max {np.max(Z_used):.6e}")

    tau_pe, tau_lw = optical_depths(
        h_cgs, rho_cgs, Z_used, path, opt.kernel_gamma, opt.dust_to_gas_ratio
    )
    ok = check_dusty(opt, mass, u_pe, u_lw, tau_pe, tau_lw, L_PE, L_LW, Delta_t)
    print("RESULT: PASS" if ok else "RESULT: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
