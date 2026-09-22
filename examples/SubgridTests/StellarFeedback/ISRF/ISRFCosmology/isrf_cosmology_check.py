################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
"""Check the ISRF module against closed-form solutions, with or without cosmology.

Every snapshot field is converted to physical CGS with its own
``a-scale exponent`` attribute. Elapsed proper time is ``Time`` minus the run's
start time. ``a0`` is the scale factor at the start (1 without cosmology) and
``H`` the Hubble rate from the snapshot's own ``Cosmology`` group,

    H(a) = H0 sqrt(Omega_r a^-4 + Omega_m a^-3 + Omega_k a^-2 + Omega_lambda) .

Integrals over time are done in ln a, ``dt = d ln a / H``, between the scale
factors SWIFT wrote, so the reference uses SWIFT's own a(t).

free_field
    No star, no dust (kappa = 0), uniform seeded field, so the flux divergence
    vanishes. The module's energy equation reduces to ``du/dt = -(c_hyp/c) H
    u`` for the mass-specific field: the reduced-speed-of-light method is
    only correct if EVERY rate is dilated by the same c_hyp/c factor, and the
    Hubble term is no exception (see radiation_isrf.c's
    radiation_end_force_propagation, fixed for this). ``c_hyp`` in this
    fixture is NOT pinned (c_hyp_pin is 0, the run.sh default for
    ``free_field``), so it is a per-particle, per-step quantity
    (ISRF_c_hyp_margin*h/dt, clamped at c), not a single constant, and the
    exact solution is the integral

        ln[Q(t)/Q(0)] = -(1/c) int_0^t c_hyp(t') H(t') dt'                  (A1)

    with Q the ledger the run's own propagation scheme conserves. Which one
    that is depends on GEARFeedback:ISRF_c_hyp_scheme, read from the run's
    used_parameters.yml:

        Q = [sum_i m_i u_i / c_hyp,i] / [sum_i m_i / c_hyp,i]   schemes 3, 4
        Q = [sum_i m_i u_i] / [sum_i m_i]                       schemes 0, 1, 2

    The consistent-variable-c schemes (3, and 4, the shipped default) rewrite
    every pairwise operator as c_hyp_i/c times the true-speed equation, so
    their transport conserves sum m u / c_hyp and NOT sum m u; the
    shared-pair-speed schemes conserve sum m u. Both statements are made by
    radiation_propagation_iact.h at the two dispatch branches that implement
    them. Measuring the second ledger on a scheme that conserves the first
    reports the receiver-weighted redistribution as an error.

    The c_hyp,i weights are the HyperbolicPropagationSpeeds snapshot field.
    A snapshot written before that field existed, or one whose c_hyp is not
    everywhere finite and positive, degrades to the sum m u ledger with a
    printed message: exact for schemes 0 to 2, approximate for 3 and 4.

    Both forms are ratios of two sums at the SAME time, so a spatially
    uniform c_hyp cancels between numerator and denominator, whether or not
    it varies from one snapshot to the next. Every pinned run (c_hyp_pin >
    0), every fixed-fraction run (scheme 2) and every scheme-0/1 run
    therefore gets the number this check reported before the ledger became
    scheme-aware, up to round-off: the weighted branch divides each mass by
    c_hyp before summing, so the two are not the same float expression.

    A1 itself is not evaluated: with the module's own light-speed clamp,
    c_hyp <= c always, and on every fixture this file runs c_hyp/c is of
    order 1e-5 (c_hyp ~ margin*h/dt is a resolved-region speed, km/s-scale,
    against c ~ 3e5 km/s in this unit system), so the integral above is many
    orders below this check's own float32/discretisation floor over the
    run's span. The practical prediction is therefore Q(t) = Q(0), with the
    un-modelled decay folded into the bar as an explicit term bounded by
    c_hyp_bound/c (c_hyp_bound an upper estimate from the run's own
    smoothing length and step size, printed below), not asserted as exactly
    0.

    The unshielded H2 photodissociation rate the module hands to Grackle is
    ``k = sigma_H2 c rho u_LW / E_LW``, with rho = rho0 (a0/a)^3 and, to the
    same leading order as A1 (u_LW ~= u_LW,0, not u_LW,0 a0/a as an undilated
    Hubble term would give), so

        ln[x_H2(t)/x_H2(0)] = -k0 int_0^t (a0/a)^3 dt'    (-k0 t without).   (A2)

    This is the more discriminating of the two: the exponent changed from 4
    (density cubed times a linearly-decaying field) to 3 (density alone, the
    field no longer decaying at leading order), a ~10% shift in the
    predicted x_H2 over this fixture's span, well above any noise floor.

dust_absorption
    Seeded field, solar metallicity, propagation speed pinned to c_pin (so,
    unlike free_field, c_hyp/c is a single run-wide constant here, not a
    per-particle/per-step quantity). The exact solution of the module's
    relaxation update is

        ln[u(t)/u0] = -c_pin kappa0 int_0^t (a0/a)^3 dt'
                      - (c_pin/c) ln[a(t)/a0] ,                              (B1)

    the Hubble term dilated by the same c_pin/c factor as the absorption
    term (see free_field's own note above); with c_pin a few km/s against
    c ~ 3e5 km/s in this unit system, that second term is ~1e-5 of what an
    undilated -ln[a(t)/a0] would give, negligible next to the dust-absorption
    term for any metal-enriched fixture. This leg is not run as part of this
    check's own verification (see the ISRF cosmological Hubble-term rescale
    fix's own log): the correction is algebraically exact given a pinned
    c_hyp (no simulation needed to derive it), and the dust-absorption term
    dominates B1 by many orders of magnitude here, so this fixture does not
    discriminate the fix either way; the formula and code below are kept
    accurate regardless, so a future run is not compared against a
    knowingly-stale reference.

    with kappa0 = sigma_d (Z/0.01295) rho0 / (1.4 m_H) the linear absorption
    coefficient at the start (sigma_d = 9e-22 and 1.5e-21 cm^2 for PE and LW).

photoelectric
    Seeded G0, solar metallicity, low pinned speed so G0 barely changes.
    Grackle's constant-efficiency photoelectric heating
    (``photoelectric_heating = 2``, cool1d_multi_g.F) is

        Gamma = 1e-24 * 0.05 * G0 * n_H * Z/0.01295   erg cm^-3 s^-1 ,       (D1)

    for T < 2e4 K. The same fixture without the field (``photoelectric_dark``)
    carries every other heating and cooling term, and expansion, so

        u_on(t) - u_dark(t) = int_0^t Gamma / rho dt' .                      (D2)

injection
    Propagation off, no dust: the injection kernel weights sum to 1, so
    ``sum_j m_j u_j = Delta_t L`` per band, Delta_t the star's step read from
    the run log.

Bars
----
Each bar is the sum of terms stated in the output, derived from the run's
discretisation: the non-cosmological run of the same configuration measures
the error that does not depend on a (float32 updates, flux divergence on the
glass, Grackle's implicit solve); ``--reference`` passes it to the
cosmological check, whose bar is the larger of the a-priori budget and twice
that measured error, plus the cosmological terms:

- H and the rates are frozen at the step end (``cosmology_update`` runs before
  the step's tasks). For a rate r(a) ~ a^-p, the ln error per step is
  (p/2) dlna_step * r dt (matter domination, d ln H/d ln a = -3/2 adds 3/4
  for the H term), summed over the run.
- Particles are updated at their step ends, so a snapshot can lag by one
  step: one step's worth of the change.
"""

import argparse
import glob
import sys
from typing import Dict, List, Optional

import h5py
import numpy as np
from scipy.integrate import quad

SIGMA_H2_LW_CGS = 2.5111667e-18
LW_PHOTON_ENERGY_CGS = 12.2 * 1.602176634e-12
HABING_FLUX_CGS = 1.6e-3
SIGMA_D_CGS = {"PE": 9e-22, "LW": 1.5e-21}
GRACKLE_DEFAULT_DUST_TO_GAS_RATIO = 0.009387
# radiation.h's own RADIATION_HYDROGEN_MASS_CGS, which the extinction chain
# uses; M_H_CGS below is the physical constant the photoelectric rate uses.
RADIATION_HYDROGEN_MASS_CGS = 1.6726219e-24
KERNEL_GAMMA_DEFAULT = 1.936492
MU_H = 1.4
GRACKLE_SOLAR_METAL_FRACTION = 0.01295
C_LIGHT_CGS = 2.99792458e10
M_H_CGS = 1.67262171e-24
HYDROGEN_MASS_FRACTION = 0.76
PHOTOELECTRIC_RATE_CGS = 1e-24 * 0.05
FLOAT32_EPS = np.finfo(np.float32).eps
# enum isrf_c_hyp_scheme values whose pairwise operators conserve
# sum m u / c_hyp rather than sum m u (feedback_properties.h).
VARIABLE_C_SCHEMES = (3, 4)


def parse_options() -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--config",
        required=True,
        choices=[
            "free_field",
            "dust_absorption",
            "photoelectric",
            "injection",
            "injection_dusty",
        ],
    )
    parser.add_argument("-s", "--snapshots", required=True, help="Snapshot glob")
    parser.add_argument(
        "--reference",
        default=None,
        help="Snapshot glob of the non-cosmological run of the same configuration",
    )
    parser.add_argument(
        "--dark", default=None, help="Snapshot glob of photoelectric_dark"
    )
    parser.add_argument(
        "--reference-dark",
        default=None,
        help="Snapshot glob of the non-cosmological photoelectric_dark run",
    )
    parser.add_argument("--log", default=None, help="Run log, for injection")
    parser.add_argument(
        "--dt-max",
        type=float,
        default=None,
        help="TimeIntegration:dt_max of the run (ln a with cosmology); "
        "read from used_parameters.yml next to snap/ when omitted",
    )
    parser.add_argument("--c-hyp-pin", type=float, default=None, help="km/s")
    parser.add_argument(
        "--dust-tol",
        type=float,
        default=1e-4,
        help="injection_dusty only: max allowed max_j |R_j| on the band-ratio "
        "gate (default: %(default)s), two decades above this reconstruction's "
        "measured float32 floor. Pass a negative value to report the residual "
        "without gating on it.",
    )
    parser.add_argument(
        "--kernel-gamma",
        type=float,
        default=KERNEL_GAMMA_DEFAULT,
        help="injection_dusty only: the binary's kernel_gamma (default: "
        "%(default)s, Wendland C2 in 3D).",
    )
    parser.add_argument(
        "--extinction-path",
        type=float,
        default=None,
        help="injection_dusty only: GEARFeedback:ISRF_extinction_path in "
        "kernel support radii (constant_kernel_path = its own float, "
        "kernel support radii). Read "
        "from the run's used_parameters.yml when omitted.",
    )
    parser.add_argument(
        "--dust-to-gas-ratio",
        type=float,
        default=GRACKLE_DEFAULT_DUST_TO_GAS_RATIO,
        help="injection_dusty only: the run's resolved Grackle "
        "chemistry_data.local_dust_to_gas_ratio (default: %(default)s, "
        "Grackle's own compiled default).",
    )
    return parser.parse_args()


def physical(dataset: h5py.Dataset, a: float, unit_cgs: float) -> np.ndarray:
    """Return a snapshot dataset in physical CGS."""
    exponent = float(np.atleast_1d(dataset.attrs["a-scale exponent"])[0])
    return dataset[:].astype(np.float64) * a**exponent * unit_cgs


def read_snapshot(filename: str) -> Dict:
    """Read one snapshot in physical CGS."""
    with h5py.File(filename, "r") as handle:
        units = handle["/Units"].attrs
        length = float(np.atleast_1d(units["Unit length in cgs (U_L)"])[0])
        mass = float(np.atleast_1d(units["Unit mass in cgs (U_M)"])[0])
        time = float(np.atleast_1d(units["Unit time in cgs (U_t)"])[0])
        header = handle["/Header"].attrs
        a = float(np.atleast_1d(header["Scale-factor"])[0])
        cosmo = handle["/Cosmology"].attrs
        is_cosmo = int(np.atleast_1d(cosmo.get("Cosmological run", [0]))[0]) == 1
        velocity = length / time
        energy = velocity**2
        gas = handle["/PartType0"]
        order = np.argsort(gas["ParticleIDs"][:])
        out = {
            "a": a,
            "cosmological": is_cosmo,
            "time": float(np.atleast_1d(header["Time"])[0]) * time,
            "cosmology": {
                key: float(np.atleast_1d(cosmo[key])[0])
                for key in [
                    "H0 [internal units]",
                    "Omega_m",
                    "Omega_r",
                    "Omega_k",
                    "Omega_lambda",
                ]
            },
            "time_unit": time,
            "mass_unit": mass,
            "energy_unit": energy,
            "density": physical(gas["Densities"], a, mass / length**3)[order],
            "mass": physical(gas["Masses"], a, mass)[order],
            "h": physical(gas["SmoothingLengths"], a, length)[order],
            "u": physical(gas["InternalEnergies"], a, energy)[order],
            "u_PE": physical(gas["PESpecificEnergies"], a, energy)[order],
            "u_LW": physical(gas["LWSpecificEnergies"], a, energy)[order],
            "c_hyp": (
                physical(gas["HyperbolicPropagationSpeeds"], a, velocity)[order]
                if "HyperbolicPropagationSpeeds" in gas
                else None
            ),
            "H2I": gas["H2I"][:].astype(np.float64)[order],
            "hydrogen": sum(
                gas[name][:].astype(np.float64)[order]
                for name in ["HI", "HII", "H2I", "H2II"]
            ),
        }
        metals = gas["MetalMassFractions"][:].astype(np.float64)
        out["Z"] = (metals[:, -1] if metals.ndim == 2 else metals)[order]
        # The extinction chain reads the SMOOTHED metal mass fraction
        # (chemistry_get_total_metal_mass_fraction_for_cooling), as a float32.
        # It coincides with the unsmoothed array only at Z = 0.
        if "SmoothedMetalMassFractions" in gas:
            smoothed = gas["SmoothedMetalMassFractions"][:]
            out["Z_smoothed"] = (smoothed[:, -1] if smoothed.ndim == 2 else smoothed)[
                order
            ].astype(np.float32)
        if "/PartType4" in handle and handle["/PartType4/Masses"].shape[0] > 0:
            out["L_PE"] = float(handle["/PartType4/PELuminosities"][0])
            out["L_LW"] = float(handle["/PartType4/LWLuminosities"][0])
            out["time_internal"] = float(np.atleast_1d(header["Time"])[0])
    return out


def load_run(pattern: str) -> List[Dict]:
    """Read every snapshot of a run, sorted by time."""
    files = sorted(glob.glob(pattern))
    if len(files) < 2:
        raise RuntimeError(f"Need at least two snapshots for {pattern!r}")
    run = []
    for name in files:
        snap = read_snapshot(name)
        # SWIFT also dumps at time_end, which can repeat the last output time.
        if run and snap["time"] == run[-1]["time"]:
            continue
        run.append(snap)
    return run


def hubble_rate_cgs(a: float, snap: Dict) -> float:
    """Return H(a) in s^-1 from the snapshot's cosmology (0 without cosmology)."""
    if not snap["cosmological"]:
        return 0.0
    c = snap["cosmology"]
    e2 = (
        c["Omega_r"] * a**-4
        + c["Omega_m"] * a**-3
        + c["Omega_k"] * a**-2
        + c["Omega_lambda"]
    )
    return c["H0 [internal units]"] * np.sqrt(e2) / snap["time_unit"]


def power_integral(run: List[Dict], power: float) -> np.ndarray:
    """Return int_0^t (a0/a)^power dt' at every snapshot, in seconds."""
    first = run[0]
    if not first["cosmological"]:
        return np.array([s["time"] - first["time"] for s in run])
    a0 = first["a"]
    return np.array(
        [
            (
                quad(
                    lambda x: np.exp(-power * x)
                    / hubble_rate_cgs(a0 * np.exp(x), first),
                    0.0,
                    np.log(s["a"] / a0),
                    epsabs=0.0,
                    epsrel=1e-12,
                )[0]
                if s["a"] > a0
                else 0.0
            )
            for s in run
        ]
    )


def friedmann_time_residual(run: List[Dict]) -> float:
    """Return max |t_Friedmann/t_SWIFT - 1| over the snapshots (0 without cosmology)."""
    if not run[0]["cosmological"]:
        return 0.0
    t_model = power_integral(run, 0.0)[1:]
    t_swift = np.array([s["time"] - run[0]["time"] for s in run])[1:]
    return float(np.max(np.abs(t_model / t_swift - 1.0)))


def read_dt_max(pattern: str, given: Optional[float]) -> float:
    """Return dt_max, from the argument or the run's used_parameters.yml."""
    if given is not None:
        return given
    import os
    import yaml

    directory = os.path.dirname(os.path.dirname(sorted(glob.glob(pattern))[0]))
    with open(os.path.join(directory, "used_parameters.yml")) as handle:
        return float(yaml.safe_load(handle)["TimeIntegration"]["dt_max"])


def read_extinction_path(pattern: str, given: Optional[float]) -> float:
    """Return the extinction path in kernel support radii, from the argument
    or the run's used_parameters.yml."""
    if given is not None:
        return given
    import os
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
    if name in ("pair_separation", "temperature_capped_jeans"):
        raise RuntimeError(
            f"GEARFeedback:ISRF_extinction_path {name!r} does not build the "
            "column from a multiple of the kernel support radius, so the "
            "closed form this check mirrors does not apply to it. Rerun the "
            "fixture with a constant_kernel_path mechanism, or extend this "
            "check to that mechanism's own length."
        )
    raise RuntimeError(f"Unknown GEARFeedback:ISRF_extinction_path {name!r}")


def optical_depths(
    snap: Dict, path_in_kernel_radii: float, kernel_gamma: float, dust_to_gas: float
) -> Dict[str, np.ndarray]:
    """Return each gas particle's PE and LW dust optical depth.

    Mirrors radiation_get_part_ISRF_extinction_factors and the chain below
    it in radiation_isrf.c, entirely in physical CGS: converting the
    comoving column to a physical one is exactly what the code's a^-2 does,
    so no scale factor appears here beyond the per-dataset ones read_snapshot
    already applied.

    Parameters
    ----------
    snap
        One snapshot as returned by read_snapshot.
    path_in_kernel_radii
        GEARFeedback:ISRF_extinction_path, in kernel support radii.
    kernel_gamma
        The binary's kernel_gamma.
    dust_to_gas
        The run's resolved Grackle chemistry_data.local_dust_to_gas_ratio.

    Returns
    -------
    dict
        The dimensionless optical depth of each particle, keyed by band.
    """
    column = path_in_kernel_radii * kernel_gamma * snap["h"] * snap["density"]
    d_relative = (
        np.maximum(snap["Z_smoothed"].astype(np.float64), 0.0)
        / GRACKLE_SOLAR_METAL_FRACTION
        * (dust_to_gas / GRACKLE_DEFAULT_DUST_TO_GAS_RATIO)
    )
    prefactor = d_relative / (MU_H * RADIATION_HYDROGEN_MASS_CGS) * column
    return {band: SIGMA_D_CGS[band] * prefactor for band in ("PE", "LW")}


def read_c_hyp_margin(pattern: str) -> float:
    """Return GEARFeedback:ISRF_c_hyp_margin from the run's used_parameters.yml."""
    import os
    import yaml

    directory = os.path.dirname(os.path.dirname(sorted(glob.glob(pattern))[0]))
    with open(os.path.join(directory, "used_parameters.yml")) as handle:
        return float(yaml.safe_load(handle)["GEARFeedback"]["ISRF_c_hyp_margin"])


def read_c_hyp_scheme(pattern: str) -> Optional[int]:
    """Return GEARFeedback:ISRF_c_hyp_scheme, or None when it is not recorded."""
    import os
    import yaml

    directory = os.path.dirname(os.path.dirname(sorted(glob.glob(pattern))[0]))
    path = os.path.join(directory, "used_parameters.yml")
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        parameters = yaml.safe_load(handle)
    try:
        return int(parameters["GEARFeedback"]["ISRF_c_hyp_scheme"])
    except (KeyError, TypeError, ValueError):
        return None


def use_c_hyp_ledger(run: List[Dict], pattern: str, label: str) -> bool:
    """Report whether the c_hyp-weighted ledger applies to this run.

    It applies only to the consistent-variable-c schemes, and only when every
    snapshot carries a finite, strictly positive HyperbolicPropagationSpeeds.
    Every rejection prints why, so a degraded run is never silently gated on
    the wrong invariant.
    """
    scheme = read_c_hyp_scheme(pattern)
    if scheme is None:
        print(f"  {label}: ISRF_c_hyp_scheme not recorded, using the sum m u " "ledger")
        return False
    if scheme not in VARIABLE_C_SCHEMES:
        return False
    if any(s["c_hyp"] is None for s in run):
        print(
            f"  {label}: scheme {scheme} conserves sum m u / c_hyp, but "
            "HyperbolicPropagationSpeeds is absent from these snapshots; "
            "falling back to the sum m u ledger, which over-reports the "
            "receiver-weighted redistribution as an error"
        )
        return False
    for snap in run:
        if not np.all(np.isfinite(snap["c_hyp"])) or np.any(snap["c_hyp"] <= 0.0):
            print(
                f"  {label}: HyperbolicPropagationSpeeds is not everywhere "
                "finite and positive (propagation off, or a pre-first-step "
                "snapshot); falling back to the sum m u ledger"
            )
            return False
    return True


def summarize(label: str, error: np.ndarray) -> float:
    """Print and return the worst per-snapshot median |error|."""
    medians = np.median(np.abs(error), axis=1)
    worst = float(np.max(medians))
    print(
        f"  {label:<34s} worst median |err| {worst:.3e}, "
        f"worst particle {np.max(np.abs(error)):.3e}"
    )
    return worst


def gate(label: str, worst: float, bar: float) -> bool:
    """Print a pass/fail line."""
    ok = bool(np.isfinite(worst)) and worst <= bar
    print(
        f"  {'PASS' if ok else 'FAIL'}: {label}: {worst:.3e} <= bar {bar:.3e}"
        if ok
        else f"  FAIL: {label}: {worst:.3e} > bar {bar:.3e}"
    )
    return ok


def step_count(run: List[Dict], dt_max: float) -> float:
    """Return the number of dt_max steps spanned by the run."""
    first, last = run[0], run[-1]
    if first["cosmological"]:
        return np.log(last["a"] / first["a"]) / dt_max
    return (last["time"] - first["time"]) / (dt_max * first["time_unit"])


def ledger_mean(snap: Dict, key: str, use_c_hyp: bool) -> float:
    """Return one band's box mean under the ledger the run's scheme conserves.

    With ``use_c_hyp`` the particle weight is ``m_i/c_hyp,i`` instead of
    ``m_i``. Numerator and denominator are summed at the same time, so a
    spatially uniform c_hyp cancels between them and the two branches then
    agree to round-off, whether or not c_hyp varies from snapshot to
    snapshot.
    """
    mass = snap["mass"]
    if not use_c_hyp:
        return float(np.sum(mass * snap[key]) / np.sum(mass))
    weight = mass / snap["c_hyp"]
    return float(np.sum(weight * snap[key]) / np.sum(weight))


def free_field_errors(run: List[Dict], use_c_hyp: bool = False) -> Dict:
    """Return the errors of Eqs. (A1) and (A2) at every snapshot.

    The transport moves energy between particles and conserves the ledger of
    the run's own c_hyp scheme (this module's docstring), so on a glass each
    particle's field departs from the uniform solution by the glass noise
    while the box mean follows (A1) exactly. The gates use the box means; the
    per-particle spread is reported.

    A1's reference is the ledger's own initial value (see this module's own
    docstring: to the precision this check can resolve, the c_hyp/c-dilated
    Hubble decay is un-modelled, not asserted as exactly 0), so `out[band]`
    is the box mean's own fractional departure from it, not from an
    a-dependent target.

    Parameters
    ----------
    run
        The run's snapshots, in time order.
    use_c_hyp
        Weight each particle by ``m_i/c_hyp,i`` rather than ``m_i``, for the
        consistent-variable-c schemes. Decided by `use_c_hyp_ledger`.
    """
    first = run[0]
    mass = first["mass"]
    out = {}
    for band in ["PE", "LW"]:
        u0 = ledger_mean(first, f"u_{band}", use_c_hyp)
        out[band] = np.array(
            [ledger_mean(s, f"u_{band}", use_c_hyp) / u0 - 1.0 for s in run]
        )
        out[f"{band}_spread"] = np.array(
            [np.median(np.abs(s[f"u_{band}"] / first[f"u_{band}"] - 1.0)) for s in run]
        )
    # Box-mean density and field: the closed form is for the uniform state.
    rho0 = np.sum(mass) / np.sum(mass / first["density"])
    u_lw0 = np.sum(mass * first["u_LW"]) / np.sum(mass)
    k0 = SIGMA_H2_LW_CGS * C_LIGHT_CGS * rho0 * u_lw0 / LW_PHOTON_ENERGY_CGS
    # Power 3, not 4: rho ~ (a0/a)^3 alone now (A2, this module's docstring).
    # u_LW no longer contributes an (a0/a)^1 factor once the Hubble term is
    # correctly dilated by c_hyp/c.
    integral = power_integral(run, 3.0)
    measured = np.array([np.mean(np.log(s["H2I"] / first["H2I"])) for s in run])
    predicted = -k0 * integral
    out["H2"] = (measured[1:] - predicted[1:]) / np.abs(predicted[1:])
    out["exponent"] = float(-predicted[-1])
    out["rate"] = k0
    out["integral"] = integral
    return out


def check_free_field(opt: argparse.Namespace) -> bool:
    """Check Eqs. (A1) and (A2)."""
    run = load_run(opt.snapshots)
    dt_max = read_dt_max(opt.snapshots, opt.dt_max)
    cosmological = run[0]["cosmological"]
    n_steps = step_count(run, dt_max)
    use_c_hyp = use_c_hyp_ledger(run, opt.snapshots, "run")
    ledger = "sum m u / c_hyp" if use_c_hyp else "sum m u"
    errors = free_field_errors(run, use_c_hyp)
    span = float(np.log(run[-1]["a"] / run[0]["a"]))
    elapsed = np.array([s["time"] - run[0]["time"] for s in run])[1:]
    print(
        f"free_field: cosmological={cosmological}, a {run[0]['a']:.6g} -> "
        f"{run[-1]['a']:.6g}, {len(run)} snapshots, {n_steps:.0f} dt_max steps, "
        f"H2 exponent {errors['exponent']:.3f}, ledger {ledger}"
    )
    print(
        f"  Friedmann time vs SWIFT time: max rel. diff {friedmann_time_residual(run):.2e}"
    )
    for band in ["PE", "LW"]:
        print(
            f"  per-particle spread of u_{band} (glass noise, not gated): worst median "
            f"{np.max(errors[f'{band}_spread']):.3e}"
        )

    reference = None
    if opt.reference:
        reference_run = load_run(opt.reference)
        use_c_hyp_reference = use_c_hyp_ledger(
            reference_run, opt.reference, "reference"
        )
        if use_c_hyp_reference != use_c_hyp:
            print(
                "  WARNING: the reference run uses the other ledger, so its "
                "measured error is not comparable with this run's"
            )
        reference = free_field_errors(reference_run, use_c_hyp_reference)

    ok = True
    # Longest step in proper time, from dt_max (ln a with cosmology).
    dt_step = (
        dt_max / hubble_rate_cgs(run[0]["a"], run[0])
        if cosmological
        else dt_max * run[0]["time_unit"]
    )
    # Upper-bound estimate of c_hyp/c, from the run's own resolution: c_hyp is
    # not a snapshot field (this fixture's c_hyp_pin is off, so c_hyp itself
    # varies per particle and per step, see this module's own docstring), but
    # it is bounded by ISRF_c_hyp_margin*h/dt_step (h the smallest physical
    # dt_step gives the largest bound) and by c. dt_step above already uses
    # H(a0), the largest H over this matter-domination run, hence the
    # smallest physical dt: a conservative (not tight) upper bound.
    c_hyp_ratio = 0.0
    if cosmological:
        margin = read_c_hyp_margin(opt.snapshots)
        h_phys_cgs = np.median(run[0]["h"])
        c_hyp_bound_cgs = min(C_LIGHT_CGS, margin * h_phys_cgs / dt_step)
        c_hyp_ratio = c_hyp_bound_cgs / C_LIGHT_CGS
        print(
            f"  c_hyp/c upper bound: {c_hyp_ratio:.3e} (margin {margin:g}, "
            f"median h {h_phys_cgs:.3e} cm, dt_step {dt_step:.3e} s)"
        )
    for band in ["PE", "LW"]:
        # Float32 round-off of each update, averaged over the particles.
        budget = FLOAT32_EPS * n_steps / np.sqrt(run[0]["mass"].size)
        measured_nc = (
            0.0 if reference is None else float(np.max(np.abs(reference[band])))
        )
        # The un-modelled dilated decay itself (span), plus its own step-end
        # H ((3/4) dlna_step per unit ln a) and one-step lag discretisation,
        # all dilated by the same c_hyp/c factor as the term itself.
        cosmo = (
            c_hyp_ratio * (span + 0.75 * dt_max * span + dt_max)
            if cosmological
            else 0.0
        )
        bar = max(budget, 2.0 * measured_nc) + cosmo
        worst = float(np.max(np.abs(errors[band])))
        print(
            f"  bar u_{band}: max(float32 {budget:.1e}, 2 x non-cosmological "
            f"{2.0 * measured_nc:.1e}) + c_hyp/c-dilated decay, step-end H "
            f"and lag {cosmo:.1e}"
        )
        ok &= gate(f"box-mean u_{band} / u0 - 1 (A1), {ledger}", worst, bar)

    # H2, per snapshot: implicit solve (k dt/2), one-step snapshot lag
    # (dt/t without cosmology; with it the rate varies as a^-3, same order),
    # float32 round-off; with cosmology the step-end rate adds 1.5 dlna_step
    # ((p/2) dlna_step at p = 3, this module's own Bars section, not 2.0 at
    # p = 4: A2's rate no longer carries u_LW's own a0/a factor).
    budget = 0.5 * errors["rate"] * dt_step + dt_step / elapsed + FLOAT32_EPS * n_steps
    cosmo = 1.5 * dt_max if cosmological else 0.0
    measured_nc = np.zeros_like(budget)
    if reference is not None:
        # Snapshot by snapshot when both runs have the same output count: the
        # lag term is largest at the first snapshots in both runs.
        if reference["H2"].size == budget.size:
            measured_nc = np.abs(reference["H2"])
        else:
            measured_nc = np.full_like(budget, np.max(np.abs(reference["H2"])))
    bar = np.maximum(budget, 2.0 * measured_nc) + cosmo
    ratio = np.abs(errors["H2"]) / bar
    k = int(np.argmax(ratio))
    print(
        f"  bar ln x_H2 (per snapshot): implicit solve {0.5 * errors['rate'] * dt_step:.1e} "
        f"+ lag dt/t {dt_step / elapsed[-1]:.1e} (end) to {dt_step / elapsed[0]:.1e} (first), "
        f"2 x non-cosmological up to {2.0 * np.max(measured_nc):.1e}, step-end rate {cosmo:.1e}"
    )
    print(
        f"  worst snapshot {k + 1}: error {errors['H2'][k]:.3e}, bar {bar[k]:.3e}; "
        f"final error {errors['H2'][-1]:.3e}, bar {bar[-1]:.3e}"
    )
    ok &= gate("box-mean ln x_H2 exponent (A2), worst error/bar", float(ratio[k]), 1.0)
    return ok


def dust_absorption_errors(run: List[Dict], c_pin_cgs: float) -> Dict:
    """Return the per-snapshot error of Eq. (B1) on the box, in ln u.

    The transport conserves sum m u and mixes the field between neighbours,
    so a particle does not decay at its own kappa (its SPH density scatters by
    a few 1e-3 on the glass) but at the neighbourhood mean. The box sum decays
    exactly at the mass-weighted mean kappa for a uniform field, which is what
    is compared. kappa is proportional to the density; the drift of the
    mass-weighted mean comoving density over the run is returned for the bar.
    """
    first = run[0]
    a0 = first["a"]
    mass = first["mass"]
    integral = power_integral(run, 3.0)
    ln_a = np.array([np.log(s["a"] / a0) for s in run])
    comoving_mean = np.array(
        [
            np.sum(s["mass"] * s["density"]) / np.sum(s["mass"]) * (s["a"] / a0) ** 3
            for s in run
        ]
    )
    out = {
        "density_drift": float(np.max(np.abs(comoving_mean / comoving_mean[0] - 1.0)))
    }
    rho0 = comoving_mean[0]
    z0 = np.sum(mass * first["Z"]) / np.sum(mass)
    for band in ["PE", "LW"]:
        kappa0 = (
            SIGMA_D_CGS[band]
            * (z0 / GRACKLE_SOLAR_METAL_FRACTION)
            * rho0
            / (MU_H * M_H_CGS)
        )
        # The Hubble term dilated by c_pin/c, same factor as the absorption
        # term (this module's own docstring, B1): negligible here (c_pin is
        # km/s-scale) but kept exact rather than dropped.
        predicted = -c_pin_cgs * kappa0 * integral - (c_pin_cgs / C_LIGHT_CGS) * ln_a
        total0 = np.sum(mass * first[f"u_{band}"])
        measured = np.array(
            [np.log(np.sum(s["mass"] * s[f"u_{band}"]) / total0) for s in run]
        )
        out[band] = (measured - predicted)[:, None]
        out[f"{band}_depth"] = float(-predicted[-1])
    return out


def check_dust_absorption(opt: argparse.Namespace) -> bool:
    """Check Eq. (B1)."""
    if opt.c_hyp_pin is None:
        raise RuntimeError("--c-hyp-pin (km/s) is required for dust_absorption")
    run = load_run(opt.snapshots)
    dt_max = read_dt_max(opt.snapshots, opt.dt_max)
    cosmological = run[0]["cosmological"]
    n_steps = step_count(run, dt_max)
    errors = dust_absorption_errors(run, opt.c_hyp_pin * 1e5)
    span = np.log(run[-1]["a"] / run[0]["a"])
    print(
        f"dust_absorption: cosmological={cosmological}, {len(run)} snapshots, "
        f"{n_steps:.0f} dt_max steps, final ln depth PE {errors['PE_depth']:.3f}, "
        f"LW {errors['LW_depth']:.3f}; median Z {np.median(run[0]['Z']):.4g}"
    )
    worst = {
        band: summarize(f"box ln sum m u_{band} (B1)", errors[band])
        for band in ["PE", "LW"]
    }
    nc = {"PE": 0.0, "LW": 0.0}
    if opt.reference:
        ref = dust_absorption_errors(load_run(opt.reference), opt.c_hyp_pin * 1e5)
        nc = {
            band: float(np.max(np.median(np.abs(ref[band]), axis=1)))
            for band in ["PE", "LW"]
        }
        print(
            f"  non-cosmological reference errors: PE {nc['PE']:.3e}, LW {nc['LW']:.3e}"
        )
    ok = True
    for band in ["PE", "LW"]:
        depth = errors[f"{band}_depth"]
        per_step = depth / max(n_steps, 1.0)
        budget = (
            FLOAT32_EPS * n_steps / np.sqrt(run[0]["mass"].size)
            + per_step
            + depth * errors["density_drift"]
        )
        # The kappa step-end term is unaffected by the Hubble-term dilation
        # (kappa was already correctly dilated); the H-alone step-end term is
        # dilated by c_pin/c, same as the B1 formula's own second term above.
        c_pin_ratio = (opt.c_hyp_pin * 1e5) / C_LIGHT_CGS if cosmological else 0.0
        cosmo = (
            1.5 * dt_max * depth + c_pin_ratio * 0.75 * dt_max * span
            if cosmological
            else 0.0
        )
        bar = max(budget, 2.0 * nc[band]) + cosmo
        print(
            f"  bar {band}: max(float32 + one-step lag + density drift "
            f"{errors['density_drift']:.1e} x depth = {budget:.1e}, 2 x reference "
            f"{2 * nc[band]:.1e}) + step-end kappa and H {cosmo:.1e}"
        )
        ok &= gate(f"box ln sum m u_{band} (B1)", worst[band], bar)
    return ok


def photoelectric_errors(on: List[Dict], dark: List[Dict]) -> Dict:
    """Return the per-snapshot relative error of Eq. (D2) and its inputs."""
    if len(on) != len(dark):
        raise RuntimeError("photoelectric and photoelectric_dark differ in snapshots")
    times = np.array([s["time"] - on[0]["time"] for s in on])
    heating = []
    for s in on:
        g0 = C_LIGHT_CGS * s["density"] * (s["u_PE"] + s["u_LW"]) / HABING_FLUX_CGS
        # Grackle's rhoH: the hydrogen species, not the primordial fraction.
        n_h = s["hydrogen"] * s["density"] / M_H_CGS
        gamma = (
            PHOTOELECTRIC_RATE_CGS * g0 * n_h * s["Z"] / GRACKLE_SOLAR_METAL_FRACTION
        )
        heating.append(np.median(gamma / s["density"]))
    heating = np.array(heating)
    predicted = np.concatenate(
        [[0.0], np.cumsum(0.5 * (heating[1:] + heating[:-1]) * np.diff(times))]
    )
    measured = np.array([np.median(s["u"] - d["u"]) for s, d in zip(on, dark)])
    measured -= measured[0]
    dark_u = np.array([np.median(d["u"]) for d in dark])
    on_u = np.array([np.median(s["u"]) for s in on])
    return {
        "times": times,
        "heating": heating,
        "predicted": predicted,
        "relative": (measured[1:] - predicted[1:]) / predicted[1:],
        "dark_u": dark_u,
        "on_u": on_u,
    }


def check_photoelectric(opt: argparse.Namespace) -> bool:
    """Check Eq. (D2)."""
    if opt.dark is None:
        raise RuntimeError("--dark is required for photoelectric")
    on = load_run(opt.snapshots)
    dark = load_run(opt.dark)
    dt_max = read_dt_max(opt.snapshots, opt.dt_max)
    cosmological = on[0]["cosmological"]
    err = photoelectric_errors(on, dark)
    dt_step = (
        dt_max / hubble_rate_cgs(on[0]["a"], on[0])
        if cosmological
        else dt_max * on[0]["time_unit"]
    )
    print(
        f"photoelectric: cosmological={cosmological}, heating "
        f"{err['heating'][0]:.4e} -> {err['heating'][-1]:.4e} erg/g/s, "
        f"u_on - u_dark at end {err['predicted'][-1] * (1 + err['relative'][-1]):.4e} erg/g, "
        f"u_dark {err['dark_u'][0]:.4e} -> {err['dark_u'][-1]:.4e} erg/g"
    )
    # Terms of the per-snapshot bar:
    # - a snapshot can lag the heating by one step: dt/t;
    # - the heated gas cools faster than the dark gas. Fine-structure cooling
    #   scales as exp(-T_line/T) with T_line = 92 K (C+), so the dark run's
    #   own net loss, scaled by exp(92/T_dark - 92/T_on) - 1, bounds the
    #   difference; any T-independent loss cancels in the dark twin;
    # - float32 storage of u relative to the difference.
    t = err["times"][1:]
    temperature_ratio = err["on_u"][1:] / err["dark_u"][1:]
    # Neutral atomic gas, mu = 4/(1 + 3 X).
    t_dark = (
        (2.0 / 3.0)
        * (4.0 / (1.0 + 3.0 * HYDROGEN_MASS_FRACTION))
        * M_H_CGS
        * err["dark_u"][1:]
        / 1.380649e-16
    )
    boost = np.expm1(92.0 / t_dark * (1.0 - 1.0 / temperature_ratio))
    cooling = (
        np.abs(err["dark_u"][1:] - err["dark_u"][0]) * boost / err["predicted"][1:]
    )
    lag = dt_step / t
    storage = 4.0 * FLOAT32_EPS * err["on_u"][1:] / err["predicted"][1:]
    bar = lag + cooling + storage
    if opt.reference:
        if opt.reference_dark is None:
            raise RuntimeError("--reference needs --reference-dark for photoelectric")
        ref = photoelectric_errors(
            load_run(opt.reference), load_run(opt.reference_dark)
        )
        # Reference error at the same elapsed times.
        ref_error = np.interp(t, ref["times"][1:], np.abs(ref["relative"]))
        bar = np.maximum(bar, 2.0 * ref_error)
        print(
            f"  non-cosmological reference errors: first {ref['relative'][0]:.3e}, "
            f"final {ref['relative'][-1]:.3e}"
        )
    ratio = np.abs(err["relative"]) / bar
    k = int(np.argmax(ratio))
    print(f"  errors: first {err['relative'][0]:.3e}, final {err['relative'][-1]:.3e}")
    print(
        f"  bar terms at the end: lag {lag[-1]:.1e}, cooling change {cooling[-1]:.1e}, "
        f"float32 {storage[-1]:.1e}; worst snapshot {k + 1}: error "
        f"{err['relative'][k]:.3e}, bar {bar[k]:.3e}"
    )
    return gate("photoelectric heating (D2), worst error/bar", float(ratio[k]), 1.0)


def check_injection(opt: argparse.Namespace) -> bool:
    """Check one injection pass on the last snapshot.

    At zero metallicity the extinction factor is 1 on every particle, so
    sum_j m_j u_j = Delta_t L exactly. With dust that identity is false, and
    config=injection_dusty gates the two weaker exact statements instead:
    the per-particle band ratio, and the bracket on the weighted mean of
    exp(-tau). See the ISRFInjectionConservation check, which carries the
    same metric and the full derivation, for what they do and do not test.
    """
    if opt.log is None:
        raise RuntimeError("--log is required for injection")
    dusty = opt.config == "injection_dusty"
    run = load_run(opt.snapshots)
    last = run[-1]
    if not dusty and np.any(last["Z"] != 0.0):
        raise RuntimeError(
            "injection needs zero metallicity; use config=injection_dusty to "
            "check the extinction identities at nonzero metallicity instead"
        )
    if dusty and "Z_smoothed" not in last:
        raise RuntimeError(
            "injection_dusty needs the SmoothedMetalMassFractions snapshot "
            "field, which is the array the extinction chain reads"
        )
    if dusty and not np.any(last["Z_smoothed"] > 0.0):
        raise RuntimeError("injection_dusty needs a nonzero metallicity")
    # The log prints Time with 7 significant digits, so take the step row
    # closest to the snapshot time and require it to lie within half a step.
    delta_t = None
    best = np.inf
    with open(opt.log) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                int(fields[0])
                t = float(fields[1])
                dt = float(fields[4])
            except ValueError:
                continue
            distance = abs(t - last["time_internal"])
            if distance < best and distance <= 0.5 * dt:
                best = distance
                delta_t = dt
    if delta_t is None:
        raise RuntimeError("No step in the log matches the last snapshot's time")
    ok = True
    print(
        f"injection: cosmological={last['cosmological']}, a {last['a']:.6g}, "
        f"Delta_t {delta_t:.6e} internal"
    )
    if not dusty:
        for band in ["PE", "LW"]:
            lhs = np.sum(last["mass"] * last[f"u_{band}"]) / (
                last["mass_unit"] * last["energy_unit"]
            )
            rhs = delta_t * last[f"L_{band}"]
            ok &= gate(f"sum m u_{band} / (Delta_t L) - 1", abs(lhs / rhs - 1.0), 1e-5)
        return ok

    path = read_extinction_path(opt.snapshots, opt.extinction_path)
    tau = optical_depths(last, path, opt.kernel_gamma, opt.dust_to_gas_ratio)
    print(
        f"  column: path {path:g} kernel radii, kernel_gamma "
        f"{opt.kernel_gamma:g}, local_dust_to_gas_ratio "
        f"{opt.dust_to_gas_ratio:g}; smoothed Z max "
        f"{np.max(last['Z_smoothed']):.6e}"
    )

    # A NaN compares false against every bound, so it would drop out of the
    # `u > 0` selection below unnoticed rather than fail the gate.
    for name, array in (
        ("u_PE", last["u_PE"]),
        ("u_LW", last["u_LW"]),
        ("tau_PE", tau["PE"]),
        ("tau_LW", tau["LW"]),
        ("masses", last["mass"]),
    ):
        bad = int(np.sum(~np.isfinite(array)))
        if bad:
            print(f"  FAIL: {bad} non-finite {name} value(s) in the snapshot")
            return False

    for name in ("L_PE", "L_LW"):
        if not np.isfinite(last[name]) or last[name] <= 0.0:
            print(f"  FAIL: {name} must be finite and positive, got {last[name]}")
            return False

    lit_pe, lit_lw = last["u_PE"] > 0.0, last["u_LW"] > 0.0
    if int(np.sum(lit_pe)) != int(np.sum(lit_lw)):
        print(
            f"  FAIL: the bands illuminate different particle counts, "
            f"{int(np.sum(lit_pe))} PE against {int(np.sum(lit_lw))} LW"
        )
        return False
    if not np.any(lit_pe):
        print("  FAIL: no illuminated gas particle")
        return False

    sigma_ratio = SIGMA_D_CGS["PE"] / SIGMA_D_CGS["LW"]
    residual = (
        np.log(last["u_PE"][lit_pe] / last["u_LW"][lit_pe])
        - np.log(last["L_PE"] / last["L_LW"])
        - (1.0 - sigma_ratio) * tau["LW"][lit_pe]
    )
    signal = abs((1.0 - sigma_ratio) * float(np.max(tau["LW"][lit_pe])))
    print(
        f"  tau_LW {np.min(tau['LW'][lit_pe]):.4f} to "
        f"{np.max(tau['LW'][lit_pe]):.4f}, band-ratio signal {signal:.4f}, "
        f"float32 budget {4.0 * FLOAT32_EPS / 2.0 + 3.0 * FLOAT32_EPS * signal:.2e}"
    )
    worst = float(np.max(np.abs(residual)))
    if opt.dust_tol < 0.0:
        print(f"  REPORT (no bar given): band-ratio residual max_j |R_j| = {worst:.3e}")
        ok &= bool(np.isfinite(worst))
    else:
        ok &= gate("band-ratio residual (G2), max_j |R_j|", worst, opt.dust_tol)

    for band in ("PE", "LW"):
        measured = float(
            np.sum(last["mass"] * last[f"u_{band}"])
            / (last["mass_unit"] * last["energy_unit"])
            / (delta_t * last[f"L_{band}"])
        )
        low = float(np.min(np.exp(-tau[band][lit_pe])))
        high = float(np.max(np.exp(-tau[band][lit_pe])))
        inside = bool(np.isfinite(measured)) and low * (
            1.0 - 1e-5
        ) <= measured <= high * (1.0 + 1e-5)
        print(
            f"  {'PASS' if inside else 'FAIL'}: weighted-mean extinction "
            f"bracket (G1) {band}: {measured:.8f} in [{low:.8f}, {high:.8f}]"
        )
        ok &= inside
    return ok


def main() -> int:
    """Run the requested check."""
    opt = parse_options()
    checks = {
        "free_field": check_free_field,
        "dust_absorption": check_dust_absorption,
        "photoelectric": check_photoelectric,
        "injection": check_injection,
        "injection_dusty": check_injection,
    }
    ok = checks[opt.config](opt)
    print("RESULT: PASS" if ok else "RESULT: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
