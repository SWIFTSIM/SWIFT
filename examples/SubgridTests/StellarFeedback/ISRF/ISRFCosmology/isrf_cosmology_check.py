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
    radiation_end_force_propagation, fixed for this). In the unpinned legs
    of this fixture ``c_hyp`` is NOT pinned (c_hyp_pin is 0, the run.sh
    default for ``free_field``), so it is a per-particle, per-step quantity
    (ISRF_c_hyp_margin*h/dt, clamped at c), not a single constant, and the
    exact solution is the integral

        ln[Q(t)/Q(0)] = -(1/c) int_0^t c_hyp(t') H(t') dt'                  (A1)

    with Q the ledger the run's own propagation scheme conserves. Both
    schemes use the same pairwise operators, each one c_hyp_i/c times the
    true-speed equation (radiation_propagation_iact.h), so the transport
    conserves sum m u / c_hyp. Which weights the check forms depends on
    GEARFeedback:ISRF_c_hyp_scheme, read from the run's used_parameters.yml:

        Q = [sum_i m_i u_i / c_hyp,i] / [sum_i m_i / c_hyp,i]   scheme 4
        Q = [sum_i m_i u_i] / [sum_i m_i]                       scheme 2

    Scheme 4 (the default) has a c_hyp that varies from particle to particle,
    so sum m u is NOT conserved and measuring it reports the receiver-weighted
    redistribution as an error. Scheme 2 has one fixed speed for the whole
    box, where the two weightings are proportional and the ratio is the same
    number; the unweighted form is used there because it needs no
    HyperbolicPropagationSpeeds, which reads 0 on a snapshot written before the
    first drift. Any other recorded scheme value is a run from a removed
    scheme: `read_c_hyp_scheme` raises on it.

    The c_hyp,i weights are the HyperbolicPropagationSpeeds snapshot field.
    A snapshot written before that field existed, or one whose c_hyp is not
    everywhere finite and positive, degrades to the sum m u ledger with a
    printed message: exact for scheme 2, approximate for scheme 4.

    Both forms are ratios of two sums at the SAME time, so a spatially
    uniform c_hyp cancels between numerator and denominator, whether or not
    it varies from one snapshot to the next. Every pinned run (c_hyp_pin >
    0) and every fixed-fraction run (scheme 2) therefore gets the number this
    check reported before the ledger became scheme-aware, up to round-off:
    the weighted branch divides each mass by c_hyp before summing, so the two
    are not the same float expression.

    A1 is gated as a TWO-SIDED RESIDUAL, like (A3) below: the measured box-
    mean drift ln[Q(t)/Q(t_ref)] MINUS the predicted drift, against a bar
    that carries this check's own errors and nothing else. Folding the
    predicted decay into the allowance instead, and comparing an absolute
    deviation against it, passes a build with no cosmological decay at all,
    and one whose decay has the wrong sign, since such a comparison never
    looks at either.

    The predicted drift is drift_coeff * int (c_hyp/c) d ln a
    (`trapezoid_c_hyp_integral`), with c_hyp read from the run's own
    HyperbolicPropagationSpeeds: a per-step rate the particles actually
    decayed under, not an a-priori bound from margin*h/dt_max. SWIFT
    quantises dt_max DOWN to a power-of-two subdivision of the run's own
    span before any particle uses it (`read_timeline_dt_max`, whose
    docstring states what the log line it reads does and does not report),
    so a bound built from the raw parameter under-states c_hyp ~ h/dt by
    that factor and is not conservative. The coefficients are SIGNED, since the
    residual is two-sided:

        d ln u_LW / d[int (c_hyp/c) d ln a] = -lambda_E(LW)
        d ln u_PE / d[int (c_hyp/c) d ln a] = [lambda_E(LW) - 1] r
                                              - lambda_E(PE) ,

    r the box-mean u_LW/u_PE at the reference snapshot. PE's own second term
    is the LW-to-PE band-edge transfer: radiation_end_force_propagation adds
    [lambda_E(LW) - 1] H_dilated u_LW into PE's `u` every step, on top of
    PE's own -lambda_E(PE) H_dilated u_PE decay. At this fixture's own table
    weights that sum is POSITIVE, so PE's box mean RISES while LW's falls,
    and a gate on |drift| alone cannot tell a rise from a fall of the same
    size.

    The drift is measured from the first snapshot whose own c_hyp is a
    genuine per-step rate, not from snapshot 0: an unpinned run's snapshot 0
    is written before the first force step and still holds the module's
    first-init light-speed clamp (`band_edge_ratio_reference_index`). That
    one interval is therefore dropped rather than integrated with the clamp,
    or with a later snapshot's value substituted for it -- a substitution
    makes the interval's endpoints equal and so credits the interval whose
    speed is least known with a discretisation bound of exactly zero.

    (A1) cannot detect a wrong propagation-speed MAGNITUDE: its prediction
    integrates the same recorded c_hyp the module decayed with, so a
    systematically wrong speed moves the measurement and the prediction
    together. What it does constrain is the band-edge coefficients, the sign
    of each band's drift, and the c_hyp/c dilation of the Hubble term,
    against the run's own speed. Constraining the speed itself needs a check
    that predicts c_hyp from the run's own resolution and step, which this
    one deliberately does not.

    With propagation off (GEARFeedback:ISRF_propagation, read by
    `read_isrf_propagation`) in a COSMOLOGICAL run the predicted drift is
    exactly zero: radiation_end_force_propagation returns before its
    relaxation update, so no cosmological term exists to predict. (A1) and
    (A3) then gate the measured drift against their error-only bar rather
    than skipping, which is still a real gate -- a module that applied the
    decay anyway fails it.

    In a NON-COSMOLOGICAL run the prediction is identically zero, since
    H = 0 leaves no decay for the module to apply, and every cosmological
    term of the bar vanishes with it. What remains is the float-divergence
    floor on its own, a derived bound, the same one (A3) builds. Whether
    that floor bounds the measured drift depends on one thing, and (A1)
    splits the leg on it.

    The ledger is a ratio of two sums taken at the SAME time, so a c_hyp
    that is UNIFORM ACROSS THE BOX cancels between numerator and
    denominator. An unpinned run's c_hyp is not: under scheme 4 it is
    margin*h/dt_max, so it carries the glass's own h spread, and
    the ledger mean then moves for a second reason that has nothing to do
    with conservation. The scheme conserves the ledger's NUMERATOR at
    fixed weights; it conserves neither the numerator once the weights
    move nor the normalised ratio. Per interval the drift splits exactly
    into a REWEIGHT term, mean(u_next, w_next) - mean(u_next, w), and a
    fixed-weight TRANSPORT term, and the first is bounded by
    spread(u)*spread(1/c_hyp) with no conservation content.

    So (A1) at H = 0:

    - With the recorded c_hyp BIT-UNIFORM across the box, every particle's
      c_hyp is one bit-identical value, the weight cancels identically, the
      metric becomes the mass-weighted box mean the pairwise exchange
      conserves, and the float-divergence floor is the whole error budget.
      The leg is GATED against that floor, and against it alone:
      `--reference` cannot inflate it, because the reference term is not in
      this bar. Uniformity is read from the field, not from whichever
      parameter produced it, because it is the condition the cancellation
      rests on: GEARFeedback:ISRF_c_hyp_pin_for_debugging gives it (the same
      mechanism `dust_absorption` uses), and so does ISRF_c_hyp_scheme 2,
      whose GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c is one speed for the
      whole box with the Courant condition imposed on the timestep. The pin
      is applied after the light-speed clamp
      (radiation_isrf.c's radiation_snapshot_part_propagation and
      radiation_end_density_propagation), and the operators are the same
      whether the speed is uniform or not (radiation_propagation_iact.h).
      What a uniform speed does cost is the variable-c coverage: this leg
      says nothing about a defect that only appears once c_hyp varies between
      neighbours.
    - With c_hyp VARYING, the drift is REPORTED and not gated, the way (A3)
      already reports and skips there. The reweighting term is not unbounded:
      spread(u)*spread(1/c_hyp) bounds it, from the run's own recorded
      spreads. What that bound covers, though, is a quantity the scheme
      never promised to conserve, so adding it to the bar would gate
      nothing, while sizing a term from the measured drift instead would
      fit the bar to the data. Neither was done. The consequence is stated
      rather than hidden: between the float floor
      and (A2)'s own bar, the unpinned leg does not bound the transport
      ledger's conservation, and the pinned leg is where that bound lives.
      The drift is printed instead, so it stays visible, and the two
      checks the unpinned leg does carry stay live: the measured drift must
      be finite, and (A2) is gated exactly as it is with cosmology. The
      unpinned non-cosmological run also remains the `--reference` of the
      cosmological one, where its drift enters that run's bar, is
      therefore still acted upon, and is held below the predicted decay by
      the resolution self-test.

    A claimed uniform speed the module did not deliver FAILS rather than
    falling back to the report. Both routes count as a claim: a positive
    ISRF_c_hyp_pin_for_debugging, and ISRF_c_hyp_scheme == 2. If either is set
    and the recorded c_hyp is not bit-uniform, the uniform-weight ledger the
    gate rests on is not live, and a silent fallback would remove the gate
    exactly when something on that path had broken. Note what this covers:
    the parameter validation in feedback_props_init() already rejects a
    scheme/fraction mismatch at start-up, so what reaches here is a module or
    plumbing failure, which is what this leg exists to catch.

    The unshielded H2 photodissociation rate the module hands to Grackle is
    ``k = (sigma_H2/E_LW) c rho u_LW``, with rho = rho0 (a0/a)^3 and u_LW
    carrying its OWN decay. Dilating the Hubble term by c_hyp/c does not remove
    that decay, it scales it: the module relaxes u_LW at depth
    ``lambda_E(LW) H_dilated``, so with no absorption
    ``u_LW = u_LW,0 (a0/a)^beta``, ``beta = lambda_E(LW) c_hyp/c``, and

        ln[x_H2(t)/x_H2(0)] = -k0 int_0^t (a0/a)^(3+beta) dt'                (A2)

    (``-k0 t`` without cosmology). **CORRECTED 2026-10-01, and the correction
    matters:** this reference previously used the exponent 3 alone, on the
    ground that the field "no longer decays at leading order". That is true only
    of a GREY band, ``lambda_E = 1``. Once the band-edge transfer gave the LW
    band ``lambda_E(LW) = 6.42`` at this fixture's own table, the neglected
    factor grew by 6.42 and the check failed a correct run: on the 1000 km/s
    cosmological leg the omission is 2.25e-03 at the last snapshot against a
    7.18e-04 bar, which is the 3.05x that leg reported. The omission is in the
    PREDICTION, so it is corrected there, not absorbed into the bar.

    The exponent beta is formed only when the recorded speed is one value for
    the whole run, where it is exact with no quadrature; a varying speed would
    need the trapezoid and is not attempted, and the omission is printed instead
    of hidden.

    This is the more discriminating of the two: the exponent is 3 + beta
    (density cubed times the field's own decay), a ~10% shift in the
    predicted x_H2 over this fixture's span, well above any noise floor.

    A THIRD, more direct probe of the band-edge transfer itself cancels
    everything the two LW moments share: LW's energy moment and its photon-number
    moment (LWSpecificEnergies/LWPhotonSpecificEnergies, ``u`` and ``n``
    below) share one transport operator and the same per-particle c_hyp and
    dt, but ``radiation_end_force_propagation`` (radiation_isrf.c) relaxes
    them at two DIFFERENT depths, ``lambda_E(LW)*H_dilated`` and
    ``lambda_N(LW)*H_dilated`` (the function's own doxygen: "the diagnostic
    d ln(U/N)/dt identity the photon moment exists for collapses to zero" if
    a shared, per-operator lambda were used instead of the per-moment lookup
    it insists on). Everything else the two moments share -- transport,
    glass noise, c_hyp itself, the float32 upstream inputs -- cancels in
    their ratio, to first order, leaving

        d ln(u/n)/dt = -(lambda_E(LW) - lambda_N(LW)) H_dilated ,

    and, since H dt = d ln a exactly (the scale factor's own definition, no
    approximation), integrating over the run's own snapshots gives

        ln[u(t)/n(t)] - ln[u(0)/n(0)]
            = -(lambda_E(LW) - lambda_N(LW)) int_0^t (c_hyp/c) d ln a .        (A3)

    ``int (c_hyp/c) d ln a`` is evaluated as a trapezoid in ln a over the
    run's own snapshots, from the c_hyp each one's own HyperbolicPropagation-
    Speeds field recorded (median over particles; the field is uniform so
    the spread is glass noise, not signal) -- NOT from a fitted or assumed
    c_hyp, and not from A1's conservative c_hyp/c upper bound, which is too
    loose by orders of magnitude to resolve (A3)'s own signal. The first
    snapshot of an unpinned run is excluded: it is written before the first
    force step, so its own HyperbolicPropagationSpeeds is still the light-
    speed clamp `c_hyp = c` the module falls back on before any step has
    run, not a rate a particle ever actually decayed under.

    lambda_E(LW) and lambda_N(LW) are read from the SAME start-up log line
    A1 already reads (radiation_set_band_edge_coefficients, radiation.c),
    computed from the table before any step runs -- a different code path
    from the per-step force update (A3) is checking, so this is not an
    algebraic identity of the module under test: a shared-lambda bug (the
    one the function's own doxygen warns against) leaves the measured ratio
    flat while (A3)'s prediction stays at its full, nonzero, mechanism-
    derived value.

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
discretisation. Every step size and step count in them is the quantised
dt_max read back from the run's own start-up log (`read_timeline_dt_max`,
which states the one way that number is approximate under cosmology), not
the raw TimeIntegration:dt_max parameter, which over-states the step and
under-states the step count at once, so using it is not conservative in
either direction. The generic terms:

- H and the rates are frozen at the step end (``cosmology_update`` runs before
  the step's tasks). For a rate r(a) ~ a^-p, the ln error per step is
  (p/2) dlna_step * r dt (matter domination, d ln H/d ln a = -3/2 adds 3/4
  for the H term), summed over the run.
- Particles are updated at their step ends, so a snapshot can lag by one
  step: one step's worth of the change.
- The floor that does not depend on a. The specific energy is DOUBLE in the
  struct (feedback_struct.h), in its update and in the snapshot field
  (feedback_io.h), so a term modelling float32 round-off OF it is not a
  legitimate noise floor and none is used. What is still float is each
  moment's own flux divergence, which the double relaxation update
  subtracts every step: bounded by FLOAT32_EPS times the step count times
  that increment's own relative size against u, read from the run's own
  snapshots (`float_divergence_pull`). A paired non-cosmological run of the
  same configuration (``--reference``) MEASURES this floor rather than
  bounding it, and twice its measured drift replaces the bound wherever it
  is larger.

(A2)'s bar is relative to the predicted exponent, like its residual. Its
float term is the species fraction's own storage: H2I is a float rewritten
every step (cooling_struct.h), so ln x_H2 carries at most eps/2 per step
taken up to that snapshot, divided by the exponent predicted up to it.

(A1)'s bar carries no term for the predicted decay itself, which is on the
other side of the residual now, only: the log line's own ``%.5g`` precision
on each lambda, propagated through the band's own coefficient and times the
run's own integral; the trapezoid's own discretisation bound; the two
generic terms above, scaled by that coefficient and by the run's own mean
c_hyp/c; and, for PE alone, the drift of r over the run, which bounds
holding it at its reference value.

(B1) is already a two-sided residual (measured minus predicted ln sum m u),
and its bar does not contain the predicted decay: the terms are the one-step
lag (depth / step count), the drift of the mass-weighted comoving density
that kappa is proportional to (times depth), the float flux-divergence floor
above, and, with cosmology, the step-end kappa and H terms. A build with no
absorption, or with the wrong sign, leaves a residual of 1 to 2 times the
depth, thousands of times the bar. B1 is written for a pinned c_hyp with
lambda = 1: the run's own band-edge weights (a factor lambda_E on the Hubble
term, plus the LW-to-PE transfer) are of order (c_pin/c) ln(a/a0), ~1e-5,
and are not modelled.

(A3)'s bar is built the same way, plus two terms of its own, neither fitted
to a measured residual:

- the log line's own ``%.5g`` precision on lambda_E(LW)/lambda_N(LW) (five
  significant digits), propagated through their difference as a half-ulp
  bound on each, times the run's own int (c_hyp/c) d ln a;
- the trapezoid's own discretisation error from treating c_hyp as piecewise
  constant between snapshots, bounded per interval by half the interval's
  own |Delta c_hyp|/c times its d ln a (a bound on the deviation from the
  true, continuously-varying c_hyp(t), not a correction to it).
"""

import argparse
import atexit
import glob
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np
from scipy.integrate import quad

RADIATION_H = (
    Path(__file__).resolve().parents[5] / "src" / "feedback" / "GEAR" / "radiation.h"
)

ELECTRON_VOLT_CGS = 1.602176634e-12


# Values of the macros below, for the case where the header is out of reach.
# This script is normally copied out of the source tree to run beside a run
# directory, where there is nothing to parse; in the tree the parse is what
# catches a copy that has drifted from the header. Falling back is announced
# on stdout, since a stale copy would otherwise print a confident number.
RADIATION_H_FALLBACK = {
    "RADIATION_SIGMA_H2_OVER_E_LW_CGS": 1.2847106348798106e-07,
    "RADIATION_LW_PHOTON_ENERGY_EV": 12.2,
}
RADIATION_H_FALLBACK_SOURCE = "radiation.h at 18bcf76cb, copied 2026-09-24"
RADIATION_H_FALLBACK_USED: List[str] = []


def _warn_fallback_used() -> None:
    """Repeat the unverified-constant note next to the verdict."""
    if RADIATION_H_FALLBACK_USED:
        print(
            "NOTE: "
            + ", ".join(RADIATION_H_FALLBACK_USED)
            + f" came from this script's own copy of {RADIATION_H_FALLBACK_SOURCE}, "
            "not from the header, and were not verified against the source. "
            "Do not cite the numbers above as checked against the code "
            "without confirming the header still carries these values."
        )


atexit.register(_warn_fallback_used)


def read_radiation_h_constant(name: str) -> float:
    """Read a #define'd float constant's value out of radiation.h.

    Parameters
    ----------
    name : str
        The macro name, e.g. ``"RADIATION_SIGMA_H2_OVER_E_LW_CGS"``.

    Returns
    -------
    float
        The macro's value, or its entry in `RADIATION_H_FALLBACK` when the
        header is not reachable from this file. The fallback prints a note
        on stdout naming the constant and the unverified value it used.
    """
    try:
        text = RADIATION_H.read_text()
    except OSError:
        value = RADIATION_H_FALLBACK[name]
        RADIATION_H_FALLBACK_USED.append(name)
        print(
            f"NOTE: {RADIATION_H} is not readable from this copy of the "
            f"script, so {name} = {value!r} is taken from the script's own "
            f"copy of the header value ({RADIATION_H_FALLBACK_SOURCE}) and "
            f"was NOT verified against the source. If the header has "
            f"changed since, every number below is stale."
        )
        return value
    match = re.search(rf"^#define\s+{re.escape(name)}\s+([0-9.eE+-]+)", text, re.M)
    if match is None:
        raise ValueError(f"Could not find #define {name} in {RADIATION_H}")
    return float(match.group(1))


LW_CALIBRATION_NOTE: "List[str]" = []

# The two lines radiation_set_lw_photon_energy_cgs() prints, anchored on the
# function name error.h's message() macro prefixes them with. A bare
# "sigma_H2/E_LW" would also match this script's own note below, so a log
# holding an earlier run of this script would classify itself as clean.
LW_QUOTIENT_LOG_RE = re.compile(
    r"radiation_set_lw_photon_energy_cgs: H2 photodissociation sigma_H2/E_LW ="
)
LW_TABLE_ENERGY_LOG_RE = re.compile(
    r"radiation_set_lw_photon_energy_cgs: Mean Lyman-Werner photon energy "
    r"from the table = (\S+) erg"
)
# engine_config()'s own quantisation of TimeIntegration:dt_max: it halves
# (time_end - time_begin) until the result is at most dt_max, so the raw
# parameter over-states the step size and under-states the step count.
# Under cosmology the halved span is PROPER TIME while dt_max is in d ln a,
# so the printed number is not the time-line's own step: see
# `read_timeline_dt_max`.
TIMELINE_DT_MAX_LOG_RE = re.compile(
    r"engine_config: Maximal timestep size \(on time-line\): (\S+)"
)
# radiation_set_band_edge_coefficients()'s own start-up announcement
# (radiation.c): the table-derived lambda_E(PE)/lambda_E(LW)/lambda_N(LW)
# this run actually propagated with.
BAND_EDGE_WEIGHTS_LOG_RE = re.compile(
    r"radiation_set_band_edge_coefficients: Band-edge weights .*"
    r"lambda_E\(PE\)=([-+0-9.eE]+), lambda_E\(LW\)=([-+0-9.eE]+), "
    r"lambda_N\(LW\)=([-+0-9.eE]+)"
)


def _note_lw_calibration(message: str) -> None:
    """Print a calibration note now, and again next to the verdict."""
    print("NOTE: " + message)
    LW_CALIBRATION_NOTE.append(message)


def _repeat_lw_calibration_notes() -> None:
    """Repeat the calibration notes next to the verdict."""
    for message in LW_CALIBRATION_NOTE:
        print("NOTE: " + message)


atexit.register(_repeat_lw_calibration_notes)


def find_run_logs(
    *snapshot_globs: "Optional[str]", log: "Optional[str]" = None
) -> "List[Path]":
    """Find the ``output.log`` of every run this check was pointed at.

    The example ``run.sh`` scripts move the log into the per-run directory
    beside ``snap``, and the READMEs then call this script from the example
    directory with a snapshot glob, so the log sits one level above the
    snapshots rather than in the working directory.

    Parameters
    ----------
    snapshot_globs : str, optional
        Snapshot globs of the runs being checked. None entries are ignored,
        so an unset option can be passed straight through.
    log : str, optional
        A run log named explicitly on the command line. Taken as given,
        whatever it is called.

    Returns
    -------
    list of pathlib.Path
        Existing logs, deduplicated, in the order they were named. The
        working directory is consulted only when nothing else named a log,
        so checking one run from inside another run's directory cannot
        report the wrong run.
    """
    candidates: "List[Path]" = []
    if log is not None:
        candidates.append(Path(log))
    for hint in snapshot_globs:
        if hint is None:
            continue
        given = Path(hint)
        candidates += [given.parent / "output.log", given.parent.parent / "output.log"]
    if not candidates:
        candidates.append(Path.cwd() / "output.log")
    logs: "List[Path]" = []
    seen = set()
    for candidate in candidates:
        if not candidate.is_file():
            continue
        key = candidate.resolve()
        if key in seen:
            continue
        seen.add(key)
        logs.append(candidate)
    return logs


def read_band_edge_weights(
    *snapshot_globs: "Optional[str]", log: "Optional[str]" = None
) -> "Optional[Dict[str, float]]":
    """Return this run's own table-derived band-edge weights, from its log.

    ``radiation_set_band_edge_coefficients()`` announces
    ``lambda_E(PE)``/``lambda_E(LW)``/``lambda_N(LW)`` once at start-up,
    computed from the same radiation table the run propagated with
    (``radiation.c``). Reading them here, rather than hardcoding a number,
    keeps the cosmological allowance in `check_free_field` correct when the
    table changes; the run's own log is the record of what it actually used.

    Parameters
    ----------
    snapshot_globs : str, optional
        Snapshot globs of the run being checked.
    log : str, optional
        A run log named explicitly on the command line.

    Returns
    -------
    dict or None
        ``{"PE": lambda_E(PE), "LW": lambda_E(LW), "N_LW": lambda_N(LW)}``,
        or None if no log announced them.
    """
    for path in find_run_logs(*snapshot_globs, log=log):
        match = BAND_EDGE_WEIGHTS_LOG_RE.search(path.read_text(errors="replace"))
        if match:
            return {
                "PE": float(match.group(1)),
                "LW": float(match.group(2)),
                "N_LW": float(match.group(3)),
            }
    return None


def check_run_lw_calibration(
    *snapshot_globs: "Optional[str]", log: "Optional[str]" = None
) -> None:
    """Report which H2 calibration each run's own binary used, from its log.

    ``radiation_set_lw_photon_energy_cgs`` announces the coefficient in force
    on every run whose radiation model is active, so the run's log is the
    witness, not this source tree and not the radiation table. Three binaries
    exist and the log tells them apart:

    * one printing ``sigma_H2/E_LW``: it multiplies the LW energy flux by the
      Sternberg-anchored quotient, which is what `unshielded_rate` reproduces.
      No bias is possible;
    * one printing only the table's mean LW photon energy: it divided the flux
      by that value while the cross section stayed pinned to 12.2 eV, so its
      rates sit below this script's by the ratio of the two energies. The
      offset is quantified here from the logged value;
    * one printing neither: it divided by the header constant, which is the
      same rate the quotient gives.

    A missing log leaves the question open and is reported as such. The table
    is deliberately not consulted: it says what a binary COULD have read, not
    what it did.

    Parameters
    ----------
    snapshot_globs : str, optional
        Snapshot globs of the runs being checked.
    log : str, optional
        A run log named explicitly on the command line.
    """
    logs = find_run_logs(*snapshot_globs, log=log)
    if not logs:
        _note_lw_calibration(
            "no output.log was found beside the runs being checked, so the H2 "
            "calibration their binaries used is UNKNOWN. A binary that divided "
            "the LW flux by the radiation table's mean photon energy produces "
            "rates about 0.4 per cent below the ones predicted here, which "
            "shows up as a residual of that size."
        )
        return
    for log in logs:
        text = log.read_text(errors="replace")
        if LW_QUOTIENT_LOG_RE.search(text):
            print(
                f"H2 calibration: {log} reports sigma_H2/E_LW, so that run used "
                "the same Sternberg-anchored quotient as this script."
            )
            continue
        match = LW_TABLE_ENERGY_LOG_RE.search(text)
        if match is None:
            print(
                f"H2 calibration: {log} reports no photon-energy line, so that "
                "run divided by RADIATION_LW_PHOTON_ENERGY_EV, which gives the "
                "same rate as the quotient this script uses."
            )
            continue
        energy_logged = float(match.group(1))
        energy_ev = read_radiation_h_constant("RADIATION_LW_PHOTON_ENERGY_EV")
        offset = energy_ev * ELECTRON_VOLT_CGS / energy_logged - 1.0
        _note_lw_calibration(
            f"{log} reports a table mean LW photon energy of "
            f"{energy_logged:.6e} erg and no sigma_H2/E_LW line, so that run's "
            "binary DIVIDED the LW flux by that energy against a cross section "
            f"pinned at {energy_ev:.1f} eV. Its H2 photodissociation rates "
            f"differ from the ones predicted here by {offset * 100.0:+.2f} per "
            "cent, and a residual of that size is the expected symptom, not a "
            "physics result."
        )


# The H2 photodissociation rate is the LW energy flux times this one
# Sternberg-anchored quotient; no cross section or photon energy enters it.
SIGMA_H2_OVER_E_LW_CGS = read_radiation_h_constant("RADIATION_SIGMA_H2_OVER_E_LW_CGS")
HABING_FLUX_CGS = 1.6e-3
SIGMA_D_CGS = {"PE": 9e-22, "LW": 1.5e-21}
GRACKLE_DEFAULT_DUST_TO_GAS_RATIO = 0.009387
# Mirrors the C code: the extinction chain takes the proton mass from the
# physical constants (const_proton_mass_cgs in physical_constants_cgs.h).
EXTINCTION_PROTON_MASS_CGS = 1.67262192369e-24
KERNEL_GAMMA_DEFAULT = 1.936492
MU_H = 1.4
GRACKLE_SOLAR_METAL_FRACTION = 0.01295
C_LIGHT_CGS = 2.99792458e10
M_H_CGS = 1.67262171e-24
HYDROGEN_MASS_FRACTION = 0.76
PHOTOELECTRIC_RATE_CGS = 1e-24 * 0.05
# float(), so every bar term built from it is evaluated in double. Left as
# the numpy float32 scalar, an expression such as n * eps / (1 - n * eps)
# rounds to float32 at each step, which is the precision of the quantity the
# bar is meant to bound.
FLOAT32_EPS = float(np.finfo(np.float32).eps)
# The step line prints its step-size field with "%14e" (src/engine.c), i.e.
# six decimals of mantissa, so a step size read back from the log carries
# half a unit in that last decimal.
LOG_STEP_MANTISSA_HALF_ULP = 0.5e-6
# enum isrf_c_hyp_scheme values whose c_hyp varies between particles, so that
# the conserved sum m u / c_hyp is not proportional to sum m u
# (feedback_properties.h). Scheme 2's operators conserve the same sum, but its
# c_hyp is one value for the whole box.
VARIABLE_C_SCHEMES = (4,)
# enum isrf_c_hyp_scheme values the module accepts (feedback_properties.h).
VALID_C_HYP_SCHEMES = (2, 4)
# HyperbolicPropagationSpeeds reads exactly c on a snapshot written before
# the first force step (the module's first-init clamp): below this fraction
# of c, (A3) takes it as a genuine per-step value instead.
C_HYP_CLAMP_FRACTION_OF_C = 0.999


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
        default=None,
        help="injection_dusty only: max allowed max_j |R_j| on the band-ratio "
        "gate. Omitted, the bar is DERIVED from the run's own signal and the "
        "float32 widths of the reconstruction, and printed with its terms. "
        "Pass a positive value to override it, or a negative one to report the "
        "residual without gating on it.",
    )
    parser.add_argument(
        "--kernel-gamma",
        type=float,
        default=KERNEL_GAMMA_DEFAULT,
        help="injection_dusty only: the binary's kernel_gamma (default: "
        "%(default)s, Wendland C2 in 3D). Only the constant_kernel_path "
        "mechanism builds its column from it.",
    )
    parser.add_argument(
        "--extinction-path",
        type=float,
        default=None,
        help="injection_dusty only: force the constant_kernel_path mirror at "
        "this many kernel support radii, whatever mechanism the run "
        "recorded. The mechanism and its length are read from the run's "
        "used_parameters.yml when this is omitted.",
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
            # LW's photon-number moment (energy-equivalent units), used only by (A3).
            "n_LW": (
                physical(gas["LWPhotonSpecificEnergies"], a, energy)[order]
                if "LWPhotonSpecificEnergies" in gas
                else None
            ),
            # (A1)'s and (A3)'s float-residual bar terms: the FLOAT inputs
            # still feeding the double relaxation update
            # (radiation_end_force_propagation, radiation_isrf.c:913-917).
            "div_PE": (
                physical(gas["PESpecificFluxDivergences"], a, energy / time)[order]
                if "PESpecificFluxDivergences" in gas
                else None
            ),
            "div_LW": (
                physical(gas["LWSpecificFluxDivergences"], a, energy / time)[order]
                if "LWSpecificFluxDivergences" in gas
                else None
            ),
            "div_LW_photon": (
                physical(gas["LWPhotonSpecificFluxDivergences"], a, energy / time)[
                    order
                ]
                if "LWPhotonSpecificFluxDivergences" in gas
                else None
            ),
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
        # Comoving code units, for the pair_separation column only: the
        # minimum image is taken in those units against the header box and
        # scaled once, coordinates carrying an a-scale exponent of 1.
        out["pos_comoving"] = gas["Coordinates"][:].astype(np.float64)[order]
        out["box_comoving"] = np.atleast_1d(header["BoxSize"]).astype(np.float64)
        out["length_unit"] = length
        if "/PartType4" in handle and handle["/PartType4/Masses"].shape[0] > 0:
            out["L_PE"] = float(handle["/PartType4/PELuminosities"][0])
            out["L_LW"] = float(handle["/PartType4/LWLuminosities"][0])
            out["n_stars"] = int(handle["/PartType4/Masses"].shape[0])
            out["star_comoving"] = handle["/PartType4/Coordinates"][:][0].astype(
                np.float64
            )
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


def read_extinction_path(pattern: str, given: Optional[float]) -> Tuple[str, float]:
    """Return the extinction path mechanism and its kernel-radii multiple.

    The multiple is meaningful for constant_kernel_path only; the other
    mechanisms carry their own length and report it as nan.

    Parameters
    ----------
    pattern
        Snapshot glob, used to locate the run's used_parameters.yml.
    given
        --extinction-path, which forces the constant_kernel_path mirror at
        that many kernel support radii whatever the run recorded.

    Returns
    -------
    mechanism, path_in_kernel_radii
        The mechanism name and its kernel-radii multiple.
    """
    if given is not None:
        return "constant_kernel_path", given
    import os
    import yaml

    directory = os.path.dirname(os.path.dirname(sorted(glob.glob(pattern))[0]))
    with open(os.path.join(directory, "used_parameters.yml")) as handle:
        used = yaml.safe_load(handle)["GEARFeedback"]
    if "ISRF_extinction_path" not in used:
        # A run archived before the key existed recorded no value at all, and
        # the path then in force was two kernel support radii.
        return "constant_kernel_path", 2.0
    name = used["ISRF_extinction_path"]
    if name == "constant_kernel_path":
        return name, float(used["ISRF_extinction_path_in_kernel_radii"])
    if name == "pair_separation":
        return name, float("nan")
    if name == "temperature_capped_jeans":
        raise RuntimeError(
            f"GEARFeedback:ISRF_extinction_path {name!r} builds the column "
            "from the gas state through cooling_get_temperature, which this "
            "check cannot reconstruct from a snapshot. Rerun the fixture "
            "with constant_kernel_path or pair_separation, or extend this "
            "check to that mechanism's own length."
        )
    raise RuntimeError(f"Unknown GEARFeedback:ISRF_extinction_path {name!r}")


def extinction_path_cgs(
    snap: Dict, mechanism: str, path_in_kernel_radii: float, kernel_gamma: float
) -> np.ndarray:
    """Return each gas particle's physical extinction path, cm.

    Parameters
    ----------
    snap
        One snapshot as returned by read_snapshot.
    mechanism
        GEARFeedback:ISRF_extinction_path, as read_extinction_path reports it.
    path_in_kernel_radii
        Kernel-radii multiple of the constant_kernel_path mechanism.
    kernel_gamma
        The binary's kernel_gamma.

    Returns
    -------
    numpy.ndarray
        Physical extinction path of each particle, cm.
    """
    if mechanism != "pair_separation":
        return path_in_kernel_radii * kernel_gamma * snap["h"]
    if snap.get("n_stars", 0) != 1:
        raise RuntimeError(
            "The pair_separation column is a per-pair length. It is a "
            "per-particle one, and these identities hold, only while a "
            f"single star illuminates the box; this snapshot has "
            f"{snap.get('n_stars', 0)}."
        )
    delta = snap["pos_comoving"] - snap["star_comoving"]
    box = snap["box_comoving"]
    delta -= box * np.round(delta / box)
    separation = np.sqrt((delta * delta).sum(axis=1))
    return separation * snap["a"] * snap["length_unit"]


def optical_depths(
    snap: Dict, path_cgs: np.ndarray, dust_to_gas: float
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
    path_cgs
        Physical extinction path of each particle, cm, as the run's own
        GEARFeedback:ISRF_extinction_path mechanism builds it.
    dust_to_gas
        The run's resolved Grackle chemistry_data.local_dust_to_gas_ratio.

    Returns
    -------
    dict
        The dimensionless optical depth of each particle, keyed by band.
    """
    column = path_cgs * snap["density"]
    d_relative = (
        np.maximum(snap["Z_smoothed"].astype(np.float64), 0.0)
        / GRACKLE_SOLAR_METAL_FRACTION
        * (dust_to_gas / GRACKLE_DEFAULT_DUST_TO_GAS_RATIO)
    )
    prefactor = d_relative / (MU_H * EXTINCTION_PROTON_MASS_CGS) * column
    return {band: SIGMA_D_CGS[band] * prefactor for band in ("PE", "LW")}


def read_isrf_propagation(pattern: str) -> bool:
    """Return whether GEARFeedback:ISRF_propagation was on for this run.

    `radiation_end_force_propagation` (radiation_isrf.c) returns before its
    own relaxation update when propagation is off, so no cosmological term
    exists to predict on such a run: (A1) and (A3) then gate the measured
    drift against their error-only bar with a zero prediction, rather than
    skipping. Missing from used_parameters.yml means the run predates the
    parameter, when propagation was unconditional.

    Parameters
    ----------
    pattern : str
        Snapshot glob of the run.

    Returns
    -------
    bool
        True when propagation was on, or when the parameter is not recorded.
    """
    import os
    import yaml

    directory = os.path.dirname(os.path.dirname(sorted(glob.glob(pattern))[0]))
    path = os.path.join(directory, "used_parameters.yml")
    if not os.path.exists(path):
        return True
    with open(path) as handle:
        parameters = yaml.safe_load(handle)
    try:
        return int(parameters["GEARFeedback"]["ISRF_propagation"]) != 0
    except (KeyError, TypeError, ValueError):
        return True


def read_timeline_dt_max(
    *snapshot_globs: "Optional[str]", log: "Optional[str]" = None
) -> "Optional[float]":
    """Return the quantised dt_max SWIFT announced at start-up.

    SWIFT quantises TimeIntegration:dt_max down to a power-of-two
    subdivision of the run's own span before any particle uses it, and
    announces one such quantisation at start-up (`engine_config.c:675-679`,
    which halves `time_end - time_begin` until the result is at most
    `dt_max`). Every bar term built from a step size or a step count needs a
    quantised value, not the raw parameter: the raw parameter over-states
    the step (inflating a one-step lag term) and under-states the step count
    (shrinking a term proportional to it), in opposite directions, so using
    it is not conservative either way.

    WITHOUT cosmology the announced number IS the time-line's own step.
    WITH cosmology it is not, and this function does not claim it is: the
    halved span `time_end - time_begin` is PROPER TIME (`engine.c:3752-3753`
    copies it from the cosmology model) while `dt_max` and the time-line
    itself are in d ln a (`cosmology.c:915` builds `time_base` from
    `log_a_end - log_a_begin`), so the message halves one quantity against a
    threshold in the other. WHEN both spans exceed `dt_max`, which is the
    only case the halving loop acts in, the true time-line step and this one
    both lie in `(dt_max/2, dt_max]` and so differ by less than a factor of
    two in either direction. A span at or below `dt_max` is halved zero
    times and carries no such bound, and under cosmology nothing rejects
    that case: `engine_config.c:688-691` exempts a cosmological run from the
    `dt_max > span` error. On this example's own z9 fixture both spans
    exceed `dt_max` by 2000x and the difference is 0.14% (5.440215e-05
    printed against ln(0.125/0.1)/4096 = 5.447841e-05).
    Every bar term here is linear in this value or in the step count derived
    from it, so that bounded factor carries straight into the bar and
    nowhere else: no residual and no prediction reads it.

    Parameters
    ----------
    snapshot_globs : str, optional
        Snapshot globs of the run being checked.
    log : str, optional
        A run log named explicitly on the command line.

    Returns
    -------
    float or None
        The quantised dt_max in the same units as the parameter (d ln a
        under cosmology), or None if no log announced it. Under cosmology
        this is the proper-time span quantised against that parameter, not
        the time-line's own step; see above.
    """
    for path in find_run_logs(*snapshot_globs, log=log):
        match = TIMELINE_DT_MAX_LOG_RE.search(path.read_text(errors="replace"))
        if match:
            return float(match.group(1))
    return None


def read_c_hyp_scheme(pattern: str) -> Optional[int]:
    """Return GEARFeedback:ISRF_c_hyp_scheme, or None when it is not recorded.

    Raises
    ------
    ValueError
        If a value is recorded and is not 2 or 4: the run used a removed
        scheme (0, 1 or 3), whose ledger this check cannot interpret.
    """
    import os
    import yaml

    directory = os.path.dirname(os.path.dirname(sorted(glob.glob(pattern))[0]))
    path = os.path.join(directory, "used_parameters.yml")
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        parameters = yaml.safe_load(handle)
    try:
        recorded = parameters["GEARFeedback"]["ISRF_c_hyp_scheme"]
    except (KeyError, TypeError):
        return None
    try:
        scheme = int(recorded)
    except (TypeError, ValueError):
        scheme = None
    if scheme not in VALID_C_HYP_SCHEMES:
        raise ValueError(
            f"{path}: GEARFeedback:ISRF_c_hyp_scheme is {recorded!r}; only 2 "
            "(fixed fraction of c) and 4 (kernel-local speed) are supported. "
            "The values 0, 1 and 3 were removed."
        )
    return scheme


def read_c_hyp_pin(pattern: str) -> Optional[float]:
    """Return GEARFeedback:ISRF_c_hyp_pin_for_debugging, or None if unrecorded.

    The parameter is a propagation speed in the run's own internal velocity
    units, and 0 disables the pin
    (`feedback_properties.h`'s ISRF_c_hyp_pin_for_debugging). It is read
    here only to decide whether the run CLAIMS a pinned speed: whether the
    module actually applied one is established from the snapshots
    themselves (`c_hyp_spatial_spread`).

    Parameters
    ----------
    pattern : str
        Snapshot glob of the run.

    Returns
    -------
    float or None
        The parameter's value, or None when used_parameters.yml is absent
        or does not carry the key.
    """
    import os
    import yaml

    directory = os.path.dirname(os.path.dirname(sorted(glob.glob(pattern))[0]))
    path = os.path.join(directory, "used_parameters.yml")
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        parameters = yaml.safe_load(handle)
    try:
        return float(parameters["GEARFeedback"]["ISRF_c_hyp_pin_for_debugging"])
    except (KeyError, TypeError, ValueError):
        return None


def c_hyp_spatial_spread(run: List[Dict], start: int = 0) -> Optional[float]:
    """Return the worst per-snapshot spatial spread of c_hyp, or None.

    (A1)'s ledger is a ratio of two sums taken at the SAME time, so a
    c_hyp that is uniform ACROSS THE BOX cancels between numerator and
    denominator whatever it does from one snapshot to the next (this
    module's docstring). What the gated leg needs is therefore this
    spatial spread, per snapshot, and exactly zero: a bit-uniform speed
    makes the weighted ledger and the mass-weighted one the same
    functional, whether a debug pin or scheme 2's fixed fraction produced
    it.

    Parameters
    ----------
    run : list of dict
        The run's snapshots, in time order, as `load_run` returns them.
    start : int, optional
        First snapshot to include.

    Returns
    -------
    float or None
        max over snapshots of (max c_hyp - min c_hyp)/max c_hyp, or None
        when any snapshot has no HyperbolicPropagationSpeeds, or one that is
        not everywhere finite, or a non-positive value in any snapshot other
        than an all-zero snapshot 0, or when the run's only snapshot is an
        all-zero snapshot 0. Snapshot 0 alone is skipped when it is all zero:
        see the body.
    """
    seen = False
    worst = 0.0
    for index, snap in enumerate(run[start:], start=start):
        c_hyp = snap["c_hyp"]
        if c_hyp is None:
            return None
        if not np.all(np.isfinite(c_hyp)):
            return None
        if index == 0 and np.all(c_hyp == 0.0):
            # Snapshot 0 is written before the first force step. Under the
            # scheme that sets c_hyp in radiation_snapshot_part_propagation
            # (2) it therefore still holds the first-init seed of exactly
            # zero (radiation_isrf.c:156); the scheme that sets it in
            # radiation_end_density_propagation (4) reads the light-speed
            # value instead, because the initial density pass does
            # run. A snapshot with no speed at all carries nothing to be
            # uniform, so it is skipped. Only index 0 qualifies: after a step
            # has run, an all-zero c_hyp is a degraded run and is rejected
            # below.
            continue
        if np.any(c_hyp <= 0.0):
            return None
        seen = True
        worst = max(worst, float((np.max(c_hyp) - np.min(c_hyp)) / np.max(c_hyp)))
    return worst if seen else None


def use_c_hyp_ledger(run: List[Dict], pattern: str, label: str) -> bool:
    """Report whether the c_hyp-weighted ledger applies to this run.

    It applies only to the scheme whose c_hyp varies between particles (4),
    and only when every snapshot carries a finite, strictly positive
    HyperbolicPropagationSpeeds. Scheme 2 has one speed for the whole box,
    where the weighted and unweighted ledgers are proportional.
    Snapshot 0 is NOT excused here, unlike in `c_hyp_spatial_spread`: this
    weight divides by c_hyp on every snapshot the ledger is evaluated on, and
    for a non-cosmological run that includes snapshot 0. Every rejection
    prints why, so a degraded run is never silently gated on the wrong
    invariant.
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
    """Print a pass/fail line.

    Both operands are tested for finiteness before the comparison: a NaN
    compares false against any bar, and an infinite bar would otherwise
    admit any residual.
    """
    ok = bool(np.isfinite(worst)) and bool(np.isfinite(bar)) and worst <= bar
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


def free_field_errors(
    run: List[Dict],
    use_c_hyp: bool = False,
    ref_index: int = 0,
    h2_field_decay: float = 0.0,
) -> Dict:
    """Return the measured drift of Eq. (A1) and the error of Eq. (A2).

    The transport moves energy between particles and conserves the ledger of
    the run's own c_hyp scheme (this module's docstring), so on a glass each
    particle's field departs from the uniform solution by the glass noise
    while the box mean follows (A1) exactly. The gates use the box means; the
    per-particle spread is reported.

    ``out[band]`` is the MEASURED quantity of (A1)'s two-sided residual,
    ``ln[Q(t)/Q(t_ref)]``, not an error: the caller subtracts (A1)'s own
    predicted drift from it (`check_free_field`). Entries before
    ``ref_index`` are NaN, so a caller cannot read a value the reference
    snapshot does not define.

    Parameters
    ----------
    run
        The run's snapshots, in time order.
    use_c_hyp
        Weight each particle by ``m_i/c_hyp,i`` rather than ``m_i``, for the
        scheme with a per-particle speed (4). Decided by `use_c_hyp_ledger`.
    ref_index
        Snapshot the log drift is measured from. (A1)'s reference is the
        first snapshot whose c_hyp is a genuine per-step rate, since its
        prediction integrates that same recorded c_hyp
        (`band_edge_ratio_reference_index`).
    """
    first = run[0]
    mass = first["mass"]
    out = {"ref_index": ref_index}
    for band in ["PE", "LW"]:
        q_ref = ledger_mean(run[ref_index], f"u_{band}", use_c_hyp)
        drift = np.full(len(run), np.nan)
        for i in range(ref_index, len(run)):
            drift[i] = np.log(ledger_mean(run[i], f"u_{band}", use_c_hyp) / q_ref)
        out[band] = drift
        out[f"{band}_spread"] = np.array(
            [np.median(np.abs(s[f"u_{band}"] / first[f"u_{band}"] - 1.0)) for s in run]
        )
    # Box-mean u_LW/u_PE, the scale factor of the LW-to-PE band-edge
    # transfer in (A1)'s PE prediction. Its own drift over the run bounds
    # the error of holding it at its reference value (`check_free_field`).
    out["r_lw_over_pe"] = np.array(
        [
            float(np.sum(s["mass"] * s["u_LW"]) / np.sum(s["mass"] * s["u_PE"]))
            for s in run
        ]
    )
    # Box-mean density and field: the closed form is for the uniform state.
    rho0 = np.sum(mass) / np.sum(mass / first["density"])
    u_lw0 = np.sum(mass * first["u_LW"]) / np.sum(mass)
    k0 = SIGMA_H2_OVER_E_LW_CGS * C_LIGHT_CGS * rho0 * u_lw0
    # rho ~ (a0/a)^3, times the LW field's OWN decay. Dilating the Hubble term
    # by c_hyp/c does not remove that decay, it scales it: the module relaxes
    # u_LW at depth lambda_E(LW)*H_dilated (radiation_isrf.c's relaxation
    # depth, whose `lambda(m) = 1` grey case is the only one in which the field
    # is constant at leading order), so u_LW ~ (a0/a)^beta with
    # beta = lambda_E(LW)*c_hyp/c, and the integrand picks up 3 + beta. The
    # caller passes beta, or 0 when it cannot be formed exactly.
    integral = power_integral(run, 3.0 + h2_field_decay)
    measured = np.array([np.mean(np.log(s["H2I"] / first["H2I"])) for s in run])
    predicted = -k0 * integral
    out["H2"] = (measured[1:] - predicted[1:]) / np.abs(predicted[1:])
    out["exponent"] = float(-predicted[-1])
    out["rate"] = k0
    out["integral"] = integral
    # Elapsed physical time of every snapshot, seconds, so a reference run on
    # its own time grid can be paired with this one by TIME.
    out["times"] = np.array([s["time"] - first["time"] for s in run])
    return out


def band_edge_weight_log_precision(value: float) -> float:
    """Return half the last-digit step of a value printed with ``%.5g``.

    ``radiation_set_band_edge_coefficients`` (radiation.c) announces each
    lambda with five significant digits, so the run's own value is known
    only to within half of that last printed digit's step. Never used to
    round a value, only to bound (A3)'s log-quantisation bar term.

    Parameters
    ----------
    value : float
        A lambda value as `read_band_edge_weights` parsed it from the log.

    Returns
    -------
    float
        Half the absolute step of its fifth significant digit, or 0 for
        ``value == 0``.
    """
    if value == 0.0:
        return 0.0
    exponent = np.floor(np.log10(abs(value)))
    return 0.5 * 10.0 ** (exponent - 4.0)


def band_edge_ratio_reference_index(run: List[Dict]) -> Optional[int]:
    """Return the first snapshot whose c_hyp is a genuine per-step value.

    Snapshot 0 of a run is written before the first force step, so its own
    HyperbolicPropagationSpeeds carries no rate any particle actually decayed
    under; (A3) needs the trajectory the particles actually experienced, so it
    starts integrating one snapshot later. Two shapes of that pre-step
    snapshot exist and both must be rejected. Under the scheme that sets
    ``c_hyp`` in the density loop (4) it reads the module's first-init
    light-speed clamp (``c_hyp = c``), because that loop does run before the
    first snapshot is written. Under the scheme that sets it in the snapshot
    hook (2) it reads exactly zero for every particle; that is MEASURED for
    scheme 2 (both cluster legs of 2026-09-30, 32768 particles, one distinct
    float32 value, 0). A median test alone passes the zero shape, because zero is below the clamp,
    so strict positivity is required as well.

    Parameters
    ----------
    run : list of dict
        The run's snapshots, in time order, as `load_run` returns them.

    Returns
    -------
    int or None
        Index of the first usable snapshot, or None if none qualifies.
    """
    for i, snap in enumerate(run):
        c_hyp = snap["c_hyp"]
        if (
            c_hyp is not None
            and np.all(np.isfinite(c_hyp))
            and np.all(c_hyp > 0.0)
            and np.median(c_hyp) < C_HYP_CLAMP_FRACTION_OF_C * C_LIGHT_CGS
        ):
            return i
    return None


def trapezoid_c_hyp_integral(
    run: List[Dict], c_hyp_at: List[float], ref_index: int
) -> Dict:
    """Return the trapezoid ``int (c_hyp/c) d ln a`` and its own error bound.

    Shared by (A1) and (A3), which need the same construction for two
    different prefactors. ``H dt = d ln a`` exactly (the scale factor's own
    definition, no approximation), so the dilated integral is a trapezoid in
    ln a over the snapshots the run itself wrote. Its discretisation error,
    from holding c_hyp piecewise constant between snapshots rather than at
    its true, continuously-varying value, is bounded per interval by half
    the interval's own |Delta c_hyp|/c times its d ln a (a bound on the
    deviation, not a signed correction to it), accumulated in absolute
    value.

    The integral starts at ``ref_index`` and is zero at and before it: no
    interval is ever credited to a snapshot whose own c_hyp is not a
    genuine per-step rate, and none is carried with substituted endpoints,
    which would give the one interval whose speed is least known a
    discretisation bound of exactly zero.

    Parameters
    ----------
    run : list of dict
        The run's snapshots, in time order.
    c_hyp_at : list of float
        One representative c_hyp per snapshot, physical cgs (the median over
        particles: the field is uniform, so the spread is glass noise). All
        zero for a run whose module applied no dilated Hubble term at all.
    ref_index : int
        Snapshot index where the accumulated integral is defined to be zero.

    Returns
    -------
    dict
        ``dilated_lna`` (the cumulative trapezoid, one entry per snapshot),
        ``trapezoid_bound`` (its cumulative discretisation bound) and
        ``span`` (total ln a from ``ref_index`` to the run's end).
    """
    dilated_lna = np.zeros(len(run))
    trapezoid_bound = np.zeros(len(run))
    for i in range(ref_index + 1, len(run)):
        dlna = float(np.log(run[i]["a"] / run[i - 1]["a"]))
        c_prev, c_here = c_hyp_at[i - 1], c_hyp_at[i]
        dilated_lna[i] = (
            dilated_lna[i - 1] + 0.5 * (c_prev + c_here) / C_LIGHT_CGS * dlna
        )
        trapezoid_bound[i] = (
            trapezoid_bound[i - 1] + 0.5 * abs(c_here - c_prev) / C_LIGHT_CGS * dlna
        )
    return {
        "dilated_lna": dilated_lna,
        "trapezoid_bound": trapezoid_bound,
        "span": float(np.log(run[-1]["a"] / run[ref_index]["a"])),
    }


def representative_c_hyp(run: List[Dict], propagation_on: bool) -> List[float]:
    """Return one c_hyp per snapshot for `trapezoid_c_hyp_integral`.

    All zero when propagation is off: `radiation_end_force_propagation`
    returns before its relaxation update then, so the module applies no
    dilated Hubble term and the predicted drift is exactly zero, whatever
    the snapshots' own HyperbolicPropagationSpeeds field happens to hold
    (the module's first-init light-speed clamp, in that case).

    Parameters
    ----------
    run : list of dict
        The run's snapshots, in time order.
    propagation_on : bool
        GEARFeedback:ISRF_propagation of the run, from
        `read_isrf_propagation`.

    Returns
    -------
    list of float
        Median c_hyp per snapshot, physical cgs, or zeros.
    """
    if not propagation_on:
        return [0.0] * len(run)
    return [float(np.median(s["c_hyp"])) for s in run]


def band_edge_ratio_errors(
    run: List[Dict], start: int, propagation_on: bool = True
) -> Dict:
    """Return (A3)'s measured ln(u_LW/n_LW) drift and its dilated-ln(a) integral.

    Both series start at `run[start]` (see `band_edge_ratio_reference_index`)
    and are indexed like `run` itself, with every entry before `start` left
    as NaN so a caller cannot silently read a meaningless value. The dilated
    integral is `trapezoid_c_hyp_integral`, shared with (A1), and is
    identically zero when propagation is off.

    Parameters
    ----------
    run : list of dict
        The run's snapshots, in time order.
    start : int
        Index of the first snapshot to integrate from.
    propagation_on : bool
        GEARFeedback:ISRF_propagation of the run. With propagation off the
        module applies no relaxation at all, so the predicted drift is zero
        (`representative_c_hyp`).

    Returns
    -------
    dict
        ``measured`` (ln(u_LW/n_LW) drift from `run[start]`), ``dilated_lna``
        (the cumulative trapezoid), ``trapezoid_bound`` (its cumulative
        discretisation bound) and ``span`` (total ln a covered).
    """
    ref = run[start]
    ln_ratio0 = float(
        np.log(ledger_mean(ref, "u_LW", False) / ledger_mean(ref, "n_LW", False))
    )
    measured = np.full(len(run), np.nan)
    for i in range(start, len(run)):
        snap = run[i]
        ratio = ledger_mean(snap, "u_LW", False) / ledger_mean(snap, "n_LW", False)
        measured[i] = np.log(ratio) - ln_ratio0
    integral = trapezoid_c_hyp_integral(
        run, representative_c_hyp(run, propagation_on), start
    )
    return {
        "measured": measured,
        "dilated_lna": integral["dilated_lna"],
        "trapezoid_bound": integral["trapezoid_bound"],
        "span": integral["span"],
    }


def float_divergence_pull(
    snap: Dict, u_key: str, div_key: str, dt_step_cgs: float
) -> Optional[float]:
    """Return one step's largest relative transport pull on a moment's ``u``.

    The relaxation update is double, but the flux divergence it subtracts is
    still float end to end (radiation_end_force_propagation,
    radiation_isrf.c:913-917), so each step's increment carries a float32
    relative error. This returns the increment's own size relative to ``u``,
    ``max_i |dt * div_i / u_i|``, which the caller multiplies by
    ``FLOAT32_EPS`` and the run's step count. None when the snapshot carries
    no such field, or when the ratio is not everywhere finite, so a caller
    never folds a NaN into a bar.

    Parameters
    ----------
    snap : dict
        One snapshot, as `read_snapshot` returns it.
    u_key, div_key : str
        The moment's specific energy and its flux divergence.
    dt_step_cgs : float
        Longest step in proper time, seconds.

    Returns
    -------
    float or None
        The largest relative pull, or None.
    """
    div = snap.get(div_key)
    if div is None:
        return None
    pull = np.abs(dt_step_cgs * div / snap[u_key])
    if not np.all(np.isfinite(pull)):
        return None
    return float(np.max(pull))


def check_band_edge_ratio(
    run: List[Dict],
    band_edge_weights: Dict[str, float],
    dt_max: float,
    n_steps: float,
    dt_step_cgs: float,
    propagation_on: bool = True,
) -> bool:
    """Check (A3): d ln(u_LW/n_LW)/dt = -(lambda_E(LW) - lambda_N(LW)) H_dilated.

    The reference is independent of the module under test: lambda_E(LW) and
    lambda_N(LW) come from `radiation_set_band_edge_coefficients`'s own
    start-up log line, computed from the table before any step runs, a
    different code path from the per-step force update
    (`radiation_end_force_propagation`) that this check exercises. A bug
    that applied one shared lambda to both moments (the failure mode that
    function's own doxygen names) would leave the measured ratio flat while
    this prediction stays at its full, nonzero value, so the two cannot
    agree by construction the way a residual built from the same quantities
    the code computes would.

    Parameters
    ----------
    run : list of dict
        The run's snapshots, in time order.
    band_edge_weights : dict
        This run's own lambda_E(PE)/lambda_E(LW)/lambda_N(LW), from
        `read_band_edge_weights`.
    dt_max : float
        TimeIntegration:dt_max of the run, ln a units.
    n_steps : float
        Number of dt_max steps the run spans, from `step_count`.
    dt_step_cgs : float
        Longest step in proper time, seconds (the run's own dt_max
        converted, as `check_free_field`'s own `dt_step` already is).
    propagation_on : bool
        GEARFeedback:ISRF_propagation of the run. With propagation off
        `radiation_end_force_propagation` returns before its relaxation
        update, so the prediction is exactly zero and the measured ratio
        must stay flat: that is gated here, not skipped, the same way (A1)
        gates its own zero prediction on such a run.

    Returns
    -------
    bool
        Whether (A3) passed.
    """
    if propagation_on:
        start = band_edge_ratio_reference_index(run)
        if start is None or start >= len(run) - 1:
            print(
                "  (A3) SKIPPED: fewer than two snapshots have a genuine "
                "(non-clamped) HyperbolicPropagationSpeeds"
            )
            return True
    else:
        start = 0
        print(
            "  (A3): propagation off, so the predicted ln(u_LW/n_LW) drift is "
            "exactly 0 and the measured drift is gated against the "
            "error-only bar alone"
        )
    if any(s["n_LW"] is None for s in run[start:]):
        raise RuntimeError(
            "cosmological free_field run but LWPhotonSpecificEnergies is "
            "absent from a snapshot: (A3) needs the run's own photon-number "
            "moment and must not silently skip the one check that field "
            "supports."
        )
    checked = ("u_LW", "n_LW", "c_hyp") if propagation_on else ("u_LW", "n_LW")
    for i in range(start, len(run)):
        snap = run[i]
        for name in checked:
            if not np.all(np.isfinite(snap[name])):
                print(f"  FAIL: (A3): non-finite {name} in snapshot {i}")
                return False
        if np.any(snap["u_LW"] <= 0.0) or np.any(snap["n_LW"] <= 0.0):
            print(f"  FAIL: (A3): non-positive u_LW or n_LW in snapshot {i}")
            return False
        # This check's whole discriminating power depends on u_LW/n_LW
        # actually carrying double precision end to end, not merely being
        # cast to float64 on read (`physical()` always upcasts, so a value
        # narrowed to float32 anywhere upstream -- storage, or an I/O path
        # that reads it through a float buffer -- looks identical to a
        # genuine double once loaded): before the double widen, both
        # bands' exponentials rounded to the same float32 value and this
        # identity read exactly 0.
        for name in ("u_LW", "n_LW"):
            arr = snap[name]
            if np.array_equal(arr.astype(np.float32).astype(np.float64), arr):
                print(
                    f"  FAIL: (A3): {name} in snapshot {i} is exactly "
                    "float32-representable end to end: the double-precision "
                    "path this check depends on is not actually live, "
                    "whatever the HDF5 dtype claims"
                )
                return False

    errors = band_edge_ratio_errors(run, start, propagation_on)
    delta_lambda = band_edge_weights["LW"] - band_edge_weights["N_LW"]
    predicted = -delta_lambda * errors["dilated_lna"]
    residual = errors["measured"] - predicted

    quantisation = (
        band_edge_weight_log_precision(band_edge_weights["LW"])
        + band_edge_weight_log_precision(band_edge_weights["N_LW"])
    ) * float(np.max(np.abs(errors["dilated_lna"])))
    mean_c_hyp_ratio = (
        errors["dilated_lna"][-1] / errors["span"] if errors["span"] > 0 else 0.0
    )
    step_end_and_lag = (
        abs(delta_lambda) * mean_c_hyp_ratio * (0.75 * dt_max * errors["span"] + dt_max)
    )
    trapezoid = abs(delta_lambda) * float(errors["trapezoid_bound"][-1])

    # c_hyp and dt_prev are shared by both moments (one operator,
    # radiation_isrf.c:753-754) and cancel in the ratio to first order, so
    # they do not enter here. What survives is the FLOAT operands that
    # differ BETWEEN the two moments feeding the double relaxation update:
    # each moment's own flux divergence and dissipation source
    # (radiation_isrf.c:913-917 is still float end to end for the
    # divergence and dissipation inputs). Bounded by
    # eps_f * n_steps * the divergence's own relative pull on u this step,
    # read from the run's own snapshots, never fitted.
    float_terms = []
    for i in range(start, len(run)):
        snap = run[i]
        for u_key, div_key in (("u_LW", "div_LW"), ("n_LW", "div_LW_photon")):
            float_terms.append(float_divergence_pull(snap, u_key, div_key, dt_step_cgs))
    float_terms = [term for term in float_terms if term is not None]
    if float_terms:
        float_residual = FLOAT32_EPS * n_steps * max(float_terms)
        float_note = "flux-divergence-scaled"
    else:
        # No flux-divergence fields on this run (older snapshot): fall back
        # to the double relaxation chain's own operand count instead of a
        # value that would need those fields to derive -- about ten double
        # FLOPs (exp, expm1, phi, multiply-adds) per moment update.
        float_residual = 10.0 * np.finfo(np.float64).eps * n_steps
        float_note = "fallback, no flux-divergence fields"
    bar = quantisation + step_end_and_lag + trapezoid + float_residual

    print(
        f"  (A3) band-edge ratio: lambda_E(LW)-lambda_N(LW)={delta_lambda:.6g}, "
        f"int(c_hyp/c) d ln a={errors['dilated_lna'][-1]:.3e}, predicted "
        f"{predicted[-1]:.3e}, measured {errors['measured'][-1]:.3e}"
    )
    print(
        f"  bar (A3): quantisation {quantisation:.2e} + step-end/lag "
        f"{step_end_and_lag:.2e} + trapezoid {trapezoid:.2e} + float "
        f"residual ({float_note}) {float_residual:.2e} = {bar:.2e}"
    )
    worst = float(np.max(np.abs(residual[start:])))
    return gate("(A3) ln(u_LW/n_LW), worst |residual|", worst, bar)


def check_free_field(opt: argparse.Namespace) -> bool:
    """Check Eqs. (A1) and (A2)."""
    run = load_run(opt.snapshots)
    dt_max_parameter = read_dt_max(opt.snapshots, opt.dt_max)
    # Every term below built from a step size or a step count uses the
    # QUANTISED dt_max from the run's own start-up log, not the raw
    # parameter: see `read_timeline_dt_max`.
    dt_max = read_timeline_dt_max(opt.snapshots, log=opt.log)
    if dt_max is None:
        dt_max = dt_max_parameter
        print(
            "  NOTE: no engine_config line announced the quantised "
            "maximal step, so the raw TimeIntegration:dt_max parameter is "
            "used: the step count is then under-counted and the one-step "
            "lag term over-counted, in opposite directions"
        )
    cosmological = run[0]["cosmological"]
    propagation_on = read_isrf_propagation(opt.snapshots)
    # A run with propagation off applies no relaxation at all
    # (radiation_end_force_propagation returns first), so its predicted
    # drift is exactly zero, whatever its cosmology.
    cosmo_decay = cosmological and propagation_on
    # A c_hyp that is bit-uniform across the box makes the ledger's own
    # 1/c_hyp weight cancel between its numerator and its denominator, so the
    # metric becomes the mass-weighted box mean the transport conserves. That
    # removes the reweighting term a varying-c_hyp run carries, of order
    # spread(u)*spread(1/c_hyp), and leaves the float divergence floor as
    # (A1)'s whole error budget. Uniformity is the physical condition the
    # gate rests on, so it is read from the recorded field itself and not
    # from whichever parameter produced it: the debug pin and
    # ISRF_c_hyp_scheme 2's fixed fraction both give a bit-uniform speed.
    # The pin parameter is still read, because a pin the module did not apply
    # must FAIL rather than fall back to the ungated report.
    c_hyp_pin = read_c_hyp_pin(opt.snapshots)
    c_hyp_scheme = read_c_hyp_scheme(opt.snapshots)
    # Two configurations CLAIM one speed for the whole box: the debug pin, and
    # scheme 2, whose ISRF_c_hyp_fixed_fraction_of_c is a single fraction of c
    # (feedback_properties.h already errors at start-up if that key is not
    # positive under scheme 2, or positive without it, so the scheme number
    # alone is the claim). Either claim must FAIL when the recorded field does
    # not honour it, rather than fall back to the ungated report: the gate
    # below is the only bound on this fixture's transport ledger, and a claim
    # the module did not deliver would otherwise remove it silently.
    pin_claimed = c_hyp_pin is not None and c_hyp_pin > 0.0
    uniform_claimed = pin_claimed or c_hyp_scheme == 2
    c_hyp_spread = c_hyp_spatial_spread(run)
    uniform_c_hyp = c_hyp_spread is not None and c_hyp_spread == 0.0
    n_steps = step_count(run, dt_max)
    use_c_hyp = use_c_hyp_ledger(run, opt.snapshots, "run")
    ledger = "sum m u / c_hyp" if use_c_hyp else "sum m u"
    # (A1) predicts from the run's OWN recorded c_hyp, so it measures the
    # drift from the first snapshot whose c_hyp is a genuine per-step rate:
    # snapshot 0 of an unpinned run still holds the module's first-init
    # light-speed clamp (`band_edge_ratio_reference_index`), and neither
    # substituting a later snapshot's value for it nor integrating the clamp
    # itself is a measurement of anything the particles experienced.
    ref_index = 0
    if cosmo_decay:
        genuine = band_edge_ratio_reference_index(run)
        if genuine is None or genuine >= len(run) - 1:
            raise RuntimeError(
                "cosmological free_field run with propagation on, but fewer "
                "than two snapshots carry a genuine (non-clamped) "
                "HyperbolicPropagationSpeeds: (A1) predicts the dilated "
                "Hubble decay from the run's own per-step c_hyp and has no "
                "a-priori substitute for it (a bound from "
                "TimeIntegration:dt_min is permissive enough on every "
                "fixture this file runs that it never fails). Rerun with a "
                "build that writes HyperbolicPropagationSpeeds."
            )
        for i in range(genuine, len(run)):
            c_hyp = run[i]["c_hyp"]
            if not np.all(np.isfinite(c_hyp)) or np.any(c_hyp <= 0.0):
                raise RuntimeError(
                    f"HyperbolicPropagationSpeeds in snapshot {i} is not "
                    "everywhere finite and positive, after snapshot "
                    f"{genuine} carried a genuine sample: (A1)'s prediction "
                    "would be built on it."
                )
        ref_index = genuine
    # (A2)'s reference needs the LW field's own decay exponent, and it is exact
    # only when the recorded speed is one value for the whole run: beta is then
    # lambda_E(LW)*c_hyp/c with no quadrature. A varying speed would need the
    # trapezoid, which must not start at a snapshot whose recorded speed is
    # zero, so it is not attempted here and the omission is printed.
    h2_field_decay = 0.0
    h2_decay_note = "not applied: non-cosmological, no expansion to decay under"
    if cosmo_decay:
        weights_for_a2 = read_band_edge_weights(opt.snapshots, log=opt.log)
        if weights_for_a2 is None:
            h2_decay_note = (
                "not applied: the run's log announced no lambda_E(LW), so the "
                "exponent cannot be formed from the run's own table"
            )
        elif not uniform_c_hyp:
            h2_decay_note = (
                "not applied: the recorded c_hyp is not one value for the run, "
                "so beta is not a constant and (A2) carries the omission"
            )
        else:
            c_hyp_cgs = float(np.median(run[ref_index]["c_hyp"]))
            h2_field_decay = weights_for_a2["LW"] * c_hyp_cgs / C_LIGHT_CGS
            h2_decay_note = (
                f"lambda_E(LW) {weights_for_a2['LW']:.5g} x c_hyp/c "
                f"{c_hyp_cgs / C_LIGHT_CGS:.4e} = {h2_field_decay:.4e}"
            )
    errors = free_field_errors(run, use_c_hyp, ref_index, h2_field_decay)
    print(f"  (A2) LW field decay exponent beta: {h2_decay_note}")
    integral = trapezoid_c_hyp_integral(
        run, representative_c_hyp(run, cosmo_decay), ref_index
    )
    dilated_lna = integral["dilated_lna"]
    span = integral["span"]
    mean_c_hyp_ratio = dilated_lna[-1] / span if span > 0.0 else 0.0
    elapsed = np.array([s["time"] - run[0]["time"] for s in run])[1:]
    print(
        f"free_field: cosmological={cosmological}, propagation={propagation_on}, "
        f"a {run[0]['a']:.6g} -> {run[-1]['a']:.6g}, {len(run)} snapshots, "
        f"{n_steps:.0f} steps of the run's quantised dt_max {dt_max:.6g} "
        f"(parameter {dt_max_parameter:.6g}), H2 exponent "
        f"{errors['exponent']:.3f}, ledger {ledger}"
    )
    print(f"  (A1) reference snapshot: {ref_index}, ln a span from it {span:.6g}")
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
        reference_ref_index = band_edge_ratio_reference_index(reference_run)
        reference = free_field_errors(
            reference_run,
            use_c_hyp_reference,
            0 if reference_ref_index is None else reference_ref_index,
            h2_field_decay,
        )

    ok = True
    spread_note = (
        "absent or non-positive" if c_hyp_spread is None else f"{c_hyp_spread:.3e}"
    )
    pin_note = (
        "unrecorded"
        if c_hyp_pin is None
        else f"{c_hyp_pin:.9g} (internal velocity units)"
    )
    print(
        f"  (A1) worst spatial spread (max-min)/max of the recorded c_hyp "
        f"{spread_note}, so the uniform-weight ledger is "
        f"{'LIVE' if uniform_c_hyp else 'not live'}; "
        f"GEARFeedback:ISRF_c_hyp_pin_for_debugging = {pin_note}, "
        f"ISRF_c_hyp_scheme = {c_hyp_scheme}"
    )
    if uniform_claimed and not uniform_c_hyp:
        source = (
            "a c_hyp pin is set"
            if pin_claimed
            else "ISRF_c_hyp_scheme is 2, which is one fixed speed for the " "whole box"
        )
        # Two different failures, and naming the wrong one sends the
        # investigation the wrong way: no snapshot carrying a usable speed is
        # usually propagation switched off, while a spread is the claim not
        # being honoured.
        cause = (
            "no snapshot records a usable speed (propagation off, the field "
            "absent, or a non-positive value after the first step)"
            if c_hyp_spread is None
            else "the recorded HyperbolicPropagationSpeeds is not bit-uniform "
            "across the box"
        )
        print(
            f"  FAIL: (A1) {source} but {cause}, so the uniform-weight ledger "
            "this leg is gated on is not live"
        )
        ok = False
    # Longest step in proper time, from the run's quantised dt_max (ln a
    # with cosmology).
    dt_step = (
        dt_max / hubble_rate_cgs(run[0]["a"], run[0])
        if cosmological
        else dt_max * run[0]["time_unit"]
    )
    if cosmo_decay:
        print(
            f"  int (c_hyp/c) d ln a from the run's own "
            f"HyperbolicPropagationSpeeds: {dilated_lna[-1]:.4e} "
            f"(mean c_hyp/c {mean_c_hyp_ratio:.4e}, trapezoid "
            f"discretisation bound {integral['trapezoid_bound'][-1]:.3e})"
        )
    # Per-band redshift depth is lambda_E(band)*H_dilated, not H_dilated
    # alone (radiation_end_force_propagation): read the run's own
    # table-derived lambda_E(PE)/lambda_E(LW) from its log rather than
    # hardcoding a number here, so this bar tracks the radiation table
    # instead of one measurement of it.
    band_edge_weights = None
    if cosmological:
        band_edge_weights = read_band_edge_weights(opt.snapshots, log=opt.log)
        if band_edge_weights is None:
            raise RuntimeError(
                "cosmological free_field run but no output.log announced "
                "radiation_set_band_edge_coefficients()'s lambda_E(PE)/"
                "lambda_E(LW): (A1) predicts each band's own drift from "
                "them and cannot default to lambda=1, which is neither "
                "band's coefficient and would make the residual meaningless "
                "rather than merely loose. Pass --log or run beside the log "
                "this fixture wrote."
            )
        print(
            f"  band-edge weights from the log: lambda_E(PE)="
            f"{band_edge_weights['PE']:.5g}, lambda_E(LW)="
            f"{band_edge_weights['LW']:.5g}"
        )
    # The LW-to-PE band-edge transfer (radiation_end_force_propagation) adds
    # (lambda_E(LW) - 1)*H_dilated*u_LW into PE's own `u` every step, on top
    # of PE's own -lambda_E(PE)*H_dilated*u_PE decay, so PE's own signed
    # drift coefficient is [lambda_E(LW) - 1]*r - lambda_E(PE), with r the
    # box-mean u_LW/u_PE at the reference snapshot (not assumed to be 1,
    # though this fixture's own IC sets u_pe = u_lw). At this fixture's own
    # weights that coefficient is POSITIVE: PE gains more from the transfer
    # than it loses to its own decay, so its box mean rises. LW receives no
    # such term and simply decays at -lambda_E(LW).
    r_reference = 0.0
    r_drift = 0.0
    if cosmo_decay:
        r_series = errors["r_lw_over_pe"]
        r_reference = float(r_series[ref_index])
        r_drift = float(np.max(np.abs(r_series[ref_index:] - r_reference)))
        print(
            f"  box-mean u_LW/u_PE (r, the LW-to-PE transfer's own scale "
            f"factor): {r_reference:.6g} at the reference snapshot, drifting "
            f"by at most {r_drift:.2e} over the run"
        )
    for band in ["PE", "LW"]:
        # SIGNED: the residual below is two-sided, so a prediction of the
        # wrong sign must not compare equal to the right one.
        if cosmo_decay:
            drift_coeff = (
                (band_edge_weights["LW"] - 1.0) * r_reference - band_edge_weights["PE"]
                if band == "PE"
                else -band_edge_weights["LW"]
            )
        else:
            drift_coeff = 0.0
        predicted = drift_coeff * dilated_lna
        residual = errors[band] - predicted
        if not np.all(np.isfinite(residual[ref_index:])):
            print(f"  FAIL: (A1) u_{band}: non-finite measured drift or residual")
            ok = False
            continue
        # Error-only bar. Nothing here is the predicted decay itself: that
        # is on the other side of the residual now.
        #
        # The log line's own %.5g precision on each lambda, propagated
        # through this band's coefficient as a half-ulp bound on each, times
        # the run's own integral.
        quantisation = 0.0
        if cosmo_decay:
            quantisation = band_edge_weight_log_precision(band_edge_weights["LW"]) * (
                r_reference if band == "PE" else 1.0
            )
            if band == "PE":
                quantisation += band_edge_weight_log_precision(band_edge_weights["PE"])
            quantisation *= abs(dilated_lna[-1])
        # The trapezoid's own discretisation bound, and the step-end freeze
        # of H ((3/4) dlna_step per unit ln a, matter domination) plus the
        # one-step snapshot lag, both scaled by this band's coefficient.
        trapezoid = abs(drift_coeff) * float(integral["trapezoid_bound"][-1])
        step_end_and_lag = (
            abs(drift_coeff) * mean_c_hyp_ratio * (0.75 * dt_max * span + dt_max)
        )
        # PE's coefficient holds r at its reference value; r itself drifts as
        # the two bands separate, which is second order in the cosmological
        # perturbation and bounded by that drift.
        transfer_linearisation = (
            abs(band_edge_weights["LW"] - 1.0) * r_drift * abs(dilated_lna[-1])
            if cosmo_decay and band == "PE"
            else 0.0
        )
        # The a-independent floor: the float32 flux divergence the double
        # relaxation update subtracts every step (`float_divergence_pull`),
        # or, when a non-cosmological run of the same configuration was
        # paired in, twice its own measured drift, which measures that floor
        # instead of bounding it.
        pull = [
            value
            for value in (
                float_divergence_pull(snap, f"u_{band}", f"div_{band}", dt_step)
                for snap in run[ref_index:]
            )
            if value is not None
        ]
        if pull:
            float_residual = FLOAT32_EPS * n_steps * max(pull)
            float_note = "flux-divergence-scaled"
        else:
            # No flux-divergence field on this run (older snapshot): fall
            # back to the double relaxation chain's own operand count, about
            # ten double FLOPs per moment update, as (A3) does.
            float_residual = 10.0 * np.finfo(np.float64).eps * n_steps
            float_note = "fallback, no flux-divergence field"
        measured_nc = (
            0.0
            if reference is None
            else float(np.max(np.abs(reference[band][reference["ref_index"] :])))
        )
        if not np.isfinite(measured_nc):
            raise RuntimeError(
                f"the reference run's u_{band} drift is not finite -- "
                "refusing to build a bar from a corrupted reference."
            )
            # The reference run's own measured error is NOT folded in. A bar
            # taken as a multiple of a measurement made on the runs it then
            # judges is fitted to the data, which this project's protocol
            # forbids and the operator has ruled against in terms: an honest
            # red is worth more than a green bought that way. The reference
            # run's error is still PRINTED, as a diagnostic.
        bar = (
            quantisation
            + trapezoid
            + step_end_and_lag
            + transfer_linearisation
            + float_residual
        )
        worst = float(np.max(np.abs(residual[ref_index:])))
        print(
            f"  (A1) u_{band}: drift coefficient {drift_coeff:+.6g}, predicted "
            f"{predicted[-1]:+.4e}, measured {errors[band][-1]:+.4e}"
        )
        if not cosmological and uniform_c_hyp:
            # GATED, on the float floor alone. With c_hyp bit-uniform the
            # ledger's 1/c_hyp weight cancels identically, so the metric is
            # the mass-weighted box mean the pairwise transport conserves,
            # and the only error left is the float32 flux divergence the
            # double relaxation update subtracts every step. That floor is
            # mechanism-derived (`float_divergence_pull`), the same
            # construction (A3) uses, and it is used here on its own rather
            # than through `bar`, so a `--reference` handed to a
            # non-cosmological run cannot inflate it.
            # Both schemes' operators conserve sum m u / c_hyp. Scheme 2's
            # speed is one value for the whole box, so there that is
            # proportional to sum m u; scheme 4's weight cancels only because
            # it is uniform here, and the two coincide exactly when it is.
            if c_hyp_scheme == 2:
                why = "this scheme's one speed makes sum m u / c_hyp proportional to sum m u"
            elif use_c_hyp:
                why = "the ledger weight cancels between the two sums"
            else:
                why = "a uniform weight makes the two ledgers the same functional"
            print(
                f"  (A1) u_{band}: c_hyp bit-uniform, so {why} and the metric "
                f"is the mass-weighted box mean; bar is the analytic float "
                f"floor ({float_note}) {float_residual:.2e} on its own"
            )
            ok &= gate(
                f"(A1) ln u_{band} drift, worst |residual|, {ledger}, "
                f"c_hyp bit-uniform",
                worst,
                float_residual,
            )
            continue
        if not cosmological:
            # Reported, not gated: see this module's docstring. Both numbers
            # are printed because they are not the same one and both are
            # used downstream -- the second is what `--reference` doubles
            # into a cosmological run's own bar.
            genuine = band_edge_ratio_reference_index(run)
            if genuine is None:
                from_genuine = worst
            else:
                # Rebaselined on the genuine snapshot, not merely sliced
                # from it: that is the drift `--reference` reads, and it
                # differs from a slice of this one, which is still measured
                # from snapshot 0.
                from_genuine = float(
                    np.max(
                        np.abs(
                            free_field_errors(run, use_c_hyp, genuine, h2_field_decay)[
                                band
                            ][genuine:]
                        )
                    )
                )
            print(
                f"  (A1) u_{band} SKIPPED: non-cosmological, so the "
                f"predicted drift is identically zero and every "
                f"cosmological bar term vanishes with it. Measured "
                f"|drift|: {worst:.3e} worst over the whole run, "
                f"{from_genuine:.3e} worst from the first genuine c_hyp "
                f"snapshot (the value --reference doubles into a "
                f"cosmological run's bar); analytic float floor "
                f"({float_note}) {float_residual:.2e}"
            )
            continue
        print(
            f"  bar u_{band}: quantisation {quantisation:.2e} + trapezoid "
            f"{trapezoid:.2e} + step-end/lag {step_end_and_lag:.2e} + "
            f"transfer linearisation {transfer_linearisation:.2e} + "
            f"float residual ({float_note}) {float_residual:.2e}"
            f" = {bar:.2e}; the reference run's own error {measured_nc:.2e} is "
            f"REPORTED and deliberately NOT in the bar"
        )
        if cosmo_decay:
            # Resolution self-test, not a bar term: nothing bounds the
            # reference drift that `2 x non-cosmological` carries into the
            # bar, so a reference run with a large drift can raise the bar
            # above the very decay this band predicts. A residual built on
            # a module that applied no decay at all is max|predicted| to
            # within the run's own floor, so once the bar reaches that
            # value the gate admits the total absence of the effect and has
            # stopped being a gate. Necessary, not sufficient: inside one
            # floor of that value the gate is still weak. Derived from the
            # no-decay failure mode, not sized from any measurement.
            signal = float(np.max(np.abs(predicted[ref_index:])))
            if not np.isfinite(signal) or not np.isfinite(bar) or bar >= signal:
                print(
                    f"  FAIL: (A1) u_{band}: bar {bar:.3e} is not below the "
                    f"predicted decay it has to discriminate {signal:.3e}, "
                    "so a run that applied no decay at all would pass here"
                )
                ok = False
                continue
        ok &= gate(f"(A1) ln u_{band} drift, worst |residual|, {ledger}", worst, bar)

    # H2, per snapshot, every term RELATIVE to the predicted exponent, which
    # is what errors["H2"] is: implicit solve (k dt/2), one-step snapshot lag
    # (dt/t without cosmology; with it the rate varies as a^-3, same order),
    # float32 storage of the species fraction; with cosmology the step-end
    # rate adds 1.5 dlna_step ((p/2) dlna_step at p = 3, this module's own
    # Bars section, not 2.0 at p = 4: A2's rate no longer carries u_LW's own
    # a0/a factor). H2I_frac is a float rewritten every step
    # (cooling_struct.h), so ln x_H2 carries at most one unit round-off
    # (eps/2) per step taken so far, an ABSOLUTE ln error divided by the
    # predicted exponent it is compared with.
    exponent_so_far = errors["rate"] * errors["integral"][1:]
    steps_so_far = np.array(
        [step_count(run[: i + 1], dt_max) for i in range(1, len(run))]
    )
    if not (
        np.all(np.isfinite(errors["H2"]))
        and np.all(np.isfinite(exponent_so_far))
        and np.all(exponent_so_far > 0.0)
    ):
        print(
            "  FAIL: (A2) the H2 error or the predicted exponent is "
            "non-finite, or the predicted exponent is not positive"
        )
        return False
    float_storage = 0.5 * FLOAT32_EPS * steps_so_far / exponent_so_far
    budget = 0.5 * errors["rate"] * dt_step + dt_step / elapsed + float_storage
    cosmo = 1.5 * dt_max if cosmological else 0.0
    measured_nc = np.zeros_like(budget)
    if reference is not None:
        # By ELAPSED TIME, never by snapshot index: cosmo_timeline.py spaces a
        # non-cosmological run uniformly in t and a cosmological one uniformly
        # in ln a, so equal snapshot counts do not mean equal elapsed times.
        # The lag term is largest at the earliest times, so pairing by index
        # would mix an early reference error into a late cosmological bar.
        ref_times = reference["times"][1:]
        ref_error = np.abs(reference["H2"])
        if not np.all(np.isfinite(ref_error)) or not np.all(np.isfinite(ref_times)):
            raise RuntimeError(
                "the reference run's H2 errors or times are not all finite -- "
                "refusing to build a bar from a corrupted reference."
            )
        # A run spaced uniformly in ln a always starts inside the reference's
        # first proper-time interval, so np.interp clamps the earliest
        # snapshots to the reference's first error. That is the strict
        # direction for a lag-dominated reference, whose true error grows
        # towards t = 0, and the run's own dt/elapsed term in `budget`
        # already covers its own early-time lag. It is reported, not fatal.
        below = int(np.sum(elapsed < ref_times[0]))
        if below:
            print(
                f"  {below} snapshot(s) start before the reference's first "
                f"output ({ref_times[0]:.4e} s); their reference term is "
                f"clamped to {ref_error[0]:.3e}, which understates it."
            )
        # Beyond the reference's last output there is no such argument, and
        # the clamp would understate an error that is still evolving.
        if elapsed[-1] > ref_times[-1] * (1.0 + 1e-6):
            raise RuntimeError(
                f"the reference run ends at {ref_times[-1]:.4e} s but this "
                f"run reaches {elapsed[-1]:.4e} s; np.interp would silently "
                "clamp there. Extend the reference run's time_end."
            )
        measured_nc = np.interp(elapsed, ref_times, ref_error)
    # See (A1) above: the reference run's measured error is reported, not
    # folded into the bar.
    bar = budget + cosmo

    # Resolution self-test, matching the one (A1) already carries: a bar that
    # has grown to or above the signal it must discriminate gates nothing, and a
    # run that applied no H2 photodissociation at all would pass it. This can only turn a
    # VACUOUS pass into a failure. It is needed because the bar above takes a
    # `max` against twice the reference run's own measured error, which the
    # reference run can make arbitrarily large.
    signal_a2 = np.abs(errors["rate"] * errors["integral"][1:])
    if np.any(bar >= signal_a2):
        worst_k = int(np.argmax(bar / np.where(signal_a2 > 0.0, signal_a2, np.inf)))
        print(
            f"  FAIL: (A2) bar {np.atleast_1d(bar)[worst_k]:.3e} is not below the "
            f"predicted exponent it has to discriminate "
            f"{signal_a2[worst_k]:.3e}, so a run that applied no H2 "
            "photodissociation at all would pass here"
        )
        return False
    ratio = np.abs(errors["H2"]) / bar
    k = int(np.argmax(ratio))
    print(
        f"  bar ln x_H2 (per snapshot): implicit solve {0.5 * errors['rate'] * dt_step:.1e} "
        f"+ lag dt/t {dt_step / elapsed[-1]:.1e} (end) to {dt_step / elapsed[0]:.1e} (first) "
        f"+ float32 H2I storage {float_storage[-1]:.1e} (end) to {float_storage[0]:.1e} (first), "
        f"2 x non-cosmological up to {2.0 * np.max(measured_nc):.1e}, step-end rate {cosmo:.1e}"
    )
    print(
        f"  worst snapshot {k + 1}: error {errors['H2'][k]:.3e}, bar {bar[k]:.3e}; "
        f"final error {errors['H2'][-1]:.3e}, bar {bar[-1]:.3e}"
    )
    ok &= gate("box-mean ln x_H2 exponent (A2), worst error/bar", float(ratio[k]), 1.0)

    if cosmological:
        ok &= check_band_edge_ratio(
            run, band_edge_weights, dt_max, n_steps, dt_step, propagation_on
        )
    else:
        print("  (A3) SKIPPED: non-cosmological, ratio preserved, lambda not tested")
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
    # Every step size and step count below uses the QUANTISED dt_max from
    # the run's own start-up log, not the raw parameter: see
    # `read_timeline_dt_max`.
    dt_max = read_timeline_dt_max(opt.snapshots, log=opt.log)
    if dt_max is None:
        dt_max = read_dt_max(opt.snapshots, opt.dt_max)
        print(
            "  NOTE: no engine_config line announced the quantised "
            "maximal step, so the raw TimeIntegration:dt_max parameter is "
            "used: the step count is then under-counted and the one-step "
            "lag term over-counted, in opposite directions"
        )
    cosmological = run[0]["cosmological"]
    n_steps = step_count(run, dt_max)
    errors = dust_absorption_errors(run, opt.c_hyp_pin * 1e5)
    dt_step = (
        dt_max / hubble_rate_cgs(run[0]["a"], run[0])
        if cosmological
        else dt_max * run[0]["time_unit"]
    )
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
        # The flux divergence the double relaxation update subtracts is
        # still float: see `float_divergence_pull`. u itself is double, so
        # it contributes no round-off term.
        pulls = [
            float_divergence_pull(snap, f"u_{band}", f"div_{band}", dt_step)
            for snap in run[1:]
        ]
        if any(pull is None for pull in pulls):
            print(f"  NOTE: no finite {band} flux divergence field, float term 0")
            float_floor = 0.0
        else:
            float_floor = FLOAT32_EPS * n_steps * max(pulls)
        budget = float_floor + per_step + depth * errors["density_drift"]
        # The kappa step-end term is unaffected by the Hubble-term dilation
        # (kappa was already correctly dilated); the H-alone step-end term is
        # dilated by c_pin/c, same as the B1 formula's own second term above.
        c_pin_ratio = (opt.c_hyp_pin * 1e5) / C_LIGHT_CGS if cosmological else 0.0
        cosmo = (
            1.5 * dt_max * depth + c_pin_ratio * 0.75 * dt_max * span
            if cosmological
            else 0.0
        )
        # See (A1): reported, not folded in.
        bar = budget + cosmo

        # Resolution self-test, matching the one (A1) already carries: a bar that
        # has grown to or above the signal it must discriminate gates nothing, and a
        # run that applied no absorption at all would pass it. This can only turn a
        # VACUOUS pass into a failure. It is needed because the bar above takes a
        # `max` against twice the reference run's own measured error, which the
        # reference run can make arbitrarily large.
        signal_b1 = abs(float(errors[f"{band}_depth"]))
        if not np.isfinite(signal_b1) or bar >= signal_b1:
            print(
                f"  FAIL: (B1) {band}: bar {bar:.3e} is not below the predicted "
                f"absorption it has to discriminate {signal_b1:.3e}, so a run "
                "that applied no absorption at all would pass here"
            )
            return False
        print(
            f"  bar {band}: max(float32 flux divergence {float_floor:.1e} + "
            f"one-step lag {per_step:.1e} + density drift "
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
    # The lag term uses the QUANTISED dt_max from the run's own start-up
    # log, not the raw parameter: see `read_timeline_dt_max`.
    dt_max = read_timeline_dt_max(opt.snapshots, log=opt.log)
    if dt_max is None:
        dt_max = read_dt_max(opt.snapshots, opt.dt_max)
        print(
            "  NOTE: no engine_config line announced the quantised "
            "maximal step, so the raw TimeIntegration:dt_max parameter is "
            "used: the one-step lag term is then over-counted"
        )
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
    # - float32 storage of the gas u (InternalEnergies is a float in struct
    #   part and in the snapshot), relative to the difference.
    t = err["times"][1:]
    if not (
        np.all(np.isfinite(err["relative"])) and np.all(err["predicted"][1:] > 0.0)
    ):
        print(
            "  FAIL: (D2) the measured difference or the predicted heating is "
            "non-finite, or the predicted heating is not positive"
        )
        return False
    temperature_ratio = err["on_u"][1:] / err["dark_u"][1:]
    # Neutral atomic gas, mu = 4/(1 + 3 X).
    t_dark = (
        (2.0 / 3.0)
        * (4.0 / (1.0 + 3.0 * HYDROGEN_MASS_FRACTION))
        * M_H_CGS
        * err["dark_u"][1:]
        / 1.380649e-16
    )
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        boost = np.expm1(92.0 / t_dark * (1.0 - 1.0 / temperature_ratio))
    cooling = (
        np.abs(err["dark_u"][1:] - err["dark_u"][0]) * boost / err["predicted"][1:]
    )
    lag = dt_step / t
    storage = 4.0 * FLOAT32_EPS * err["on_u"][1:] / err["predicted"][1:]
    bar = lag + cooling + storage
    # A non-finite bar term makes every ratio below 0 or NaN, and 0 passes.
    if not np.all(np.isfinite(bar)):
        print(
            "  FAIL: (D2) a bar term is not finite (cold dark gas overflows "
            "the C+ cooling boost): the bar is meaningless"
        )
        return False
    if opt.reference:
        if opt.reference_dark is None:
            raise RuntimeError("--reference needs --reference-dark for photoelectric")
        ref = photoelectric_errors(
            load_run(opt.reference), load_run(opt.reference_dark)
        )
        # The reference enters the bar, so a non-finite reference error
        # would raise the bar to infinity and pass every residual.
        if not (
            np.all(np.isfinite(ref["relative"]))
            and np.all(np.isfinite(ref["times"]))
            and np.all(ref["predicted"][1:] > 0.0)
        ):
            print(
                "  FAIL: (D2) the reference run's difference or predicted "
                "heating is non-finite, or the predicted heating is not "
                "positive: refusing to build a bar from it"
            )
            return False
        # Reference error at the same elapsed times.
        ref_error = np.interp(t, ref["times"][1:], np.abs(ref["relative"]))
        # See (A1): reported, not folded in. `ref_error` is printed below.

        # Resolution self-test, matching the one (A1) already carries: a bar that
        # has grown to or above the signal it must discriminate gates nothing, and a
        # run that applied no photoelectric heating at all would pass it. This can only turn a
        # VACUOUS pass into a failure. It is needed because the bar above takes a
        # `max` against twice the reference run's own measured error, which the
        # reference run can make arbitrarily large.
        signal_d2 = np.abs(err["predicted"][1:])
        if np.any(bar >= signal_d2):
            worst_k = int(np.argmax(bar / np.where(signal_d2 > 0.0, signal_d2, np.inf)))
            print(
                f"  FAIL: (D2) bar {np.atleast_1d(bar)[worst_k]:.3e} is not below "
                f"the predicted heating it has to discriminate "
                f"{signal_d2[worst_k]:.3e}, so a run that applied no "
                "photoelectric heating at all would pass here"
            )
            return False
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


def log_step_quantisation(token: str) -> float:
    """Return the relative quantisation of a step size read from the log text.

    Parameters
    ----------
    token
        The step-size field exactly as the log printed it.

    Returns
    -------
    float
        Half a unit in the last printed decimal of the mantissa, relative,
        dimensionless.

    Raises
    ------
    RuntimeError
        When the token is not in ``%14e`` form with a mantissa of 1 or more.
        Returning the format's worst case instead would LOOSEN the bar, by up
        to 10x on this term, in exactly the situation where the reference is
        not understood. The parse branch cannot fire, because the caller has
        already read the same field with ``float()``; the mantissa branch CAN,
        because step 0 prints ``0.000000e+00`` in that field
        (``src/engine.c``'s step line, from ``e->time_step`` before the first
        step). The caller rejects a zero step before ever calling this, so the
        raise is a backstop and not a live path.
    """
    try:
        mantissa = abs(float(token.split("e")[0]))
    except (ValueError, IndexError):
        mantissa = float("nan")
    if not np.isfinite(mantissa) or mantissa < 1.0:
        raise RuntimeError(
            f"the step size was printed as {token!r}, which is not the "
            "%14e form this term's quantisation is derived from: refusing "
            "to substitute the format's worst case, which would loosen the "
            "bar"
        )
    return LOG_STEP_MANTISSA_HALF_ULP / mantissa


def dust_band_ratio_bar(signal: float) -> float:
    """Return the float-arithmetic bar on the dusty band-ratio residual.

    THE CANONICAL DERIVATION. `ISRFInjectionConservation` carries the same
    metric and must keep the same term list.

    The residual ``R_j = ln(u_PE/u_LW) - ln(L_PE/L_LW) - (1 - rho) tau_LW``,
    with ``rho = sigma_PE/sigma_LW = 0.6`` exactly
    (``radiation.h``'s RADIATION_SIGMA_D_PE_CGS and _LW_CGS), is pinned to zero
    by an identity with no free parameter, so the only admissible discrepancy is
    float32 rounding. Writing each band's computed depth as
    ``tau_b (1 + e_s + e_b)``, with ``e_s`` the error of the prefactor the two
    bands SHARE and ``e_b`` the error of the per-band tail,

        R = e_s (tau_LW - tau_PE) - tau_PE e_PE + tau_LW e_LW
            + (the two expf errors) + (the two snapshot-store errors)

    so the terms split by how they enter, not merely by how many they are.

    CONSTANT, 2.04 u32, and it is the two ``expf`` extinction factors ALONE
    (``radiation_get_dust_extinction_factor``). An ``expf`` RELATIVE error is an
    ABSOLUTE error in the logarithm, which is why they do not scale with the
    signal. MEASURED, not assumed from the library's documented bound: this
    build's ``expf`` is accurate to at most 0.51 float32 ulps over the tau range
    these fixtures occupy, under the build's own flags, and one ulp is 2 u32, so
    each call contributes 1.02 u32.

    **CORRECTED 2026-10-01, and the first version of this bar was 4 u32.** It
    attributed two of those four to "the two band specific energies as stored in
    the snapshot". There is no such term: `PESpecificEnergies` and
    `LWSpecificEnergies` are written DOUBLE (`feedback_io.h`), the struct field is
    double, and THIS MODULE'S OWN docstring already says so in terms, that "a term
    modelling float32 round-off OF it is not a legitimate noise floor and none is
    used". So the constant was about twice what its stated mechanism supports. The
    direction was conservative, so it cannot have caused a false failure, but it
    mattered most at low signal where the constant is effectively the whole bar.

    SIGNAL-SCALED. With ``S`` the printed signal, ``tau_LW = S/(1 - rho) =
    2.5 S`` and ``tau_PE = rho tau_LW = 1.5 S``. So a shared-prefactor rounding
    enters with weight ``tau_LW - tau_PE = S``, while a per-band rounding enters
    with weight ``1.5 S`` or ``2.5 S``. On the shipped
    ``ISRF_extinction_path: pair_separation`` path the column is the
    star-to-particle separation ``r`` and the smoothing length never enters
    (``radiation_get_comoving_extinction_path``). The per-band tail carries
    three roundings in each band, at ``radiation_isrf.c``'s
    ``sigma_d_band_cgs * D_relative``, its divide by the shared denominator, and
    ``-kappa_eff * Sigma_gas_p``, so the per-band contribution alone is
    ``3 u32 (1.5 + 2.5) S = 12 u32 S``. The shared roundings are the column's
    ``path * rho_gas``, its comoving factor (exact at ``a = 1``), the
    metallicity normalisation, and the separation chain, whose weight is NOT one
    half-ulp but ``u32 |x - c->loc| / r``, because the injection differences the
    cell-relative offsets in float while this check rebuilds the separation in
    double.

    THE SCALED COEFFICIENT, 16.5, and which part of it is counted. Fifteen are
    counted: ``3 u32 (1.5 + 2.5) = 12`` per-band, and three shared at weight one
    each, being the column ``extinction_path * rho``, the metallicity
    normalisation ``max(Z,0)/RADIATION_GRACKLE_SOLAR_METAL_FRACTION``, and the
    comoving-to-physical factor, which is exact at ``a = 1`` and is not at the
    high-redshift leg. Two candidates are exactly zero rather than small: the
    dust-to-gas self-ratio is ``1.0f`` at the default, so that multiply is exact,
    and the two cross-section literals have bit-identical relative error on their
    float32 cast, so it acts as a common prefactor and cancels in the band
    difference. The remaining 1.5 is an ALLOWANCE and not a count: the separation
    chain's error is absolute at the cell-offset scale, so its weight is
    ``|x - c->loc| / r`` rather than one, and that ratio is a property of the
    cell layout the snapshot does not record. It is the one soft term here.

    Parameters
    ----------
    signal : float
        ``|(1 - sigma_PE/sigma_LW) * max_j tau_LW,j|``, dimensionless.

    Returns
    -------
    float
        The bar, dimensionless.
    """
    u32 = FLOAT32_EPS / 2.0
    # 0.51 float32 ulps per expf, measured on this build's own flags; one ulp
    # is 2 u32, so 1.02 u32 per call and two calls in the residual.
    expf_u32 = 2.0 * 1.02
    return expf_u32 * u32 + 16.5 * u32 * float(signal)


def injection_identity_bar(delta_t_token: str, n_lit: int) -> Dict:
    """Return the float-arithmetic bar on the injection identity's residual.

    At zero metallicity the injected energy is
    ``sum_j m_j u_j = Delta_t L sum_j weight_j``, with
    ``weight_j = m_j W(r_j, h_i) / enrichment_weight``
    (radiation_iact.h:227,297), so the residual is bounded by the roundings
    that stand between the two sides of that identity, and by nothing else.
    These candidate terms are zero by construction, not by measurement:

    - only one step accumulates, because the field is reset on the step's
      first touch by any star (radiation_iact.h:320-326);
    - no quadrature of the source rate enters, because the star's step is
      cached once per step (radiation_iact.h:124) and its band luminosity
      once per stellar-evolution update
      (src/feedback/GEAR/stellar_evolution.c:1394 and :1668), and both are
      held constant across its neighbours. This needs the snapshot's
      ``L_band`` to be the value the injection used, which holds while the
      luminosity is constant over the run;
    - the extinction factor is exactly ``1.0f`` at ``Z = 0``, because
      ``kappa_eff`` carries ``Z`` as a factor and ``expf(-0.f)`` is exact
      (radiation_isrf.c:1571 and radiation_get_dust_extinction_factor);
    - ``u`` carries no float32 term, being a double in ``struct part`` and in
      the snapshot (src/feedback/GEAR_thermal/feedback_struct.h:165);
    - this check's own float64 summation of ``n_lit`` terms costs
      ``(n_lit - 1) * 2**-53``, eleven orders below the terms kept below.

    Two premises the bar does not cover, because neither is a rounding:

    - the two loops must evaluate the same ``r``. The injection floors it at
      ``1e-3 h_i`` (radiation_iact.h:214) and the density loop does not
      (GEAR_thermal/feedback_iact.h:58), so a pair inside that radius has
      the two sides reading different kernel arguments. The effect is small
      rather than absent, but it is a statement about THIS kernel and not
      about every SWIFT kernel. Wendland C2 (``{4,-15,20,-10,0,1}``,
      ``kernel_gamma = 1.936492``) has a vanishing gradient at zero
      separation and a leading ``-10 x**2``, so at the floor the relative
      kernel change is ``10 * (1e-3 / 1.936492)**2 = 2.667e-06``. That is
      then weighted by the coincident particle's own share of the normalised
      weight sum, 0.196 to 0.225 over both legs' measured geometry, so a
      particle ADDED at the floor costs 5.2e-07 to 6.0e-07 across that range.
      A particle that REPLACES the closest neighbour costs more, because the
      share it takes is also the share the neighbour gave up: with
      ``s_rep = 1/(1/s_add - W_c/W(0))`` and ``W_c/W(0) = 0.651`` at the
      fixture's own closest pair, the same share range gives 6.0e-07 to
      7.0e-07. The figure carried here is 6.5e-07, which sits inside that
      range rather than bounding it. Both framings fit inside the margin the
      terms below leave, but neither is one of them, and the cost is not
      uniform in the neighbour count: as ``n_lit -> 1`` the share tends to 1 and the cost
      rises toward the full 2.7e-06. On both ISRFCosmology injection
      fixtures the floor is never approached, the closest pair sitting at
      ``r / h_star = 0.4675``;
    - both call sites must compile ``W`` to the same operations.
      ``kernel_eval`` and ``kernel_deval`` build it from the same
      coefficients by the same Horner recurrence, and under ``-ffast-math``
      that is a premise rather than an identity: in ``kernel_eval`` only the
      final ``w`` is live, so associative reassociation may reorder the
      chain, while in ``kernel_deval`` every intermediate ``w`` feeds
      ``dw_dx`` and pins it. No allowance is kept for it, because an
      allowance of this size would not bound it: the monomial Horner form is
      badly conditioned near the support edge, so a real divergence could
      reach about 1.7e-04, some 50 times the bar. It can therefore only
      produce a loud failure and never a silent pass, and it is measured to
      hold in the canonical build.

    The weight-sum closure below keeps the worst-case ``gamma_{n-1}`` growth
    and not a probabilistic ``sqrt(n) u32`` form. The square-root form is
    about 7x tighter, but it is a statistic whose spread is the size of its
    own bar, so it can fail a healthy run; a bound that holds for every
    summation order cannot.

    Parameters
    ----------
    delta_t_token
        The star's step size exactly as the run's log printed it.
    n_lit
        Number of illuminated gas particles. This is the neighbour count the
        star's own weight normalisation summed over.

    Returns
    -------
    dict
        Key "bar" and one key per term, all dimensionless.
    """
    u32 = FLOAT32_EPS / 2.0
    # The star's step is narrowed to a float before any neighbour reads it
    # (radiation_iact.h:124), while the reference side uses the double the
    # log printed.
    dt_float32 = u32
    dt_text = log_step_quantisation(delta_t_token)
    # enrichment_weight accumulates the m_j W_j products in float32 over the
    # star's neighbours (GEAR_thermal/feedback_iact.h:71) while the injection
    # re-sums the same products in double. For any summation order a float32
    # running sum of n positive terms carries at most gamma_{n-1}.
    n_round = max(int(n_lit) - 1, 0)
    growth = 1.0 - n_round * u32
    closure = np.inf if growth <= 0.0 else n_round * u32 / growth
    # Roundings that differ between the two loops. Every weight is positive,
    # so the normalised sum is a convex combination of the per-term relative
    # perturbations and each of these counts once, not n_lit times: the
    # injection's per-neighbour hi_inv_dim scaling (radiation_iact.h:222),
    # then its m_j * w_j product, the density loop's own m_j * w_j product
    # over different operands, and the single hi_inv_dim scaling of
    # enrichment_weight (GEAR_thermal/feedback.c:366). The kernel evaluation
    # itself adds no term here; see this function's docstring for why
    # identical compilation of W is a premise and cannot be given one.
    reconstruction = 4.0 * u32
    return {
        "dt_float32": dt_float32,
        "dt_text": dt_text,
        "closure": closure,
        "reconstruction": reconstruction,
        "bar": dt_float32 + dt_text + closure + reconstruction,
    }


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
    delta_t_token = ""
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
            # A zero step is the step-0 row, whose printed size is not a step
            # any particle took; matching it would put a zero in the reference
            # and turn the identity into 0/0.
            if not np.isfinite(dt) or dt <= 0.0:
                continue
            distance = abs(t - last["time_internal"])
            if distance < best and distance <= 0.5 * dt:
                best = distance
                delta_t = dt
                delta_t_token = fields[4]
    if delta_t is None:
        raise RuntimeError("No step in the log matches the last snapshot's time")
    ok = True
    print(
        f"injection: cosmological={last['cosmological']}, a {last['a']:.6g}, "
        f"Delta_t {delta_t:.6e} internal"
    )
    if not dusty:
        if last.get("n_stars", 0) != 1:
            print(
                f"  FAIL: the identity is written for a single illuminating "
                f"star, and reads PELuminosities[0] alone; this snapshot has "
                f"{last.get('n_stars', 0)}"
            )
            return False
        for name in ("L_PE", "L_LW"):
            if not np.isfinite(last[name]) or last[name] <= 0.0:
                print(f"  FAIL: {name} must be finite and positive, got {last[name]}")
                return False
        for band in ["PE", "LW"]:
            u_band = last[f"u_{band}"]
            bad = int(np.sum(~np.isfinite(u_band)))
            if bad:
                print(f"  FAIL: {bad} non-finite u_{band} value(s) in the snapshot")
                return False
            n_lit = int(np.sum(u_band > 0.0))
            if n_lit == 0:
                print(f"  FAIL: no gas particle carries {band}-band energy")
                return False
            lhs = np.sum(last["mass"] * u_band) / (
                last["mass_unit"] * last["energy_unit"]
            )
            rhs = delta_t * last[f"L_{band}"]
            terms = injection_identity_bar(delta_t_token, n_lit)
            print(
                f"  {band}: {n_lit} illuminated, Delta_t printed as "
                f"'{delta_t_token}'; bar terms: float32 Delta_t "
                f"{terms['dt_float32']:.2e}, log text {terms['dt_text']:.2e}, "
                f"weight-sum closure {terms['closure']:.2e}, reconstruction "
                f"{terms['reconstruction']:.2e}"
            )
            residual = abs(lhs / rhs - 1.0)
            ok &= gate(
                f"sum m u_{band} / (Delta_t L) - 1",
                residual,
                terms["bar"],
            )
            # Not gated, and deliberately worded so that a PASS/FAIL grep
            # cannot pick it up. The residual sits far inside the bar, so its
            # fraction of the bar is what makes a drift visible across runs.
            fraction = "n/a"
            if (
                np.isfinite(residual)
                and np.isfinite(terms["bar"])
                and terms["bar"] > 0.0
            ):
                fraction = f"{residual / terms['bar']:.4g}"
            print(f"  report only, not gated: residual/bar for {band} = {fraction}")
        return ok

    mechanism, path_in_kernel_radii = read_extinction_path(
        opt.snapshots, opt.extinction_path
    )
    path_cgs = extinction_path_cgs(
        last, mechanism, path_in_kernel_radii, opt.kernel_gamma
    )
    tau = optical_depths(last, path_cgs, opt.dust_to_gas_ratio)
    column_label = (
        "star-to-particle separation"
        if mechanism == "pair_separation"
        else f"{path_in_kernel_radii:g} kernel radii, kernel_gamma "
        f"{opt.kernel_gamma:g}"
    )
    print(
        f"  column: {mechanism} ({column_label}), path {np.min(path_cgs):.6e} "
        f"to {np.max(path_cgs):.6e} cm, local_dust_to_gas_ratio "
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
        f"float32 budget {dust_band_ratio_bar(signal):.2e}"
    )
    worst = float(np.max(np.abs(residual)))
    derived_bar = dust_band_ratio_bar(signal)
    if opt.dust_tol is None:
        ok &= gate(
            "band-ratio residual (G2), max_j |R_j|, derived bar",
            worst,
            derived_bar,
        )
    elif opt.dust_tol < 0.0:
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
    check_run_lw_calibration(
        opt.snapshots, opt.reference, opt.dark, opt.reference_dark, log=opt.log
    )
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
