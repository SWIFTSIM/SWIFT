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
Compare the simulated growth of the HII region radius against the classical
analytic Stromgren-sphere solution (equilibrium radius + Spitzer 1978 D-type
expansion).

This is a physics-correctness check, not just a "no crash" check: a
task-graph bug that silently truncates the ionization search radius still
produces a smooth, monotonic, plausible-looking r_hii(t) curve -- it would
just be quantitatively wrong. Comparing against the analytic solution catches
that a raw "did it crash" check cannot.

The ionizing photon rate Q_H(mass) is read from this run's own
Data/Radiation table (GEARFeedback:yields_table) via
../radiation_table_reader.py, using the same log-log
interpolate_1d()-based scheme src/feedback/GEAR/radiation.c's
radiation_get_ionization_rate_from_raw() actually runs, so this script
isn't guessing at a different photon budget than the code itself used.

The ever-tagged r_hii is outlier-rejected, not the max, of ever-tagged gas
radii: HIIStarIDs is a one-way stamp (never cleared on tag expiry), so a
single particle that lapsed early and then advected outward for several
Myr can inflate a plain max() by double digits percent.
robust_ever_tagged_radius() rejects such a particle only when it is
detached from the population both radially and in its own direction, which
leaves a smooth shell edge exact and keeps a genuinely anisotropic front
(a leading HEALPix cone at HII_angular_nside > 0, a breakout finger)
intact. The true max is still reported alongside it as a secondary
diagnostic.

Composition is treated the same way as T_i. The Spitzer solution is
parameterised by n_H, T_i, Q_H and the ionized mean molecular weight, none
of which is tied to the classical setup, so c_s is evaluated at the run's
OWN measured mu_i rather than at the pure-hydrogen mu_i=0.5 that gives the
familiar 12.85 km/s. The neutral mu is measured and printed too, but this
curve has no T_o term and does not use it. See resolve_composition().

Usage:
    python3 stromgren_analytic_check.py [-s snap/snapshot] [-o out.png]
"""

import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np
from astropy import constants as const
from astropy import units as u

try:
    import swiftsimio as sw
except ImportError:
    sw = None

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from radiation_table_reader import ionizing_photon_rate, load_radiation_table


# -----------------------------------------------------------------------------
# Analytic Stromgren-sphere solution
# -----------------------------------------------------------------------------
# This code's own ionized-gas temperature floor is metallicity-dependent
# (min(Delta_u_ionized, u_collisional), src/cooling/grackle/
# cooling_gear_subgrid.h), unlike the classical D-type solution's flat
# 1e4 K assumption. At this example's default Z=0, the applied temperature
# is ~47,500 K (Delta_u_ionized binds; see theory/GEAR/Radiation/,
# Section "Metallicity is not the explanation..."), not 1e4 K -- both
# alpha_B (now temperature-dependent since the Hui & Gnedin 1997 fix,
# src/feedback/GEAR/radiation.c) and the isothermal sound speed driving
# the D-type expansion should be evaluated there, not at the classical
# 1e4 K reference point, or this analytic comparison is systematically
# biased low.
def alpha_b_hui_gnedin(T):
    """Case-B recombination coefficient (Hui & Gnedin 1997, MNRAS 292, 27,
    Appendix A), the same fit now used in
    radiation_get_case_b_recombination_coefficient_cgs
    (src/feedback/GEAR/radiation.c)."""
    lam = 315614.0 / T.to(u.K).value
    return (
        2.753e-14 * lam**1.5 / (1.0 + (lam / 2.740) ** 0.407) ** 2.242 * u.cm**3 / u.s
    )


def sound_speed_ionized(T_ionized, mu_ionized):
    """Isothermal sound speed of the photoionized gas, c_s = sqrt(kT/(mu m_p)).

    Evaluated at the run's own ionized mean molecular weight rather than at
    the classical pure-hydrogen mu = 0.5 (which gives the familiar 12.85
    km/s at 1e4 K), for the reason given in resolve_composition()."""
    return (np.sqrt(const.k_B * T_ionized / (mu_ionized * const.m_p))).to(u.km / u.s)


def stromgren_radius(Q_H, n_H, alpha_B):
    """Equilibrium (Stromgren) radius: R_st = (3 Q_H / (4 pi alpha_B n_H^2))^(1/3)."""
    R_st3 = 3.0 * Q_H / (4.0 * np.pi * alpha_B * n_H**2)
    return R_st3.to(u.cm**3) ** (1.0 / 3.0)


def dtype_expansion_radius(t, R_st, c_s):
    """Spitzer (1978) D-type expansion: R(t) = R_st (1 + 7 c_s t / (4 R_st))^(4/7)."""
    return R_st * (1.0 + 7.0 * c_s * t / (4.0 * R_st)) ** (4.0 / 7.0)


# -----------------------------------------------------------------------------
# Robust ever-tagged front radius. This block is duplicated verbatim in
# StromgrenSphere/stromgren_analytic_check.py,
# Starbench/starbench_analytic_check.py,
# Hu2017/hu_smith_analytic_check.py and
# StromgrenSphereCosmo/cosmo_stromgren_analytic_check.py. The examples share
# no import path, so the four copies are kept identical by hand -- change one,
# change all four.
# -----------------------------------------------------------------------------

# A candidate is radially detached only if the gap down to the next-lower
# ever-tagged distance exceeds this many times the local radial spacing.
GAP_REJECTION_FACTOR = 5.0
# Fraction of the ever-tagged particles, counted in from the front, whose
# radial gaps set that local spacing. It is measured at the front rather than
# over the whole population because an anisotropic region's outer radii are
# populated by its leading directions alone, and are therefore genuinely
# sparser than its interior.
GAP_SCALE_TOP_FRACTION = 0.05
# ... but never from fewer than this many gaps. The median of two gaps is
# their mean, so "one gap exceeds GAP_REJECTION_FACTOR times the median of
# the top two" cannot be satisfied by non-negative gaps, and the whole test
# switches itself off below ~40 ever-tagged particles.
GAP_SCALE_MIN_GAPS = 5
# Half-angle of the cone searched for angular support. 30 deg is roughly one
# HEALPix nside=1 pixel (4*pi/12 sr is a 33.6 deg cap), the coarsest angular
# structure HII_angular_nside can impose on the front.
ANGULAR_SUPPORT_HALFANGLE_DEG = 30.0
# A candidate is angularly supported when another retained ever-tagged
# particle in its cone reaches at least this fraction of its own radius.
ANGULAR_SUPPORT_RADIUS_FRACTION = 0.95
# Ceiling on how many particles may be rejected, so a mis-estimated spacing
# can never consume a large part of the front.
GAP_MAX_REJECT_FRACTION = 0.02
GAP_MIN_REJECT = 3
GAP_MAX_REJECT = 32
# Below this many ever-tagged particles there are too few points to estimate
# a spacing scale at all; trust the max unmodified.
GAP_MIN_POINTS = 10


def robust_ever_tagged_radius(offsets_kpc):
    """Robust "ever-tagged" HII front radius, rejecting particles that are
    isolated both radially and in their own direction.

    HIIStarIDs is a one-way stamp (radiation_reset_part_ionized_tag never
    clears it, the stamp intentionally marks the swept shell), so a particle
    whose tag lapsed early and then advected outward for several Myr can
    dominate a plain max() and inflate the reported extent by double digits
    percent. A blind percentile trim removes such stragglers too, but also
    cuts several percent into a perfectly smooth, gap-free shell edge where
    every particle is legitimate.

    A purely radial gap test cannot tell a straggler from the leading edge of
    an anisotropic front: at HII_angular_nside > 0 the region advances at a
    different rate in each HEALPix cone, so the outermost cone's tip is
    radially separated from every other cone's gas and a radial test alone
    eats it. The two cases differ in direction, not in radius. A leading edge
    has companions at a comparable radius inside its own cone; an advected
    straggler has none anywhere near it. Both conditions are therefore
    required before anything is rejected:

      1. radial detachment: the gap down to the next ever-tagged particle
         exceeds GAP_REJECTION_FACTOR times the local radial spacing (the
         median gap over the front-most GAP_SCALE_TOP_FRACTION of the
         points, taken over no fewer than GAP_SCALE_MIN_GAPS gaps);
      2. angular isolation: no retained ever-tagged particle within
         ANGULAR_SUPPORT_HALFANGLE_DEG of the candidate's direction reaches
         ANGULAR_SUPPORT_RADIUS_FRACTION of the candidate's radius.

    The largest cut satisfying both is taken, so a small clump of stragglers
    that advected together, and that would shield each other from a
    stop-at-the-first-survivor walk, is rejected as a group. At most
    max(GAP_MIN_REJECT, GAP_MAX_REJECT_FRACTION * n) particles are removed,
    and never more than GAP_MAX_REJECT.

    Assumes the offsets are already minimum-imaged about the source, so the
    directions and radii are the physical ones.

    @param offsets_kpc (N, 3) ever-tagged particle offsets from the source
    (kpc), minimum-imaged.
    @return (front radius after outlier rejection, true max), both in kpc.
    """
    offsets = np.asarray(offsets_kpc, dtype=np.float64)
    radii = np.linalg.norm(offsets, axis=1)
    order = np.argsort(radii)[::-1]
    r_sorted = radii[order]
    true_max = float(r_sorted[0])
    n = len(r_sorted)
    if n < GAP_MIN_POINTS:
        return true_max, true_max

    gaps = -np.diff(r_sorted)
    n_top = min(
        max(GAP_SCALE_MIN_GAPS, int(np.ceil(GAP_SCALE_TOP_FRACTION * n))), n - 1
    )
    local_scale = float(np.median(gaps[:n_top]))
    if local_scale <= 0.0:
        positive = gaps[gaps > 0.0]
        local_scale = float(np.median(positive)) if len(positive) else 0.0
    if local_scale <= 0.0:
        return true_max, true_max
    threshold = GAP_REJECTION_FACTOR * local_scale

    max_reject = min(
        max(GAP_MIN_REJECT, int(np.ceil(GAP_MAX_REJECT_FRACTION * n))),
        GAP_MAX_REJECT,
        n - 1,
    )

    # Directions, in the same descending-radius order. A particle sitting
    # exactly on the source has no direction; leave its unit vector at zero so
    # it never enters anyone's cone.
    unit = np.zeros((n, 3))
    has_direction = r_sorted > 0.0
    unit[has_direction] = offsets[order][has_direction] / r_sorted[has_direction, None]
    cos_halfangle = np.cos(np.radians(ANGULAR_SUPPORT_HALFANGLE_DEG))
    # cone_ranks[k]: descending-radius ranks of the particles inside candidate
    # k's cone, ascending, so the first entry past a cut is the largest radius
    # that cut retains in that direction.
    cone_ranks = []
    for k in range(max_reject):
        in_cone = (unit @ unit[k]) >= cos_halfangle
        in_cone[k] = False
        cone_ranks.append(np.flatnonzero(in_cone))

    cut = -1
    for i in range(max_reject):
        if gaps[i] <= threshold:
            continue
        supported = False
        for k in range(i + 1):
            ranks = cone_ranks[k]
            first_retained = np.searchsorted(ranks, i, side="right")
            if (
                first_retained < len(ranks)
                and r_sorted[ranks[first_retained]]
                >= ANGULAR_SUPPORT_RADIUS_FRACTION * r_sorted[k]
            ):
                supported = True
                break
        if not supported:
            cut = i
    if cut < 0:
        return true_max, true_max
    return float(r_sorted[cut + 1]), true_max


def snapshot_index(filename):
    """Trailing integer of a snapshot filename (snapshot_0007.hdf5 -> 7), or
    None when the stem carries no such index."""
    stem = os.path.splitext(os.path.basename(filename))[0]
    match = re.search(r"(\d+)$", stem)
    return int(match.group(1)) if match else None


def print_n_H_with_source(n_H, source, files, label="n_H"):
    """Print n_H together with the snapshot it was measured in, and warn when
    that snapshot is not the run's own first one.

    n_H is read once, from the first snapshot of the glob that holds a live
    star, and never revisited. A glob that skips the run's early snapshots
    samples a post-expansion, shell-biased density instead of the initial
    uniform one, which shifts every reference curve computed from it -- so
    the source has to be visible in the output rather than implicit.
    """
    print(f"{label:<19}: {n_H:.4g}  (measured in {source})")
    if source != files[0]:
        print(
            f"  NOTE: not the earliest snapshot in the glob ({files[0]}); "
            f"that one holds no live star yet."
        )
    first_index = snapshot_index(files[0])
    if first_index is not None and first_index != 0:
        print(
            f"  WARNING: the snapshot glob starts at index {first_index}, not "
            f"0. n_H is then measured after the region has already expanded, "
            f"so it is a shell-biased density and every reference curve below "
            f"is shifted. Re-run against the complete snapshot series."
        )


# -----------------------------------------------------------------------------
# Read the simulated r_hii(t) from the snapshots
# -----------------------------------------------------------------------------
def read_simulated_r_hii(snapshot_glob):
    if sw is None:
        raise RuntimeError("swiftsimio is required to read the snapshots.")

    files = sorted(glob.glob(snapshot_glob))
    if not files:
        raise RuntimeError(f"No snapshots found matching {snapshot_glob!r}.")

    times = []
    r_hii = []
    r_hii_max = []
    r_now = []
    star_mass_msun = None
    n_H_atom_cc = None
    n_H_source = None
    boxsize_kpc = None

    for f in files:
        data = sw.load(f)
        if boxsize_kpc is None:
            # Smallest box dimension -- the periodic image distance that
            # actually bounds how large r_hii can grow before the HII region
            # starts interacting with its own periodic copies.
            boxsize_kpc = float(np.min(data.metadata.boxsize.to("kpc").value))
        if len(data.stars.hiiregion_radii) == 0:
            continue

        # swiftsimio/unyt quantities can't be combined directly with
        # astropy.units quantities -- convert to plain floats in a fixed
        # unit at this boundary, then re-wrap as astropy Quantities.
        times.append(data.metadata.time.to("Myr").value)

        # h_hii of the star. Unlike the ever-tagged extent this returns to
        # exactly 0 when the star dies or ages past HII_max_age
        # (feedback_common.c retires it to the tracers there and nowhere
        # else), so it, not r_hii, is what tells the comparison window when
        # the source went dark.
        r_now.append(float(np.max(data.stars.hiiregion_radii).to("kpc").value))

        # r_hii from the star's own HIIRegionRadii (h_hii) is frozen at the
        # moment each gas particle was tagged -- it never updates as that
        # particle keeps moving afterwards. Since ionized gas has substantial
        # outward velocity (confirmed directly: tagged particles show mean
        # v_r several times larger than untagged gas, growing with time),
        # h_hii systematically undercounts the true current extent of the
        # ionized/expanding shell -- by up to ~60-80% at t~2-3 Myr in the
        # e10 check. Instead, use the current position of every gas particle
        # this star has *ever* tagged: HIIStarIDs is set once when a
        # particle is first ionized and, unlike IsIonizedFlags, is never
        # reset when that tag later expires (see
        # radiation_reset_part_ionized_tag, src/feedback/GEAR/radiation.c),
        # so it survives the rebuild-cadence tag/re-tag cycle. A plain max()
        # over these positions is fragile to a single lapsed-tag particle
        # that keeps advecting outward, so robust_ever_tagged_radius gap-
        # rejects it and reports the front position as the primary value,
        # with the true max as a secondary diagnostic.
        star_id = int(data.stars.particle_ids[0])
        ever_tagged = data.gas.hiistar_ids == star_id
        if np.any(ever_tagged):
            star_pos = data.stars.coordinates[0].to("kpc").value
            gas_pos = data.gas.coordinates[ever_tagged].to("kpc").value
            box_kpc = data.metadata.boxsize.to("kpc").value
            # Minimum-image convention: a periodic box means the naive
            # difference can be wrong by a whole box length once gas has
            # wrapped around, which would corrupt the extent with a single
            # spurious huge distance -- match the simulation's own
            # nearest-image distance convention here.
            dx = gas_pos - star_pos
            dx -= box_kpc * np.round(dx / box_kpc)
            r_robust, r_max = robust_ever_tagged_radius(dx)
            r_hii.append(r_robust)
            r_hii_max.append(r_max)
        else:
            r_hii.append(r_now[-1])
            r_hii_max.append(r_now[-1])

        if star_mass_msun is None:
            n_stars = len(data.stars.masses)
            if n_stars != 1:
                raise RuntimeError(
                    f"Expected exactly 1 star (this check assumes "
                    f"star_type=single_star, matching the Q_H(mass) fit "
                    f"reimplemented below), found {n_stars}. Re-run with "
                    f"star_type=single_star, n_stars=1, or fix this script's "
                    f"Q_H fit for whatever star_type was actually used before "
                    f"trusting this comparison."
                )
            star_mass_msun = float(np.max(data.stars.masses).to("Msun").value)
        if n_H_atom_cc is None:
            # Hydrogen number density from the gas density and the
            # Grackle-reported hydrogen mass fraction. Read the actual value
            # this run used (GrackleCooling:HydrogenFractionByMass, stored in
            # every snapshot's Parameters group) rather than assuming a
            # fixed constant -- a previous version hardcoded 0.716 here,
            # silently disagreeing with e.g. this example's actual 0.76 and
            # biasing n_H (and R_st ~ n_H^{-2/3}) by a few percent. Only fall
            # back to the primordial default if the run used a different
            # cooling backend that doesn't expose this parameter.
            rho_g_cm3 = float(np.mean(data.gas.densities).to("g/cm**3").value)
            X_H_raw = data.metadata.parameters.get(
                "GrackleCooling:HydrogenFractionByMass"
            )
            X_H = float(X_H_raw) if X_H_raw is not None else 0.716
            n_H_atom_cc = (rho_g_cm3 * u.g / u.cm**3 * X_H / const.m_p).to(1 / u.cm**3)
            n_H_source = f

    return (
        u.Quantity(times, u.Myr),
        u.Quantity(r_hii, u.kpc),
        u.Quantity(r_hii_max, u.kpc),
        u.Quantity(r_now, u.kpc),
        star_mass_msun,
        n_H_atom_cc,
        n_H_source,
        boxsize_kpc * u.kpc,
        files,
    )


def measure_ionized_temperature_K(files, min_count=20):
    """Median temperature of the currently-tagged (ionized) gas, from the
    latest snapshot holding at least min_count tagged particles.

    The run's actual T_i depends on build flags
    (IONIZATION_FEEDBACK_DEBUG_FIXED_IONIZED_TEMPERATURE_K) and
    metallicity; measuring it removes the silent tens-of-percent error of
    computing the reference at a temperature the run never used.

    mu is returned alongside so the sound speeds can be evaluated at the same
    composition this temperature was derived from.

    @return (median T [K], median mu, mu came from species arrays, source
    filename), or (None, None, False, None) if no snapshot holds enough
    tagged gas.
    """
    for filename in reversed(files):
        with h5py.File(filename, "r") as h:
            gas = h["PartType0"]
            ionized = gas["IsIonizedFlags"][:].astype(bool)
            if np.sum(ionized) < min_count:
                continue
            u_L = h["Units"].attrs["Unit length in cgs (U_L)"][0]
            u_t = h["Units"].attrs["Unit time in cgs (U_t)"][0]
            u_to_cgs = (u_L / u_t) ** 2  # internal energy unit -> erg/g
            if "HI" in gas:
                inv_mu = (
                    gas["HI"][ionized]
                    + 2.0 * gas["HII"][ionized]
                    + 0.25 * gas["HeI"][ionized]
                    + 0.5 * gas["HeII"][ionized]
                    + 0.75 * gas["HeIII"][ionized]
                )
                mu = 1.0 / np.clip(inv_mu, 1e-10, None)
                mu_measured = True
            else:
                # No species arrays (COOLING_GRACKLE_MODE 0): this selection
                # is fully ionized by construction, so use the fully-ionized
                # primordial mu.
                mu = MU_IONIZED_FALLBACK
                mu_measured = False
            u_cgs = gas["InternalEnergies"][ionized] * u_to_cgs
            gamma_m1 = 2.0 / 3.0
            T = u_cgs * gamma_m1 * mu * const.m_p.cgs.value / const.k_B.cgs.value
            return (
                float(np.median(T)),
                float(np.median(mu)) if mu_measured else mu,
                mu_measured,
                filename,
            )
    return None, None, False, None


# -----------------------------------------------------------------------------
# Composition the reference curves are evaluated at. Duplicated verbatim in all
# four analytic checks, same rule as the estimator block above: change one,
# change all four.
#
# The Raga/Spitzer/Bisbas solutions are parameterised by n_H, T_i, T_o, Q_H and
# the mean molecular weights. None of that is tied to the source papers' own
# setup, so the reference is evaluated at THIS run's own composition rather
# than at the papers' pure-hydrogen convention (mu_i = 0.5 for fully ionized
# hydrogen, mu_o = 1 for neutral atomic hydrogen), which this run does not
# have: it carries helium, and Grackle sets the ionization state.
#
# Fallbacks, used only when a run carries no species arrays
# (COOLING_GRACKLE_MODE 0) and marked in the verdict token when they are:
# primordial X=0.76, Y=0.24, fully ionized (1/mu = 2X + 3Y/4) and fully
# neutral (1/mu = X + Y/4) respectively.
# -----------------------------------------------------------------------------
MU_IONIZED_FALLBACK = 0.6
MU_NEUTRAL_FALLBACK = 1.0 / (0.76 + 0.24 / 4.0)


def measure_neutral_mu(files, min_count=20):
    """Median mean molecular weight of the gas no source has ever tagged.

    That selection is the undisturbed ambient medium the front expands into,
    which is exactly what c_o describes. Species mass fractions are
    dimensionless, so no comoving-to-physical conversion applies and this
    reads the same in a cosmological snapshot as in a plain one.

    @param files Snapshot filenames, in time order.
    @param min_count Fewest never-tagged particles a snapshot must hold.
    @return (median mu, mu came from species arrays, source filename). Falls
    back to (MU_NEUTRAL_FALLBACK, False, None) when no snapshot qualifies or
    the run carries no species arrays.
    """
    for filename in reversed(files):
        with h5py.File(filename, "r") as h:
            gas = h["PartType0"]
            if "HI" not in gas:
                return MU_NEUTRAL_FALLBACK, False, None
            star_ids = (
                h["PartType4"]["ParticleIDs"][:] if "PartType4" in h else np.array([])
            )
            never = ~np.isin(gas["HIIStarIDs"][:], star_ids)
            if np.sum(never) < min_count:
                continue
            inv_mu = (
                gas["HI"][:][never]
                + 2.0 * gas["HII"][:][never]
                + 0.25 * gas["HeI"][:][never]
                + 0.5 * gas["HeII"][:][never]
                + 0.75 * gas["HeIII"][:][never]
            )
            mu = 1.0 / np.clip(inv_mu, 1e-10, None)
            return float(np.median(mu)), True, filename
    return MU_NEUTRAL_FALLBACK, False, None


def resolve_composition(files, uses_neutral):
    """Return (mu_i, mu_o, verdict marker) the sound speeds are evaluated at.

    Both are measured from the run's own species arrays, the same arrays the
    temperature measurement already uses to convert internal energy to T, so
    c_i, c_o and T are all read off one composition instead of mixing a
    measured T with an assumed mu. See the block comment above.

    @param files Snapshot filenames, in time order.
    @param uses_neutral Whether this check's curve contains a T_o term. The
    Spitzer D-type solution does not, so mu_o is reported there for
    reference but must not reach the verdict token.
    @return (mu_i, mu_o, "" or a " (COMPOSITION FALLBACK: ...)" marker).
    """
    _, mu_i, mu_i_measured, source_i = measure_ionized_temperature_K(files)
    mu_o, mu_o_measured, source_o = measure_neutral_mu(files)
    if mu_i is None:
        mu_i, mu_i_measured, source_i = MU_IONIZED_FALLBACK, False, None

    def describe(measured, source):
        if measured:
            return f"measured from species, {source}"
        return "primordial fallback, no species arrays in this run"

    print(f"mu_ionized         : {mu_i:.4f} ({describe(mu_i_measured, source_i)})")
    unused = "" if uses_neutral else ", not used: this curve has no T_o term"
    print(
        f"mu_neutral         : {mu_o:.4f} "
        f"({describe(mu_o_measured, source_o)}{unused})"
    )

    fell_back = []
    if not mu_i_measured:
        fell_back.append("mu_i")
    if uses_neutral and not mu_o_measured:
        fell_back.append("mu_o")
    if not fell_back:
        return mu_i, mu_o, ""
    return (
        mu_i,
        mu_o,
        f" (COMPOSITION FALLBACK: {', '.join(fell_back)} assumed, not measured)",
    )


def resolve_T_ionized_K(requested_K, files):
    """Return (T_i [K] the reference is computed at, verdict marker).

    None (the default) means "measure it from the run's own tagged gas". An
    explicit value is honoured, but when it disagrees with the measured one
    by more than 10% the reference describes a temperature the run never
    held. The marker returned alongside is appended to every verdict token
    below: a banner printed at the top of the output does not survive a
    harness grepping for PASS/FAIL, so without it such a harness records a
    clean verdict reached against the wrong temperature.

    @param requested_K Explicit T_i [K], or None to measure it.
    @param files Snapshot filenames, in time order.
    @return (T_i [K], "" or a " (T_i MISMATCH: ...)" verdict marker).
    """
    T_measured_K, _, _, source = measure_ionized_temperature_K(files)
    if requested_K is None:
        if T_measured_K is None:
            raise RuntimeError(
                "Could not measure T_i (no snapshot holds enough tagged "
                "gas); pass --T-ionized-K explicitly."
            )
        print(
            f"T_ionized          : {T_measured_K:.4g} K "
            f"(measured from tagged gas, {source})"
        )
        return T_measured_K, ""
    if T_measured_K is not None and abs(T_measured_K / requested_K - 1.0) > 0.1:
        print("\n" + "!" * 72)
        print(
            f"WARNING: --T-ionized-K {requested_K:.4g} K disagrees with the "
            f"temperature this run actually holds its tagged gas at "
            f"({T_measured_K:.4g} K, measured in {source})."
        )
        print(
            "Every reference curve below is computed at a temperature the "
            "run did not use, so the verdict is not meaningful. Omit "
            "--T-ionized-K to use the measured value."
        )
        print("!" * 72 + "\n")
        return requested_K, (
            f" (T_i MISMATCH: reference {requested_K:.4g} K vs measured "
            f"{T_measured_K:.4g} K)"
        )
    return requested_K, ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--snapshot-glob",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for the snapshots to read (default: snap/snapshot_*.hdf5)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="stromgren_analytic_check.png",
        help="Output plot filename.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=0.3,
        help="Relative-error tolerance at the final time for a pass (default: 0.3).",
    )
    parser.add_argument(
        "--T-ionized-K",
        type=float,
        default=None,
        dest="T_ionized_K",
        help="Ionized-gas temperature alpha_B and c_s are evaluated at. "
        "Default: measured from the run's own tagged gas, which always "
        "matches whatever build flag/metallicity convention the run used. "
        "An explicit value overrides the measurement but warns loudly if "
        "the two disagree by >10%%.",
    )
    args = parser.parse_args()
    T_ionized_K, T_i_marker = resolve_T_ionized_K(
        args.T_ionized_K, sorted(glob.glob(args.snapshot_glob))
    )

    T_IONIZED_K = T_ionized_K * u.K
    ALPHA_B = alpha_b_hui_gnedin(T_IONIZED_K)
    mu_i, mu_o, mu_marker = resolve_composition(
        sorted(glob.glob(args.snapshot_glob)), uses_neutral=False
    )
    verdict_marker = T_i_marker + mu_marker
    C_S_IONIZED = sound_speed_ionized(T_IONIZED_K, mu_i)

    (
        t_sim,
        r_sim,
        r_sim_max,
        r_now,
        star_mass_msun,
        n_H,
        n_H_source,
        boxsize,
        files,
    ) = read_simulated_r_hii(args.snapshot_glob)

    if star_mass_msun is None or n_H is None:
        raise RuntimeError("Could not infer star mass / n_H from the snapshots.")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    radiation_table = load_radiation_table(files, script_dir)
    Q_H = ionizing_photon_rate(radiation_table, star_mass_msun)
    R_st = stromgren_radius(Q_H, n_H, ALPHA_B)
    box_half_width = 0.5 * boxsize

    print(f"Star mass          : {star_mass_msun:.3f} Msun")
    print_n_H_with_source(n_H, n_H_source, files)
    print(f"Q_H                : {Q_H:.4g}")
    print(
        f"T_ionized (applied): {T_IONIZED_K:.4g}  ->  alpha_B = {ALPHA_B:.4g}, "
        f"c_s = {C_S_IONIZED:.4g} (mu_i={mu_i:.4f})"
    )
    print(f"Equilibrium R_st   : {R_st.to(u.pc):.4g} = {R_st.to(u.kpc):.4g}")
    print(f"Box half-width     : {box_half_width.to(u.pc):.4g}")
    if R_st > box_half_width:
        print(
            "  WARNING: R_st exceeds the box half-width -- the periodic box "
            "is too small to contain a full infinite-medium Stromgren "
            "sphere. The analytic curve below is an upper bound, not a "
            "value this finite, periodic setup can actually reach; a "
            "mismatch alone does not indicate a code bug."
        )

    r_analytic = dtype_expansion_radius(t_sim, R_st, C_S_IONIZED).to(u.kpc)

    # The comparison is only meaningful while (a) the analytic curve is
    # still predicted to fit inside the box -- beyond that, periodic
    # self-interaction/wrapping makes the simulated and infinite-medium
    # analytic solutions incomparable, not wrong -- and (b) the star is
    # still alive and actively tracked: h_hii resets to 0 once a star
    # dies or ages past HII_max_age (feedback_will_do_feedback),
    # which is bookkeeping, not the region physically vanishing. Picking
    # the LAST snapshot that satisfies both, rather than blindly using
    # the final simulated time, makes this check robust to runs that
    # continue well past either limit (e.g. an acceptance run taken to
    # time_end long after the star has died).
    # Liveness comes from h_hii (r_now), not from the ever-tagged extent:
    # the HIIStarIDs stamp is never cleared, so r_sim keeps growing by
    # advection after the star dies and never returns to 0, which made the
    # old (r_sim > 0) test unable to ever detect death.
    valid = (r_analytic <= box_half_width) & (r_now > 0)
    box_exceeded_at = t_sim[r_analytic > box_half_width]
    # r_now == 0 before the first HII rebuild is not "death" -- only flag a
    # zero that follows a nonzero sample as the star going dark.
    nonzero_idx = np.where(r_now > 0)[0]
    dead_mask = np.zeros(len(r_now), dtype=bool)
    if len(nonzero_idx) > 0:
        dead_mask[nonzero_idx[0] + 1 :] = r_now[nonzero_idx[0] + 1 :] == 0
    star_dead_at = t_sim[dead_mask]
    if not np.any(valid):
        raise RuntimeError(
            "No snapshot has both an analytic radius within the box and an "
            "alive, actively-tracked r_sim > 0 -- no meaningful comparison "
            "point exists in this run. Try a bigger boxsize or a higher "
            "gas_density (shrinks R_st ~ n_H^-2/3), or check the run "
            "completed at least one HII rebuild before the star died."
        )
    last_valid = np.where(valid)[0][-1]
    t_sim_final = t_sim[last_valid]
    r_sim_final = r_sim[last_valid]
    r_sim_max_final = r_sim_max[last_valid]
    r_analytic_final = r_analytic[last_valid]
    rel_error = float(
        (np.abs(r_sim_final - r_analytic_final) / r_analytic_final).decompose()
    )
    verdict = "PASS" if rel_error <= args.tol else "FAIL"
    print(
        f"Comparison time: t={t_sim_final:.4g}  r_sim={r_sim_final:.4g}  "
        f"r_analytic={r_analytic_final:.4g}  rel_error={rel_error:.2%}  "
        f"[{verdict}{verdict_marker}]"
    )
    print(
        f"  (secondary diagnostic, true max ever-tagged extent: "
        f"{r_sim_max_final:.4g})"
    )
    if last_valid < len(t_sim) - 1:
        reason = []
        if len(box_exceeded_at) > 0:
            reason.append(f"box half-width exceeded from t={box_exceeded_at[0]:.4g}")
        if len(star_dead_at) > 0:
            reason.append(f"star inactive/dead from t={star_dead_at[0]:.4g}")
        print(
            f"  (not the final simulated time t={t_sim[-1]:.4g} -- "
            f"{'; '.join(reason)}; later snapshots are not comparable to "
            f"the infinite-medium analytic solution.)"
        )

    fig, ax = plt.subplots()
    ax.plot(t_sim.to(u.Myr), r_sim.to(u.pc), "o-", label="Simulation ($r_{hii}$)")
    ax.plot(
        t_sim.to(u.Myr),
        r_analytic.to(u.pc),
        "--",
        label="Analytic (Spitzer 1978 D-type)",
    )
    ax.axhline(
        R_st.to(u.pc).value, color="grey", linestyle=":", label="Equilibrium $R_{st}$"
    )
    ax.axhline(
        box_half_width.to(u.pc).value,
        color="firebrick",
        linestyle=":",
        label="Box half-width",
    )
    ax.axvline(
        t_sim_final.to(u.Myr).value,
        color="black",
        linestyle="-.",
        alpha=0.6,
        label="Comparison point",
    )
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("HII region radius [pc]")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
