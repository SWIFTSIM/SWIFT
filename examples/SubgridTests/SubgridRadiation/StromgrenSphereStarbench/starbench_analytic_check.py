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
Compare the simulated growth of the HII region radius against the STARBENCH
semi-empirical D-type expansion formula (Bisbas et al. 2015, MNRAS 453, 1324,
arXiv:1507.05621), Eqns 8, 11, 28 and 29.

Unlike the classical Spitzer (1978) solution used by
stromgren_analytic_check.py, STARBENCH is a numerically-integrated blend of
two ODEs, tuned against 12 independent radiation-hydrodynamics codes to
agree with high-resolution simulations to <~2%:

  Eqn 8  (Raga-I):  the exact pressure-driven thin-shell ODE (Spitzer's
                     closed-form solution is what you get by dropping the
                     small mu_i*T_o/(mu_o*T_i) term from this).
  Eqn 11 (Raga-II): Hosokawa & Inutsuka (2006)'s ODE, including the inertia
                     of the swept-up shell (2nd order in R).
  Eqn 28: R_SB = R_II + f_SB * (R_I - R_II)
  Eqn 29: f_SB = 1 - 0.733 * exp(-t / 1 Myr)

This is the STARBENCH paper's own late-phase test density and source
(rho_o=5.21e-21 g/cm3, Q_H=1e49/s) at this code's own Z=0 ionized/neutral
temperature floors (T_i~47500 K, T_o~1e3 K via SPH:minimal_temperature),
not the papers' flat T_i=1e4 K -- the comparison is evaluated at whatever
T_i/T_o this run actually used (--T-ionized-K/--T-neutral-K), not
hardcoded to the papers' own numbers. See this example's README for the
full parameter derivation and why Z=0 (not the Z/Zsun~0.231 workaround
used elsewhere in this project) is the right choice here.

Two r_hii conventions are reported, because the photon-count budget can
genuinely recede a region and the two then diverge:

  ever-tagged extent: the current position of every gas particle the star
    has *ever* tagged (HIIStarIDs) -- a historical-extent measure. The
    right front proxy for the default flag-and-floor scheme, where the
    dense swept shell keeps its ownership stamp after tag expiry.
  instantaneous front: the star's own current HIIRegionRadii (h_hii),
    per snapshot -- the distance the budget maintains *right now*. The
    right measure for rate-coupled (HII_couple_ionization_rate) runs,
    where it coincides with the actual x_HII front; the ever-tagged
    extent instead freezes at the historical maximum when the region
    recedes.

The PASS/FAIL verdict is attached to each line separately; when the two
conventions diverge by more than 10% the region receded and the
scheme-appropriate line is the meaningful one (see above).

The ever-tagged extent itself is outlier-rejected, not the max, of
ever-tagged gas radii: HIIStarIDs is a one-way stamp (never cleared on tag
expiry), so a single particle that lapsed early and advected outward for
several Myr can inflate a plain max() by double digits percent.
robust_ever_tagged_radius() rejects such a particle only when it is
detached from the population both radially and in its own direction, which
leaves a smooth shell edge exact and keeps a genuinely anisotropic front
(a leading HEALPix cone at HII_angular_nside > 0, a breakout finger)
intact. The true max is still reported as a secondary diagnostic.

Usage:
    python3 starbench_analytic_check.py [-s snap/snapshot] [-o out.png]
"""

import argparse
import glob
import os
import re

import numpy as np
from astropy import constants as const
from astropy import units as u
from scipy.integrate import solve_ivp

try:
    import swiftsimio as sw
except ImportError:
    sw = None

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def alpha_b_hui_gnedin(T):
    """Case-B recombination coefficient (Hui & Gnedin 1997, MNRAS 292, 27,
    Appendix A), the same fit used in
    radiation_get_case_b_recombination_coefficient_cgs
    (src/feedback/GEAR/radiation.c). Copied verbatim from
    stromgren_analytic_check.py."""
    lam = 315614.0 / T.to(u.K).value
    return (
        2.753e-14 * lam**1.5 / (1.0 + (lam / 2.740) ** 0.407) ** 2.242 * u.cm**3 / u.s
    )


# -----------------------------------------------------------------------------
# STARBENCH semi-empirical formula (Bisbas et al. 2015)
# -----------------------------------------------------------------------------
def starbench_curve(
    t_myr, R_St_pc, c_i_km_s, c_o_km_s, mu_i=0.5, mu_o=1.0, T_i_K=1e4, T_o_K=1e3
):
    """Solve Eqns 8 and 11 numerically and blend via Eqns 28-29.

    mu_i=0.5 (fully ionized pure H, proton+electron) and mu_o=1.0 (neutral
    atomic H) are the paper's implicit values -- confirmed by cross-checking
    against its own stated "ratio ~ 1/200" for the early-phase test
    (T_o=100 K) and its quoted c_i/c_o for the late-phase test; do not
    change these without re-deriving them the same way.

    @param t_myr Array of times (Myr) to evaluate at, must start > 0.
    @param R_St_pc Initial Stromgren radius (pc).
    @param c_i_km_s, c_o_km_s Ionized/neutral isothermal sound speeds (km/s).
    @return (R_I, R_II, R_SB) arrays in pc, same shape as t_myr.
    """
    pc_to_km = 3.0857e13
    Myr_to_s = 3.15576e13
    c_i = c_i_km_s * Myr_to_s / pc_to_km  # pc / Myr
    c_o = c_o_km_s * Myr_to_s / pc_to_km

    def rhs_raga1(t, y):
        R = y[0]
        ratio = (R_St_pc / R) ** 0.75
        term2 = (mu_i * T_o_K) / (mu_o * T_i_K) * (R_St_pc / R) ** (-0.75)
        return [c_i * (ratio - term2)]

    def rhs_raga2_conservative(t, y):
        # w = R^3 * dR/dt (proportional to the shell's physical momentum);
        # avoids the (3/R)(dR/dt)^2 term's numerical stiffness in the naive
        # R'' + (3/R)R'^2 = ... form.
        R, w = y
        dRdt = w / R**3
        dwdt = 3.0 * R_St_pc**1.5 * c_i**2 * np.sqrt(R) - 3.0 * c_o**2 * R**2
        return [dRdt, dwdt]

    t0 = min(t_myr.min(), 1e-4)
    t_span = [t0, t_myr.max()]

    sol1 = solve_ivp(rhs_raga1, t_span, [R_St_pc], t_eval=t_myr, rtol=1e-10, atol=1e-12)
    if sol1.status != 0:
        raise RuntimeError(f"Raga-I (Eqn 8) integration failed: {sol1.message}")

    # v(0) = c_i matches Raga-I's own initial expansion rate at R=R_St; a
    # larger initial velocity (e.g. sqrt(4/3)*c_i, Eqn 13's implied slope)
    # sends this ODE into an unstable oscillation that crashes through R=0
    # well before 3 Myr -- verified empirically, not just in principle.
    w0 = R_St_pc**3 * c_i
    sol2 = solve_ivp(
        rhs_raga2_conservative,
        t_span,
        [R_St_pc, w0],
        t_eval=t_myr,
        method="Radau",
        rtol=1e-10,
        atol=1e-12,
    )
    if sol2.status != 0:
        raise RuntimeError(f"Raga-II (Eqn 11) integration failed: {sol2.message}")

    R_I = sol1.y[0]
    R_II = sol2.y[0]
    f_SB = 1.0 - 0.733 * np.exp(-t_myr / 1.0)  # Eqn 29
    R_SB = R_II + f_SB * (R_I - R_II)  # Eqn 28
    return R_I, R_II, R_SB


# -----------------------------------------------------------------------------
# Robust ever-tagged front radius. This block is duplicated verbatim in
# StromgrenSphere/stromgren_analytic_check.py,
# StromgrenSphereStarbench/starbench_analytic_check.py,
# StromgrenSphereHuSmith2017/hu_smith_analytic_check.py and
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
# Q_H(mass): direct port of radiation_get_individual_star_ionizing_photon_
# emission_rate_fit() and the radius/luminosity fits it calls
# (src/feedback/GEAR/radiation.c). Copied verbatim from
# stromgren_analytic_check.py (not imported, to keep this example
# self-contained) -- keep both in sync if the fit ever changes.
# -----------------------------------------------------------------------------
def star_radius(mass_msun):
    """Empirical mass-radius relation (radiation_get_individual_star_radius)."""
    if mass_msun < 1.0:
        return (mass_msun**0.8) * const.R_sun
    elif mass_msun < 8.0:
        return (mass_msun**0.57) * const.R_sun
    else:
        return (mass_msun**0.5) * const.R_sun


def star_luminosity(mass_msun):
    """Empirical mass-luminosity relation (radiation_get_individual_star_luminosity)."""
    if mass_msun < 0.43:
        lum_sol = 0.185 * mass_msun**2
    elif mass_msun < 2.0:
        lum_sol = mass_msun**4
    elif mass_msun < 54.0:
        lum_sol = 1.5 * mass_msun**3.5
    else:
        lum_sol = 32000.0 * mass_msun
    return lum_sol * const.L_sun


def ionizing_photon_rate(mass_msun):
    """
    Q_H, the ionizing photon emission rate for a single star of the given
    mass, from a blackbody fit to its (R, L) -- direct port of
    radiation_get_individual_star_ionizing_photon_emission_rate_fit().
    """
    R = star_radius(mass_msun)
    L = star_luminosity(mass_msun)
    if R <= 0 or L <= 0:
        return 0.0 / u.s

    R_in_Rsun = (R / const.R_sun).decompose().value
    L_in_Lsun = (L / const.L_sun).decompose().value

    T = 5780.0 * (L_in_Lsun / R_in_Rsun**2) ** 0.25 * u.K

    E_threshold = 13.605 * u.eV
    x_0 = (E_threshold / (const.k_B * T)).decompose().value

    if x_0 > 45.0:
        return 0.0 / u.s

    photon_integral_sum = 0.0
    for n in range(1, 6):
        exp_term = np.exp(-n * x_0)
        if exp_term < 1e-10:
            break
        photon_integral_sum += (exp_term / n) * (x_0**2 + (2.0 * x_0) / n + 2.0 / n**2)

    prefactor = 15.0 / np.pi**4
    N_dot_ion = (L / (const.k_B * T)) * prefactor * photon_integral_sum
    return N_dot_ion.to(1 / u.s)


# -----------------------------------------------------------------------------
# Read the simulated r_hii(t) from the snapshots, in both conventions (see
# the module docstring): ever-tagged gas extent and the star's own current
# HIIRegionRadii (instantaneous front).
# -----------------------------------------------------------------------------
def read_simulated_r_hii(snapshot_glob):
    if sw is None:
        raise RuntimeError("swiftsimio is required to read the snapshots.")

    files = sorted(glob.glob(snapshot_glob))
    if not files:
        raise RuntimeError(f"No snapshots found matching {snapshot_glob!r}.")

    times, r_hii, r_hii_max, r_now = [], [], [], []
    star_mass_msun = None
    n_H_atom_cc = None
    n_H_source = None
    boxsize_kpc = None

    for f in files:
        data = sw.load(f)
        if boxsize_kpc is None:
            boxsize_kpc = float(np.min(data.metadata.boxsize.to("kpc").value))
        if len(data.stars.hiiregion_radii) == 0:
            continue

        times.append(data.metadata.time.to("Myr").value)
        r_now.append(float(np.max(data.stars.hiiregion_radii).to("kpc").value))

        star_id = int(data.stars.particle_ids[0])
        ever_tagged = data.gas.hiistar_ids == star_id
        if np.any(ever_tagged):
            star_pos = data.stars.coordinates[0].to("kpc").value
            gas_pos = data.gas.coordinates[ever_tagged].to("kpc").value
            box_kpc = data.metadata.boxsize.to("kpc").value
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
                raise RuntimeError(f"Expected exactly 1 star, found {n_stars}.")
            star_mass_msun = float(np.max(data.stars.masses).to("Msun").value)
        if n_H_atom_cc is None:
            rho_g_cm3 = float(np.mean(data.gas.densities).to("g/cm**3").value)
            X_H_raw = data.metadata.parameters.get(
                "GrackleCooling:HydrogenFractionByMass"
            )
            X_H = float(X_H_raw) if X_H_raw is not None else 0.76
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

    @return (median T [K], source filename), or (None, None) if no
    snapshot holds enough tagged gas.
    """
    import h5py

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
            else:
                # No species arrays (COOLING_GRACKLE_MODE 0): this selection
                # is fully ionized by construction, so use the fully-ionized
                # primordial mu.
                mu = 0.6
            u_cgs = gas["InternalEnergies"][ionized] * u_to_cgs
            gamma_m1 = 2.0 / 3.0
            T = u_cgs * gamma_m1 * mu * const.m_p.cgs.value / const.k_B.cgs.value
            return float(np.median(T)), filename
    return None, None


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
    T_measured_K, source = measure_ionized_temperature_K(files)
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
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-s", "--snapshot-glob", default="snap/snapshot_*.hdf5")
    parser.add_argument("-o", "--output", default="starbench_check.png")
    parser.add_argument(
        "--tol",
        type=float,
        default=0.1,
        help="Relative-error tolerance for a pass (default: 0.1, "
        "tighter than the Spitzer check's 0.3 since STARBENCH "
        "agrees with high-resolution sims to <~2%%).",
    )
    parser.add_argument(
        "--T-ionized-K",
        type=float,
        default=None,
        dest="T_ionized_K",
        help="Ionized-gas temperature the reference is computed at. Default: "
        "measured from the run's own tagged gas, which always matches "
        "whatever build flag/metallicity convention the run used. An "
        "explicit value overrides the measurement but warns loudly if the "
        "two disagree by >10%%.",
    )
    parser.add_argument("--T-neutral-K", type=float, default=1e3, dest="T_neutral_K")
    args = parser.parse_args()
    T_ionized_K, T_i_marker = resolve_T_ionized_K(
        args.T_ionized_K, sorted(glob.glob(args.snapshot_glob))
    )

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
    box_half_width = 0.5 * boxsize

    Q_H = ionizing_photon_rate(star_mass_msun)
    alpha_B = alpha_b_hui_gnedin(T_ionized_K * u.K)
    R_St = ((3 * Q_H / (4 * np.pi * alpha_B * n_H**2)) ** (1 / 3.0)).to(u.pc)
    c_i = (np.sqrt(const.k_B * T_ionized_K * u.K / (0.5 * const.m_p))).to(u.km / u.s)
    c_o = (np.sqrt(const.k_B * args.T_neutral_K * u.K / const.m_p)).to(u.km / u.s)

    print(f"Star mass          : {star_mass_msun:.3f} Msun")
    print_n_H_with_source(n_H, n_H_source, files)
    print(f"Q_H                : {Q_H.to(1/u.s).value:.3e} 1/s")
    print(f"R_St               : {R_St:.4g}")
    print(f"c_i = {c_i:.4g}, c_o = {c_o:.4g}")
    print(f"Box half-width     : {box_half_width.to(u.pc):.4g}")

    t_grid = np.linspace(max(t_sim.min().value, 1e-4), t_sim.max().value, 400) * u.Myr
    R_I, R_II, R_SB = starbench_curve(
        t_grid.to(u.Myr).value,
        R_St.to(u.pc).value,
        c_i.value,
        c_o.value,
        T_i_K=T_ionized_K,
        T_o_K=args.T_neutral_K,
    )

    # Interpolate the STARBENCH curve onto the simulation's own snapshot
    # times for the comparison-point error (avoid re-solving the ODE once
    # per snapshot).
    R_SB_at_t_sim = np.interp(t_sim.to(u.Myr).value, t_grid.to(u.Myr).value, R_SB)

    r_sim_pc = r_sim.to(u.pc).value
    r_sim_max_pc = r_sim_max.to(u.pc).value
    r_now_pc = r_now.to(u.pc).value
    box_valid = R_SB_at_t_sim <= box_half_width.to(u.pc).value
    # Liveness must come from h_hii, not from the ever-tagged extent: the
    # HIIStarIDs stamp is never cleared, so r_sim keeps growing by advection
    # after the star dies and never returns to 0, which would silently
    # compare pure advection against a live-source analytic curve.
    alive = r_now_pc > 0
    ok = box_valid & alive
    if not np.any(ok):
        raise RuntimeError("No box-valid, alive snapshot to compare at.")
    last = np.where(ok)[0][-1]
    t_c = t_sim[last]
    R_SB_c = R_SB_at_t_sim[last]

    print(f"\nComparison time: t={t_c:.4g}  R_SB={R_SB_c:.4g} pc")
    if last < len(t_sim) - 1:
        reason = []
        if np.any(~box_valid):
            reason.append(
                f"reference curve leaves the box from "
                f"t={t_sim[np.where(~box_valid)[0][0]]:.4g}"
            )
        first_live = np.where(alive)[0][0]
        dead = np.where(~alive[first_live:])[0]
        if len(dead) > 0:
            reason.append(
                f"star inactive/dead from t={t_sim[first_live + dead[0]]:.4g}"
            )
        print(
            f"  (not the final simulated time t={t_sim[-1]:.4g}: "
            f"{'; '.join(reason) if reason else 'no reason recorded'})"
        )
    for label, r_c in [
        ("ever-tagged extent ", r_sim_pc[last]),
        ("instantaneous front", r_now_pc[last]),
    ]:
        rel_error = abs(r_c - R_SB_c) / R_SB_c
        verdict = "PASS" if rel_error <= args.tol else "FAIL"
        print(
            f"  {label}: r_sim={r_c:.4g} pc  "
            f"rel_error={rel_error:.2%}  [{verdict}{T_i_marker}]"
        )
    print(
        f"  (secondary diagnostic, true max ever-tagged extent: "
        f"{r_sim_max_pc[last]:.4g} pc)"
    )
    if abs(r_sim_pc[last] - r_now_pc[last]) / max(r_now_pc[last], 1e-10) > 0.1:
        print(
            "  NOTE: the two conventions diverge by >10% -- the region "
            "receded. The scheme-appropriate line is the meaningful one: "
            "ever-tagged for flag-and-floor, instantaneous for rate-coupled "
            "(HII_couple_ionization_rate) runs. See the module docstring."
        )

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(
        t_sim.to(u.Myr),
        r_sim_pc,
        "o-",
        color="#d62728",
        markersize=4,
        label="Simulation (ever-tagged $r_{\\rm HII}$)",
    )
    ax.plot(
        t_sim.to(u.Myr),
        r_now_pc,
        "s-",
        color="#1f77b4",
        markersize=3,
        alpha=0.8,
        label="Simulation (instantaneous $h_{\\rm HII}$)",
    )
    ax.plot(t_grid, R_I, ":", color="grey", label="Raga-I (Eqn 8)")
    ax.plot(t_grid, R_II, "-.", color="grey", label="Raga-II (Eqn 11)")
    ax.plot(t_grid, R_SB, "--", color="black", label="STARBENCH (Eqns 28-29)")
    ax.axhline(
        box_half_width.to(u.pc).value,
        color="firebrick",
        linestyle=":",
        alpha=0.6,
        label="Box half-width",
    )
    ax.axvline(
        t_c.to(u.Myr).value,
        color="black",
        linestyle="-.",
        alpha=0.4,
        label="Comparison point",
    )
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel(r"HII region radius $r_{\rm HII}$ [pc]")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
