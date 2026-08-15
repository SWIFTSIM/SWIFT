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
Cosmological counterpart of ../StromgrenSphere/stromgren_analytic_check.py:
compares the simulated HII region growth against the classical Stromgren
equilibrium radius + Spitzer (1978) D-type expansion, evaluated at the
INSTANTANEOUS PHYSICAL density.

Every length/density SWIFT writes to a cosmological snapshot is COMOVING.
This script converts to physical explicitly, per field, using that field's
own "a-scale exponent" HDF5 attribute and the snapshot's Header/Scale-factor
(physical = comoving * a**exponent) -- it does not assume a fixed convention
per field name. That conversion is the entire point of this script: for
params_identity.yml (a~=1 throughout) it should barely matter, but for
params_highz.yml (a_begin=0.1) skipping it would be a ~1000x density error.

Usage:
    python3 cosmo_stromgren_analytic_check.py [-s snap/snapshot] [-o out.png]
"""

import argparse
import glob
import os
import re

import h5py
import numpy as np
from astropy import constants as const
from astropy import units as u

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def a_scale_exponent(dataset):
    """The "a-scale exponent" n such that physical = comoving * a**n for
    this HDF5 dataset (written by SWIFT on every field it outputs)."""
    return float(dataset.attrs["a-scale exponent"][0])


# -----------------------------------------------------------------------------
# Q_H(mass): direct port of radiation_get_individual_star_ionizing_photon_
# emission_rate_fit() and the radius/luminosity/temperature fits it calls
# (src/feedback/GEAR/radiation.c). All in CGS/astropy units here; the C code
# works in SWIFT internal units but the physics is unit-system independent.
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
    """Q_H, the ionizing photon emission rate for a single star of the given
    mass, from a blackbody fit to its (R, L) -- direct port of
    radiation_get_individual_star_ionizing_photon_emission_rate_fit()."""
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
# Analytic Stromgren-sphere solution, evaluated at the run's own applied
# T_ionized (see ../StromgrenSphere/stromgren_analytic_check.py for why this
# departs from the classical flat 1e4 K assumption).
# -----------------------------------------------------------------------------
def alpha_b_hui_gnedin(T):
    """Case-B recombination coefficient (Hui & Gnedin 1997, MNRAS 292, 27,
    Appendix A), the same fit used in
    radiation_get_case_b_recombination_coefficient_cgs
    (src/feedback/GEAR/radiation.c)."""
    lam = 315614.0 / T.to(u.K).value
    return (
        2.753e-14 * lam**1.5 / (1.0 + (lam / 2.740) ** 0.407) ** 2.242 * u.cm**3 / u.s
    )


def sound_speed_ionized(T_ionized):
    """Isothermal sound speed of photoionized hydrogen gas, scaled from the
    standard 1e4 K reference value (12.85 km/s) to the actual applied
    ionized temperature (c_s ~ sqrt(T) at fixed composition)."""
    return 12.85 * u.km / u.s * np.sqrt(T_ionized / (1e4 * u.K))


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
# Read the simulated, PHYSICAL r_hii(t) from the snapshots
# -----------------------------------------------------------------------------
def read_simulated_r_hii(snapshot_glob):
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
    unit_time_to_Myr = None
    time_begin_internal = None

    for f in files:
        with h5py.File(f, "r") as h:
            header = h["Header"]
            a = float(header.attrs["Scale-factor"][0])
            u_L = h["Units"].attrs["Unit length in cgs (U_L)"][0]
            u_M = h["Units"].attrs["Unit mass in cgs (U_M)"][0]
            u_t = h["Units"].attrs["Unit time in cgs (U_t)"][0]
            if unit_time_to_Myr is None:
                unit_time_to_Myr = (u_t * u.s).to(u.Myr).value

            gas = h["PartType0"]
            stars = h["PartType4"]

            # Header/Time is the ABSOLUTE physical cosmic age at this
            # snapshot's scale factor (src/engine.c: e->time_begin =
            # cosmology->time_begin, itself the age of the Universe at
            # a_begin via the Friedmann integral) -- not elapsed run time,
            # and not a comoving quantity needing an a-factor. Subtract the
            # first snapshot's cosmic age to get elapsed time since a_begin,
            # matching the non-cosmological check's t=0-at-start convention.
            if time_begin_internal is None:
                time_begin_internal = float(header.attrs["Time"][0])
            times.append(
                (float(header.attrs["Time"][0]) - time_begin_internal)
                * unit_time_to_Myr
            )

            boxsize_comoving = header.attrs["BoxSize"][0]  # cubic box
            boxsize_physical_kpc = (
                boxsize_comoving * a ** a_scale_exponent(gas["Coordinates"])
            ) * (u_L * u.cm).to(u.kpc).value

            # r_hii from the star's own HIIRegionRadii (h_hii, comoving) is
            # frozen at the moment each gas particle was tagged and never
            # updates as that particle keeps moving -- see the non-cosmological
            # check for the full rationale. Use the current PHYSICAL position
            # of every gas particle this star has *ever* tagged instead
            # (HIIStarIDs, dimensionless, never reset when a tag expires).
            # h_hii, in physical units. Unlike the ever-tagged extent this
            # returns to exactly 0 when the star dies or ages past
            # HII_max_age, so it, not r_hii, is what tells the comparison
            # window when the source went dark.
            r_now.append(
                float(
                    np.max(
                        stars["HIIRegionRadii"][...]
                        * a ** a_scale_exponent(stars["HIIRegionRadii"])
                    )
                )
                * (u_L * u.cm).to(u.kpc).value
            )

            star_id = int(stars["ParticleIDs"][0])
            hii_star_ids = gas["HIIStarIDs"][...]
            ever_tagged = hii_star_ids == star_id
            if np.any(ever_tagged):
                star_pos_comoving = stars["Coordinates"][0]
                gas_pos_comoving = gas["Coordinates"][...][ever_tagged]
                a_L = a_scale_exponent(gas["Coordinates"])
                star_pos_kpc = (star_pos_comoving * a**a_L) * (u_L * u.cm).to(
                    u.kpc
                ).value
                gas_pos_kpc = (gas_pos_comoving * a**a_L) * (u_L * u.cm).to(u.kpc).value
                # Minimum-image convention in PHYSICAL coordinates, using the
                # PHYSICAL box size at this snapshot's scale factor.
                dx = gas_pos_kpc - star_pos_kpc
                dx -= boxsize_physical_kpc * np.round(dx / boxsize_physical_kpc)
                r_robust, r_max = robust_ever_tagged_radius(dx)
                r_hii.append(r_robust)
                r_hii_max.append(r_max)
            else:
                r_hii.append(r_now[-1])
                r_hii_max.append(r_now[-1])

            if star_mass_msun is None:
                masses = stars["Masses"][...]
                if len(masses) != 1:
                    raise RuntimeError(
                        f"Expected exactly 1 star (this check assumes "
                        f"star_type=single_star, matching the Q_H(mass) fit "
                        f"reimplemented below), found {len(masses)}."
                    )
                star_mass_msun = float(
                    (masses[0] * a ** a_scale_exponent(stars["Masses"]))
                    * (u_M * u.g).to(u.Msun).value
                )

            if n_H_atom_cc is None:
                # Hydrogen number density from the gas PHYSICAL density and
                # the Grackle-reported hydrogen mass fraction (see the
                # non-cosmological check for why this is read, not assumed).
                rho_comoving = gas["Densities"][...]
                a_rho = a_scale_exponent(gas["Densities"])
                rho_g_cm3 = (
                    float(np.mean(rho_comoving * a**a_rho))
                    * ((u_M * u.g) / (u_L * u.cm) ** 3).to(u.g / u.cm**3).value
                )
                X_H_raw = h["Parameters"].attrs.get(
                    "GrackleCooling:HydrogenFractionByMass"
                )
                # Parameters-group attributes are scalar byte-strings (the
                # parsed parameter text), not numeric arrays -- indexing
                # with [0] would silently take the first ASCII byte instead
                # of decoding the value.
                X_H = float(X_H_raw.decode()) if X_H_raw is not None else 0.716
                n_H_atom_cc = (rho_g_cm3 * u.g / u.cm**3 * X_H / const.m_p).to(
                    1 / u.cm**3
                )
                n_H_source = f

    return (
        u.Quantity(times, u.Myr),
        u.Quantity(r_hii, u.kpc),
        u.Quantity(r_hii_max, u.kpc),
        u.Quantity(r_now, u.kpc),
        star_mass_msun,
        n_H_atom_cc,
        n_H_source,
        files,
    )


def measure_ionized_temperature_K(files, min_count=20):
    """Median PHYSICAL temperature of the currently-tagged (ionized) gas,
    from the latest snapshot holding at least min_count tagged particles.

    @return (median T [K], source filename), or (None, None) if no
    snapshot holds enough tagged gas.
    """
    for filename in reversed(files):
        with h5py.File(filename, "r") as h:
            gas = h["PartType0"]
            ionized = gas["IsIonizedFlags"][...].astype(bool)
            if np.sum(ionized) < min_count:
                continue
            a = float(h["Header"].attrs["Scale-factor"][0])
            u_L = h["Units"].attrs["Unit length in cgs (U_L)"][0]
            u_t = h["Units"].attrs["Unit time in cgs (U_t)"][0]
            u_to_cgs_comoving = (
                u_L / u_t
            ) ** 2  # comoving internal energy unit -> erg/g
            a_u = a_scale_exponent(gas["InternalEnergies"])
            if "HI" in gas:
                inv_mu = (
                    gas["HI"][...][ionized]
                    + 2.0 * gas["HII"][...][ionized]
                    + 0.25 * gas["HeI"][...][ionized]
                    + 0.5 * gas["HeII"][...][ionized]
                    + 0.75 * gas["HeIII"][...][ionized]
                )
                mu = 1.0 / np.clip(inv_mu, 1e-10, None)
            else:
                # No species arrays (COOLING_GRACKLE_MODE 0): this selection
                # is fully ionized by construction, so use the fully-ionized
                # primordial mu.
                mu = 0.6
            u_cgs_physical = (
                gas["InternalEnergies"][...][ionized] * a**a_u * u_to_cgs_comoving
            )
            gamma_m1 = 2.0 / 3.0
            T = (
                u_cgs_physical
                * gamma_m1
                * mu
                * const.m_p.cgs.value
                / const.k_B.cgs.value
            )
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
        default="cosmo_stromgren_analytic_check.png",
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
        "Default: measured from the run's own tagged gas.",
    )
    args = parser.parse_args()
    files = sorted(glob.glob(args.snapshot_glob))
    T_ionized_K, T_i_marker = resolve_T_ionized_K(args.T_ionized_K, files)

    T_IONIZED_K = T_ionized_K * u.K
    ALPHA_B = alpha_b_hui_gnedin(T_IONIZED_K)
    C_S_IONIZED = sound_speed_ionized(T_IONIZED_K)

    (
        t_sim,
        r_sim,
        r_sim_max,
        r_now,
        star_mass_msun,
        n_H,
        n_H_source,
        files,
    ) = read_simulated_r_hii(args.snapshot_glob)

    if star_mass_msun is None or n_H is None:
        raise RuntimeError("Could not infer star mass / n_H from the snapshots.")

    Q_H = ionizing_photon_rate(star_mass_msun)
    R_st = stromgren_radius(Q_H, n_H, ALPHA_B)

    with h5py.File(files[0], "r") as h:
        a0 = float(h["Header"].attrs["Scale-factor"][0])
        u_L = h["Units"].attrs["Unit length in cgs (U_L)"][0]
        boxsize_comoving0 = h["Header"].attrs["BoxSize"][0]
        a_L = a_scale_exponent(h["PartType0/Coordinates"])
        box_half_width0 = (
            0.5 * boxsize_comoving0 * a0**a_L * (u_L * u.cm).to(u.kpc).value
        ) * u.kpc

    print(f"Star mass          : {star_mass_msun:.3f} Msun")
    print_n_H_with_source(n_H, n_H_source, files, label="n_H (physical, t=0)")
    print(f"Q_H                : {Q_H:.4g}")
    print(
        f"T_ionized (applied): {T_IONIZED_K:.4g}  ->  alpha_B = {ALPHA_B:.4g}, "
        f"c_s = {C_S_IONIZED:.4g}"
    )
    print(f"Equilibrium R_st   : {R_st.to(u.pc):.4g} = {R_st.to(u.kpc):.4g}")
    print(f"Physical box half-width (t=0): {box_half_width0.to(u.pc):.4g}")
    if R_st > box_half_width0:
        print(
            "  WARNING: R_st exceeds the physical box half-width -- the "
            "periodic box is too small to contain a full infinite-medium "
            "Stromgren sphere. The analytic curve below is an upper bound."
        )

    r_analytic = dtype_expansion_radius(t_sim, R_st, C_S_IONIZED).to(u.kpc)

    # Same validity window logic as the non-cosmological check: only compare
    # while the analytic radius still fits the (physical) box and the star
    # is alive/tracked (r_sim > 0). box_half_width0 is used throughout since
    # its scale-factor drift over these runs' short a-span is negligible
    # (<2.4% even for params_highz.yml) compared to this check's tolerance.
    # Liveness comes from h_hii (r_now), not from the ever-tagged extent:
    # the HIIStarIDs stamp is never cleared, so r_sim keeps growing by
    # advection after the star dies and never returns to 0.
    valid = (r_analytic <= box_half_width0) & (r_now > 0)
    box_exceeded_at = t_sim[r_analytic > box_half_width0]
    nonzero_idx = np.where(r_now > 0)[0]
    dead_mask = np.zeros(len(r_now), dtype=bool)
    if len(nonzero_idx) > 0:
        dead_mask[nonzero_idx[0] + 1 :] = r_now[nonzero_idx[0] + 1 :] == 0
    star_dead_at = t_sim[dead_mask]
    if not np.any(valid):
        raise RuntimeError(
            "No snapshot has both an analytic radius within the physical "
            "box and an alive, actively-tracked r_sim > 0 -- no meaningful "
            "comparison point exists in this run."
        )
    last_valid = np.where(valid)[0][-1]
    t_sim_final = t_sim[last_valid]
    r_sim_final = r_sim[last_valid]
    r_analytic_final = r_analytic[last_valid]
    rel_error = float(
        (np.abs(r_sim_final - r_analytic_final) / r_analytic_final).decompose()
    )
    verdict = "PASS" if rel_error <= args.tol else "FAIL"
    print(
        f"Comparison time: t={t_sim_final:.4g}  r_sim={r_sim_final:.4g}  "
        f"r_analytic={r_analytic_final:.4g}  rel_error={rel_error:.2%}  "
        f"[{verdict}{T_i_marker}]"
    )
    print(
        f"  (secondary diagnostic, true max ever-tagged extent: "
        f"{r_sim_max[last_valid]:.4g})"
    )
    if last_valid < len(t_sim) - 1:
        reason = []
        if len(box_exceeded_at) > 0:
            reason.append(f"box half-width exceeded from t={box_exceeded_at[0]:.4g}")
        if len(star_dead_at) > 0:
            reason.append(f"star inactive/dead from t={star_dead_at[0]:.4g}")
        print(
            f"  (not the final simulated time t={t_sim[-1]:.4g} -- "
            f"{'; '.join(reason)}.)"
        )

    fig, ax = plt.subplots()
    ax.plot(
        t_sim.to(u.Myr), r_sim.to(u.pc), "o-", label="Simulation ($r_{hii}$, physical)"
    )
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
        box_half_width0.to(u.pc).value,
        color="firebrick",
        linestyle=":",
        label="Physical box half-width",
    )
    ax.axvline(
        t_sim_final.to(u.Myr).value,
        color="black",
        linestyle="-.",
        alpha=0.6,
        label="Comparison point",
    )
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("HII region radius [pc] (physical)")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
