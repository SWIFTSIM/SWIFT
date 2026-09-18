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
Checks Grackle's photoelectric heating response to the ISRF field this
codebase injects, against a `with_photoelectric_heating: 1` run compared
to an otherwise-identical `with_photoelectric_heating: 0` run (same IC,
same density/metallicity, same `time_end`).

Sign/magnitude check (Tier 1, original)
----------------------------------------
Bins gas by radius from the star at each run's LAST snapshot and checks,
in every radial bin where the injected field is non-negligible, that (1)
the heating-on run is hotter than the heating-off run (sign check) and (2)
the heating-on temperature stays under a physically-motivated ceiling
rather than showing signs of a numerical blow-up (magnitude check).

Quantitative rate check (Tier 2)
---------------------------------
Measures the actual photoelectric heating RATE and compares it against
Wolfire et al. (1995) Eq. 1 with the metallicity factor Grackle's
`photoelectric_heating_efficiency: constant` (SWIFT's default, igammah=2)
path actually implements,

    Gamma_PE = gamma_ha * epsilon * G0 * n_H * Z' ,

with ``gamma_ha = 1e-24 erg/s``, ``epsilon = 0.05`` (Grackle's own
hardcoded constants for igammah=2, `rate_functions.c`'s `gammah_rate()`),
``G0`` the Habing-unit ISRF strength, ``n_H`` the TOTAL hydrogen number
density and ``Z' = Z/Zsun``. This is the Z'-corrected form, not the
commonly-cited metallicity-blind one: Grackle's `cool1d_multi_g.F` (its
`igammah=2` branch, `edot += gammaha_eff * rhoH * dom_inv * dust2gas /
fgr`) folds in `dust2gas = fgr * metallicity` in the default
(`idustfield=0`, confirmed for this codebase: it never sets a dust
density field) configuration, so `fgr` cancels and the metallicity factor
survives. Verified to agree with a compiled `libgrackle` call-through to
better than 1e-6 relative error at Z' in {1, 0.1, 3}:
`theory/GEAR/Radiation/verify_photoelectric_heating_rate.py`, whose
constants this script reuses rather than retyping.

For `ispecies=0` (COOLING_GRACKLE_MODE 0, the tabulated path this
fixture's `grackle_0` build uses), Grackle forms `n_H` as
``HydrogenFractionByMass * (rho - rho*Z) / m_H`` (`cool1d_multi_g.F`,
the `imetal=1` branch of the `ispecies=0` temperature block), not the
naive `HydrogenFractionByMass * rho / m_H`; this script reproduces that
exact expression. Note `n_H / rho` is then a pure constant
(`HydrogenFractionByMass * (1-Z) / m_H`, independent of `rho`), so
`Gamma_PE / rho` depends only on `G0` and `Z'`, not on the gas density at
all -- relevant to the error budget below, since it means a particle's
own density does not need to be constant across the measurement window
for the comparison to be valid.

The equilibrium temperature Tier 1 checks mixes Gamma_PE with every other
cooling/heating term Grackle applies, so it cannot isolate the rate. This
check instead measures ``du/dt`` as a plain forward difference over one
inter-snapshot interval, matched particle-by-particle by `ParticleIDs`,
while the on/off gas states still coincide (checked below, not assumed):
subtracting the two runs' `du/dt` cancels every process identical between
them (cooling, any work term), isolating the ON/OFF difference, which is
Gamma_PE/rho to leading order. `G0`, `n_H` and `Z'` are read from the
heating-ON run's own snapshot at the start of that interval (the
heating-off run's ISRF fields are zero: `with_photoelectric_heating`
gates the whole ISRF injection, not just its consumption by Grackle,
`radiation_iact.h`) -- an input to this test, not one of its claims,
matching `isrf_h2_photodissociation_check.py`'s own convention of taking
`u` from the snapshot rather than predicting it from a transport
solution.

Why this needs its own, separately-cadenced run pair
------------------------------------------------------
`ISRF_propagation: 0` stores an "instantaneous field", not a time-averaged
one (`radiation_iact.h`): `u` is reset to 0 on a gas particle's first
touch by any star THIS GLOBAL STEP and accumulates only touches within
that same step, so what a snapshot records is essentially the star's MOST
RECENT touch's own deposit, `Delta_t_star * L * weight * extinction / m`,
not a quantity that has converged to a stable field as the timestep
refines uniformly. Measured directly on this fixture (two short runs,
`dt_max` 1e-7 vs 1e-8 internal, `G0` compared at the same physical time):
the coarser run reads `G0=0` where the finer one reads `G0~1.4e3`, and at
the fixture's own default cadence (`dt_max=1e-6`, `delta_time=1e-5`) the
near-star gas has already been heated from its 1000 K initial condition
to 11,000-20,000 K by the FIRST snapshot with any nonzero field at all --
Gamma_PE*dt/u there is 1e2-1e4, i.e. the "small perturbation" this
method needs does not hold anywhere in that run. Neither is a bug in this
script; both are consequences of `ISRF_propagation: 0`'s reset-on-touch
design measured at a cadence too coarse for it, and are reported here as
a real characteristic of that mode rather than worked around silently.

The fix is cadence, not a different formula: at `delta_time=1e-9`,
`dt_max=1e-10` internal units (`time_end=2e-9`, two snapshots), the same
measurement gives `T` still within a few K of the 1000 K initial
condition, `Gamma_PE*dt/u` of order 1e-4 (median, 42 illuminated
sub-cutoff particles) and `G0` unchanged between the two snapshots to
float32 precision -- a genuinely linear, cadence-converged regime. This
is a SEPARATE run pair from Tier 1's equilibrium `heating_on`/
`heating_off` (`--rate-on-snapshot`/`--rate-off-snapshot`; Tier 2 is
skipped, not failed, when these are not given -- see README), since
Tier 1 needs a run long enough to reach equilibrium and Tier 2 needs one
short enough to still be in the linear regime; no single run pair serves
both.

Pass bar
--------
``--rate-tol`` (default 0.12) is a FIXED number, derived once from the
fine-cadence reference run above and frozen, not recomputed from whatever
run is being checked (a self-referential bar can drift with a regression
it should catch). At that cadence, measured on the 42 illuminated,
`T < 2e4 K` particles:

* `G0` drift between the two snapshots: exactly 0 (float32-identical) --
  the window-averaging/truncation term this method would otherwise carry
  is negligible at this cadence.
* `Gamma_PE * dt / u` (linearity): median 3e-4, 16th-84th percentile
  4e-5 to 2.5e-3 -- comfortably inside the perturbative regime the on/off
  cancellation needs; also a hard precondition gate (`--perturbation-tol`,
  default 0.05), since outside this regime the whole method is invalid,
  not merely imprecise (see "Why this needs its own..." above).
* float32 snapshot storage: ~2e-7 relative, negligible.
* Measured vs predicted relative error, the residual after the above:
  median 7.3 percent, 16th-84th percentile 6.1 to 15.9 percent. With `G0`
  drift and the linearity term both negligible at this cadence, this
  residual is per-particle sampling noise in the injection sum (each
  particle's own SPH kernel weight to the star), not finite-difference
  truncation.

``--rate-tol 0.12`` sits above the measured median with headroom for
run-to-run variation, without being loose enough to miss a real
factor-of-few formula bug (e.g. a dropped `Z'` factor, or a wrong `n_H`
definition, both of which would show as an order-of-magnitude offset, not
a 10-15 percent one).
"""

import argparse
import glob
import sys
from typing import Dict, Optional

import h5py
import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

K_B_CGS = 1.380649e-16
M_P_CGS = 1.67262192e-24
GAMMA_M1 = 2.0 / 3.0  # monatomic ideal gas, gamma = 5/3
C_LIGHT_CGS = 2.99792458e10

# src/feedback/GEAR/radiation.h
RADIATION_HABING_FLUX_CGS = 1.6e-3
RADIATION_GRACKLE_SOLAR_METAL_FRACTION = 0.01295

# Grackle's own hardcoded constants for photoelectric_heating_efficiency:
# constant (igammah=2), src/clib/rate_functions.c's gammah_rate(); these
# are Grackle's compiled behaviour being checked, not values this codebase
# owns, so they are not read from radiation.h.
GRACKLE_GAMMA_HA_CGS = 1.0e-24  # erg/s
GRACKLE_EPSILON = 0.05

# Neither run is expected to sustain a temperature above this: known
# photoelectric-/collisional-ionization-equilibrium temperatures for
# irradiated atomic gas top out around 1e4-1e5 K (Lyman-alpha and
# collisional-ionization cooling both strengthen sharply above ~1e4 K).
# 1e6 K would require a channel this setup does not have and is a much
# more plausible signature of a numerical blow-up than of real heating.
IMPLAUSIBLE_TEMPERATURE_K = 1e6

# igammah=2 is pinned to zero above this temperature (cool1d_multi_g.F);
# the rate check restricts itself to particles below it.
GAMMAH_TEMPERATURE_CEILING_K = 2.0e4


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--on-snapshot",
        default="snap_on/snapshot_*.hdf5",
        help="Glob pattern for the with_photoelectric_heating=1 run's "
        "snapshots (default: %(default)s)",
    )
    parser.add_argument(
        "--off-snapshot",
        default="snap_off/snapshot_*.hdf5",
        help="Glob pattern for the with_photoelectric_heating=0 run's "
        "snapshots (default: %(default)s)",
    )
    parser.add_argument("--n-bins", type=int, default=12, help="Number of radial bins.")
    parser.add_argument(
        "--field-threshold",
        type=float,
        default=1e-3,
        help="A radial bin is considered field-illuminated if the on-run's "
        "mean (u_FUV+u_LW) there exceeds this fraction of the innermost "
        "bin's value (default: %(default)s).",
    )
    parser.add_argument(
        "--output",
        default="isrf_photoelectric_heating_check.png",
        help="Output plot filename.",
    )
    parser.add_argument(
        "--rate-on-snapshot",
        default=None,
        help="Glob pattern for the FINE-CADENCE with_photoelectric_heating=1 "
        "run's snapshots, for Tier 2 (the quantitative rate check). A "
        "SEPARATE run from --on-snapshot: see docstring for why. Tier 2 is "
        "skipped, not failed, if this (and --rate-off-snapshot) are not "
        "given.",
    )
    parser.add_argument(
        "--rate-off-snapshot",
        default=None,
        help="Glob pattern for the FINE-CADENCE with_photoelectric_heating=0 "
        "run's snapshots, for Tier 2. See --rate-on-snapshot.",
    )
    parser.add_argument(
        "--rate-tol",
        type=float,
        default=0.12,
        help="Max median relative error on the measured photoelectric "
        "heating rate vs Wolfire et al. (1995) Eq. 1, over illuminated, "
        "T < 2e4 K particles (default: %(default)s; see docstring for the "
        "error budget this is derived from).",
    )
    parser.add_argument(
        "--perturbation-tol",
        type=float,
        default=0.05,
        help="Max median Gamma_PE*dt/u over the rate window: the on/off "
        "cancellation this check relies on is only valid in this linear "
        "regime, so this is a hard precondition, not a precision knob "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--hydrogen-fraction",
        type=float,
        default=0.76,
        help="GrackleCooling:HydrogenFractionByMass the run used "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--precondition-tol",
        type=float,
        default=1e-4,
        help="Max relative difference allowed between the on- and off-run "
        "density and temperature at the first snapshot, for the rate "
        "check's on/off subtraction to be valid (default: %(default)s; "
        "the reference run's own threaded SPH density sum differs "
        "run-to-run at the few-1e-6 level even from a bit-identical IC, "
        "so this is set well above that floor).",
    )
    return parser.parse_args()


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
        unit_time_cgs = float(np.asarray(units.attrs["Unit time in cgs (U_t)"]).flat[0])
        u_to_cgs = (unit_length_cgs / unit_time_cgs) ** 2
        density_to_cgs = unit_mass_cgs / unit_length_cgs**3

        gas = f["/PartType0"]
        ids = gas["ParticleIDs"][:]
        pos = gas["Coordinates"][:, :]
        u = gas["InternalEnergies"][:] * u_to_cgs
        # Raw (internal-unit) values, used only for the Tier 1 radial-bin
        # threshold; the Tier 2 rate check uses the cgs versions below.
        u_pe = gas["FUVSpecificEnergies"][:]
        u_lw = gas["LWSpecificEnergies"][:]
        u_pe_cgs = u_pe * u_to_cgs
        u_lw_cgs = u_lw * u_to_cgs
        density_cgs = gas["Densities"][:] * density_to_cgs

        metals = gas["MetalMassFractions"][:]
        named_columns = "/SubgridScheme/NamedColumns/MetalMassFractions"
        if metals.ndim == 2:
            metal_names = [name.decode() for name in f[named_columns][:]]
            metallicity = metals[:, metal_names.index("Metals")]
        else:
            metallicity = metals

        # Species-based mean molecular weight when available
        # (COOLING_GRACKLE_MODE >= 1); otherwise fall back to the
        # primordial-neutral value (X=0.76, Y=0.24, 1/mu = X + Y/4).
        if "HI" in gas:
            inv_mu = (
                gas["HI"][:]
                + 2.0 * gas["HII"][:]
                + 0.25 * gas["HeI"][:]
                + 0.5 * gas["HeII"][:]
                + 0.75 * gas["HeIII"][:]
            )
            mu = 1.0 / np.clip(inv_mu, 1e-10, None)
        else:
            mu = np.full(pos.shape[0], 1.0 / 0.82)

        temperature = u * GAMMA_M1 * mu * M_P_CGS / K_B_CGS

        star = f["/PartType4"]
        star_pos = star["Coordinates"][0, :]

    return dict(
        time=time,
        time_cgs=time * unit_time_cgs,
        boxsize=boxsize,
        ids=ids,
        pos=pos,
        u=u,
        u_pe=u_pe,
        u_lw=u_lw,
        u_pe_cgs=u_pe_cgs,
        u_lw_cgs=u_lw_cgs,
        density_cgs=density_cgs,
        metallicity=metallicity,
        temperature=temperature,
        star_pos=star_pos,
    )


def align_by_id(a: Dict[str, np.ndarray], b: Dict[str, np.ndarray]) -> None:
    """Sort two snapshot dicts in place by ``ParticleIDs`` and check they match.

    Parameters
    ----------
    a, b : dict of str to numpy.ndarray
        Snapshot dicts from :func:`load_snapshot`, modified in place so
        every array-valued entry is reordered to ascending particle ID.

    Raises
    ------
    RuntimeError
        If the two snapshots do not contain the same particle set.
    """
    for snapshot in (a, b):
        order = np.argsort(snapshot["ids"])
        for key, value in snapshot.items():
            if (
                isinstance(value, np.ndarray)
                and value.shape[:1] == snapshot["ids"].shape
            ):
                snapshot[key] = value[order]
    if not np.array_equal(a["ids"], b["ids"]):
        raise RuntimeError("Snapshots do not contain the same particle set")


def gamma_pe_over_rho_predicted(
    density_cgs: np.ndarray,
    u_pe_cgs: np.ndarray,
    u_lw_cgs: np.ndarray,
    metallicity: np.ndarray,
    hydrogen_fraction: float,
) -> np.ndarray:
    """Wolfire et al. (1995) Eq. 1, Z'-corrected, per unit mass.

    Reproduces Grackle's `igammah=2` path exactly (module docstring): `n_H`
    is `HydrogenFractionByMass * (rho - rho*Z) / m_H`, not the naive
    `HydrogenFractionByMass * rho / m_H`.

    Parameters
    ----------
    density_cgs : numpy.ndarray
        Gas mass density, g cm^-3.
    u_pe_cgs, u_lw_cgs : numpy.ndarray
        FUV- and LW-band specific energy, erg g^-1.
    metallicity : numpy.ndarray
        Metal mass fraction (dimensionless).
    hydrogen_fraction : float
        GrackleCooling:HydrogenFractionByMass.

    Returns
    -------
    numpy.ndarray
        Gamma_PE / rho, erg g^-1 s^-1.
    """
    g0 = C_LIGHT_CGS * density_cgs * (u_pe_cgs + u_lw_cgs) / RADIATION_HABING_FLUX_CGS
    n_h_cgs = hydrogen_fraction * density_cgs * (1.0 - metallicity) / M_P_CGS
    z_prime = metallicity / RADIATION_GRACKLE_SOLAR_METAL_FRACTION
    return GRACKLE_GAMMA_HA_CGS * GRACKLE_EPSILON * g0 * n_h_cgs * z_prime / density_cgs


def check_gate(name: str, value: float, bad: bool, message: str) -> Optional[str]:
    """Build a failure string for one pass/fail gate, failing closed on NaN/inf.

    Mirrors `isrf_h2_photodissociation_check.py`'s helper of the same name:
    a bare ``value > bar`` comparison lets a non-finite ``value`` slip
    through silently (NaN compares False against every bound), so
    finiteness is checked first.

    Parameters
    ----------
    name : str
        Short name of the gated quantity, used only in the non-finite
        message.
    value : float
        The measured value being gated.
    bad : bool
        Whether the value fails its bound; only consulted when ``value``
        is finite.
    message : str
        Failure message to use when ``value`` is finite and ``bad``.

    Returns
    -------
    str or None
        A failure message, or ``None`` if the gate passes.
    """
    if not np.isfinite(value):
        return f"{name} is non-finite ({value!r})"
    if bad:
        return message
    return None


def radial_distance(pos, star_pos, boxsize):
    """Minimum-image radial distance from the star, for a periodic box."""
    dx = pos - star_pos
    dx -= boxsize * np.round(dx / boxsize)
    return np.sqrt(np.sum(dx**2, axis=1))


def bin_by_radius(r, values, edges):
    n_bins = len(edges) - 1
    binned = np.full(n_bins, np.nan)
    for i in range(n_bins):
        sel = (r >= edges[i]) & (r < edges[i + 1])
        if sel.sum() > 0:
            binned[i] = np.mean(values[sel])
    return binned


def main():
    opt = parse_options()

    on_files = sorted(glob.glob(opt.on_snapshot))
    off_files = sorted(glob.glob(opt.off_snapshot))
    if not on_files:
        raise RuntimeError(f"No snapshots match {opt.on_snapshot!r}")
    if not off_files:
        raise RuntimeError(f"No snapshots match {opt.off_snapshot!r}")

    # Use the last snapshot of each run: both are designed to reach thermal
    # equilibrium well before time_end (see README).
    on = load_snapshot(on_files[-1])
    off = load_snapshot(off_files[-1])

    r_on = radial_distance(on["pos"], on["star_pos"], on["boxsize"])
    r_off = radial_distance(off["pos"], off["star_pos"], off["boxsize"])

    # Common radial binning, out to the extent of whichever run's
    # illuminated gas has spread furthest (the heating-on run's near-star
    # gas mildly expands as it heats: see README).
    r_max = max(r_on.max(), r_off.max()) * 0.5
    edges = np.linspace(0.0, r_max, opt.n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])

    field_on = bin_by_radius(r_on, on["u_pe"] + on["u_lw"], edges)
    T_on = bin_by_radius(r_on, on["temperature"], edges)
    T_off = bin_by_radius(r_off, off["temperature"], edges)

    # A bin counts as "illuminated" if the on-run's own field there is a
    # non-negligible fraction of the peak (innermost, non-nan) bin value --
    # avoids comparing temperatures out where neither run's gas is affected
    # by the star at all (both would trivially sit at the same floor).
    valid = ~np.isnan(field_on) & ~np.isnan(T_on) & ~np.isnan(T_off)
    peak_field = np.nanmax(field_on[valid]) if valid.any() else 0.0
    illuminated = valid & (field_on > opt.field_threshold * peak_field)

    print(f"On-run snapshot: {on_files[-1]} (t={on['time']:.4e})")
    print(f"Off-run snapshot: {off_files[-1]} (t={off['time']:.4e})")
    print(
        f"{illuminated.sum()}/{opt.n_bins} radial bins are illuminated "
        f"(field > {opt.field_threshold:.1e} x peak)."
    )

    ok_sign = True
    ok_magnitude = True
    for i in np.where(illuminated)[0]:
        sign_ok = T_on[i] > T_off[i]
        magnitude_ok = T_on[i] < IMPLAUSIBLE_TEMPERATURE_K
        ok_sign &= sign_ok
        ok_magnitude &= magnitude_ok
        print(
            f"  bin {i}: r={centres[i]:.4e}, T_on={T_on[i]:.4e} K, "
            f"T_off={T_off[i]:.4e} K, "
            f"{'PASS' if sign_ok else 'FAIL'} (sign), "
            f"{'PASS' if magnitude_ok else 'FAIL'} (magnitude)"
        )

    if illuminated.sum() == 0:
        print("No illuminated bins found -- cannot check anything.")
        ok_sign = ok_magnitude = False

    print(
        f"Sign check (T_on > T_off in every illuminated bin): "
        f"{'PASS' if ok_sign else 'FAIL'}"
    )
    print(
        f"Magnitude check (T_on < {IMPLAUSIBLE_TEMPERATURE_K:.1e} K in "
        f"every illuminated bin): {'PASS' if ok_magnitude else 'FAIL'}"
    )

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.semilogy(centres[valid], T_on[valid], "o-", label="heating on")
    ax.semilogy(centres[valid], T_off[valid], "s-", label="heating off")
    ax.axvspan(
        centres[illuminated].min() if illuminated.any() else 0,
        centres[illuminated].max() if illuminated.any() else 0,
        color="grey",
        alpha=0.15,
        label="illuminated bins",
    )
    ax.set_xlabel("r (internal length units)")
    ax.set_ylabel("Temperature (K)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(opt.output, dpi=150)
    print(f"Plot saved to {opt.output}")

    print()
    print("--- Tier 2: quantitative photoelectric heating rate ---")
    rate_failures = []
    tier2_ran = False
    if opt.rate_on_snapshot is None or opt.rate_off_snapshot is None:
        print(
            "Skipped: pass --rate-on-snapshot/--rate-off-snapshot (a "
            "SEPARATE, fine-cadence run pair -- see docstring) to run it."
        )
    else:
        rate_on_files = sorted(glob.glob(opt.rate_on_snapshot))
        rate_off_files = sorted(glob.glob(opt.rate_off_snapshot))
        if len(rate_on_files) < 2 or len(rate_off_files) < 2:
            rate_failures.append(
                "need at least 2 snapshots per run for the rate check, found "
                f"{len(rate_on_files)} (on) and {len(rate_off_files)} (off)"
            )
            tier2_ran = True
        else:
            tier2_ran = True
            on_snaps = [load_snapshot(f) for f in rate_on_files]

            # Precondition: the on- and off-run start from the same gas
            # state (same IC, before any injection). Checked at snapshot 0
            # specifically (never the rate window below, which can start
            # later -- see field-arrival search) since that is the one
            # point both runs are guaranteed identical at. Verified, not
            # assumed: the on/off subtraction the rate measurement relies
            # on is only valid if this holds.
            on0_for_precondition = load_snapshot(rate_on_files[0])
            off0 = load_snapshot(rate_off_files[0])
            align_by_id(on0_for_precondition, off0)
            density_reldiff = np.abs(
                on0_for_precondition["density_cgs"] - off0["density_cgs"]
            ) / np.abs(off0["density_cgs"])
            temperature_reldiff = np.abs(
                on0_for_precondition["temperature"] - off0["temperature"]
            ) / np.abs(off0["temperature"])
            max_density_reldiff = float(np.max(density_reldiff))
            max_temperature_reldiff = float(np.max(temperature_reldiff))
            print(
                f"Precondition: max on/off relative difference at "
                f"snapshot 0, density {max_density_reldiff:.4g}, "
                f"temperature {max_temperature_reldiff:.4g} "
                f"(bar {opt.precondition_tol:.4g})"
            )
            for message in (
                check_gate(
                    "precondition density",
                    max_density_reldiff,
                    max_density_reldiff > opt.precondition_tol,
                    f"on/off density differs by {max_density_reldiff:.4g} "
                    f"at snapshot 0, above {opt.precondition_tol:.4g}: the "
                    "two runs did not start from the same state",
                ),
                check_gate(
                    "precondition temperature",
                    max_temperature_reldiff,
                    max_temperature_reldiff > opt.precondition_tol,
                    "on/off temperature differs by "
                    f"{max_temperature_reldiff:.4g} at snapshot 0, above "
                    f"{opt.precondition_tol:.4g}: the two runs did not "
                    "start from the same state",
                ),
            ):
                if message is not None:
                    rate_failures.append(message)

            # Field arrival: with time_first=0, snapshot 0 is written
            # before any injection (the ISRF module has not run a step
            # yet), so G0 = 0 everywhere there. Use the first snapshot at
            # which some particle's G0 is already non-negligible as the
            # window's start, rather than assuming index 0.
            start_index = None
            for index, snapshot in enumerate(on_snaps):
                g0 = (
                    C_LIGHT_CGS
                    * snapshot["density_cgs"]
                    * (snapshot["u_pe_cgs"] + snapshot["u_lw_cgs"])
                    / RADIATION_HABING_FLUX_CGS
                )
                if np.nanmax(g0) > 0.0:
                    start_index = index
                    break
            on1 = off1 = None
            if start_index is None:
                rate_failures.append(
                    "the on-run's ISRF field is zero in every snapshot"
                )
                start_index = 0
            elif start_index + 1 >= len(on_snaps) or start_index + 1 >= len(
                rate_off_files
            ):
                rate_failures.append(
                    f"the field only arrives at snapshot {start_index}, "
                    "leaving no further snapshot to difference against"
                )
            print(
                f"Field-arrival snapshot: {start_index} (on-run's "
                "earliest snapshot with a non-negligible G0 anywhere; "
                f"window is [{start_index}, {start_index + 1}])"
            )

            on0 = on_snaps[start_index]
            off0 = load_snapshot(rate_off_files[start_index])
            if start_index + 1 < len(on_snaps) and start_index + 1 < len(
                rate_off_files
            ):
                on1 = on_snaps[start_index + 1]
                off1 = load_snapshot(rate_off_files[start_index + 1])
                align_by_id(on0, on1)
                align_by_id(off0, off1)
                align_by_id(on0, off0)
                align_by_id(on1, off1)

            if on1 is not None:
                dt = on1["time_cgs"] - on0["time_cgs"]
                off_dt = off1["time_cgs"] - off0["time_cgs"]
                print(
                    f"Window: t0={on0['time_cgs']:.4e} s, "
                    f"t1={on1['time_cgs']:.4e} s, dt={dt:.4e} s "
                    f"(off-run dt={off_dt:.4e} s)"
                )
            else:
                print("Window: unavailable (see failure above)")

            g0_0 = (
                C_LIGHT_CGS
                * on0["density_cgs"]
                * (on0["u_pe_cgs"] + on0["u_lw_cgs"])
                / RADIATION_HABING_FLUX_CGS
            )

            # Illuminated, sub-cutoff particles only: same field-threshold
            # convention as the Tier 1 bins, applied per particle, plus
            # the igammah=2 temperature cutoff.
            illuminated_particles = g0_0 > opt.field_threshold * np.nanmax(g0_0)
            below_cutoff = on0["temperature"] < GAMMAH_TEMPERATURE_CEILING_K
            selected = illuminated_particles & below_cutoff
            n_selected = int(selected.sum())
            print(
                f"{n_selected}/{g0_0.size} particles are illuminated "
                f"(G0 > {opt.field_threshold:.1e} x peak) and below the "
                f"{GAMMAH_TEMPERATURE_CEILING_K:.1e} K igammah=2 cutoff"
            )

            if on1 is None:
                pass  # already recorded as a failure above
            elif n_selected == 0:
                rate_failures.append(
                    "no particle is both illuminated and below the "
                    "igammah=2 temperature cutoff -- cannot measure a rate"
                )
            else:
                g0_1 = (
                    C_LIGHT_CGS
                    * on1["density_cgs"]
                    * (on1["u_pe_cgs"] + on1["u_lw_cgs"])
                    / RADIATION_HABING_FLUX_CGS
                )
                g0_drift = np.abs(g0_1[selected] - g0_0[selected]) / g0_0[selected]
                median_g0_drift = float(np.median(g0_drift))
                print(
                    f"G0 drift over the window (truncation-error "
                    f"diagnostic), median {median_g0_drift:.4g}, 16th to "
                    f"84th percentile {np.percentile(g0_drift, 16):.4g} to "
                    f"{np.percentile(g0_drift, 84):.4g}"
                )

                predicted = gamma_pe_over_rho_predicted(
                    on0["density_cgs"][selected],
                    on0["u_pe_cgs"][selected],
                    on0["u_lw_cgs"][selected],
                    on0["metallicity"][selected],
                    opt.hydrogen_fraction,
                )

                # Hard precondition, not a precision knob: the on/off
                # cancellation below is only valid while the heating
                # perturbation is still small relative to u itself (see
                # docstring, "Why this needs its own...").
                perturbation = np.abs(predicted * dt / on0["u"][selected])
                median_perturbation = float(np.median(perturbation))
                print(
                    f"Perturbation Gamma_PE*dt/u over the window "
                    f"(linearity precondition), median "
                    f"{median_perturbation:.4g} (bar {opt.perturbation_tol:.4g})"
                )
                message = check_gate(
                    "linearity precondition",
                    median_perturbation,
                    median_perturbation > opt.perturbation_tol,
                    f"median Gamma_PE*dt/u {median_perturbation:.4g} above "
                    f"{opt.perturbation_tol:.4g}: the window is not in the "
                    "linear regime this method needs (use a finer "
                    "delta_time/dt_max)",
                )
                if message is not None:
                    rate_failures.append(message)

                if message is None:
                    dudt_on = (on1["u"][selected] - on0["u"][selected]) / dt
                    dudt_off = (off1["u"][selected] - off0["u"][selected]) / off_dt
                    measured = dudt_on - dudt_off

                    relative_error = np.abs(measured - predicted) / np.abs(predicted)
                    median_relative_error = float(np.median(relative_error))
                    print(
                        f"Measured Gamma_PE/rho, median "
                        f"{np.median(measured):.4e} erg/g/s; predicted, "
                        f"median {np.median(predicted):.4e} erg/g/s"
                    )
                    print(
                        f"Relative error, median {median_relative_error:.4g}, "
                        "16th to 84th percentile "
                        f"{np.percentile(relative_error, 16):.4g} to "
                        f"{np.percentile(relative_error, 84):.4g} "
                        f"(bar {opt.rate_tol:.4g})"
                    )

                    message = check_gate(
                        "rate relative error",
                        median_relative_error,
                        median_relative_error > opt.rate_tol,
                        f"relative error {median_relative_error:.4g} above "
                        f"{opt.rate_tol:.4g}",
                    )
                    if message is not None:
                        rate_failures.append(message)

    if tier2_ran:
        if rate_failures:
            for failure in rate_failures:
                print(f"FAIL: {failure}")
        else:
            print("PASS")

    if not (ok_sign and ok_magnitude) or rate_failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
