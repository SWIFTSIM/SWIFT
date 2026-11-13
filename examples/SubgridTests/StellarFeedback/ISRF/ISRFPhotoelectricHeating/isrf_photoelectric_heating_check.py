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
codebase injects, against a `with_interstellar_radiation_field: 1` run compared
to an otherwise-identical `with_interstellar_radiation_field: 0` run (same IC,
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

Without a species network (COOLING_GRACKLE_MODE 0, `ispecies=0`), Grackle forms `n_H` as
``HydrogenFractionByMass * (rho - rho*Z) / m_H`` (`cool1d_multi_g.F`,
the `imetal=1` branch of the `ispecies=0` temperature block), not the
naive `HydrogenFractionByMass * rho / m_H`; this script uses that
expression. The library's heating at `ispecies=0` is nevertheless 0.9418 of
the Eq. 1 value for any Z (0.001 to 0.04) and T (100 and 1000 K), in the
official library and in the fork; the Fortran line that sets this factor was
not found. Note `n_H / rho` is then a pure constant
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
heating-off run's ISRF fields are zero: `with_interstellar_radiation_field`
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
condition, `Gamma_PE*dt/u` of order 4e-3 (median, 39 illuminated
sub-cutoff particles) and `G0` unchanged between the two snapshots to
float32 precision -- a genuinely linear, cadence-converged regime. This
is a SEPARATE run pair from Tier 1's equilibrium `heating_on`/
`heating_off` (`--rate-on-snapshot`/`--rate-off-snapshot`; Tier 2 is
skipped, not failed, when these are not given -- see README), since
Tier 1 needs a run long enough to reach equilibrium and Tier 2 needs one
short enough to still be in the linear regime; no single run pair serves
both.

Heating models
--------------
``--pe-model`` selects the reference formula; each model is compared with
its own formula only. The heating run minus the reference run isolates the
heating only if both runs carry the same non-heating terms of the dust
chemistry that SWIFT switches on with the ISRF (grain recombination cooling
and a small dust cooling that depends little on G0). ``--reference-field``
says which reference run is used:

* ``on`` (default for the geajay models): a run with the ISRF on and the
  heating switched off or negligible. For geajay it is ``pe_model=none``. For
  Grackle it can be a ``wolfire1995`` run on neutral gas, whose heating is
  negligible. Both runs carry the same terms and none is subtracted.
* ``off`` (default for the Grackle models): the ISRF-off run. It lacks those
  terms, so the reference subtracts the Wolfire et al. (1995) Eq. 9 grain
  recombination cooling. The small dust cooling is not modelled: it is
  1.08e-3 of the heating at G0 = 520 in a one-zone call of the Grackle
  library, so this choice cannot pass a bar of that size.

The electron density of the models that read it is the snapshot's ``e``
field.

Pass bar
--------
The bar is a judgment bound, not a derived error. It has the form

    bar = 0.5 * P + 4 * 2^-24 / P + ``--rate-extra-tol``

with ``P`` the median ``Gamma_PE * dt / u`` of the window. The first term
allows for the response of the cooling to the higher energy of the heating
run; the coefficient 0.5 is assumed, not derived. The second term assumes
the rounding of four float32 energies, each 2^-24 of ``u``. SPHENIX holds
``u`` in float32 and updates it at every step, so the real rounding is larger
and grows with the number of steps in the window; a run with ten steps per
window gave an error 2.4 times this term. The measured errors of the
geajay models stay at about 0.1 per cent while ``P`` changes by a factor 8,
which neither term describes. ``--rate-extra-tol`` (default 0) is for terms
the reference does not model. The bar was written before the runs it judges
and is not changed after them. The check also needs a finite measurement, a
bar below 0.5, ``P`` above 0 and below ``--perturbation-tol``, and identical
on/off starts.

``--reference-scale`` and ``--inject-nan`` are positive controls: a scaled
reference, a reference of another model and a NaN in a measured value must
each fail.
"""

import argparse
import glob
import sys
from typing import Dict, Optional, Tuple

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

# Heating models the check knows a reference formula for. Each model is
# compared with its own formula only.
PE_MODELS = (
    "grackle-constant",
    "grackle-wolfire1995",
    "geajay-fixed-rate",
    "geajay-constant-efficiency",
    "geajay-grain-charge",
    "geajay-density-epsilon",
)
GEAJAY_FIXED_RATE_CGS = 8.5e-26  # erg/s per H atom at the reference dust-to-gas ratio
FLOAT32_EPS = 2.0**-24  # relative spacing of the snapshot floats, half ulp
RESOLVABLE_BAR_MAX = 0.5  # a bar this large cannot tell a model from its absence


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--on-snapshot",
        default="snap_on/snapshot_*.hdf5",
        help="Glob pattern for the with_interstellar_radiation_field=1 run's "
        "snapshots (default: %(default)s)",
    )
    parser.add_argument(
        "--off-snapshot",
        default="snap_off/snapshot_*.hdf5",
        help="Glob pattern for the with_interstellar_radiation_field=0 run's "
        "snapshots (default: %(default)s)",
    )
    parser.add_argument("--n-bins", type=int, default=12, help="Number of radial bins.")
    parser.add_argument(
        "--field-threshold",
        type=float,
        default=1e-3,
        help="A radial bin is considered field-illuminated if the on-run's "
        "mean (u_PE+u_LW) there exceeds this fraction of the innermost "
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
        help="Glob pattern for the FINE-CADENCE with_interstellar_radiation_field=1 "
        "run's snapshots, for Tier 2 (the quantitative rate check). A "
        "SEPARATE run from --on-snapshot: see docstring for why. Tier 2 is "
        "skipped, not failed, if this (and --rate-off-snapshot) are not "
        "given.",
    )
    parser.add_argument(
        "--rate-off-snapshot",
        default=None,
        help="Glob pattern for the FINE-CADENCE with_interstellar_radiation_field=0 "
        "run's snapshots, for Tier 2. See --rate-on-snapshot.",
    )
    parser.add_argument(
        "--pe-model",
        choices=PE_MODELS,
        default="grackle-constant",
        help="Heating model of the run; its own reference formula is used "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--reference-field",
        choices=("on", "off"),
        default=None,
        help="Whether the reference run has the ISRF on (heating off or "
        "negligible) or off. Default: on for the geajay models, off for the "
        "Grackle models.",
    )
    parser.add_argument(
        "--rate-extra-tol",
        type=float,
        default=0.0,
        help="Relative allowance for terms the reference does not model, "
        "added to the bar derived from the window (default: %(default)s).",
    )
    parser.add_argument(
        "--reference-scale",
        type=float,
        default=1.0,
        help="Positive control: multiply the reference by this factor; a "
        "factor other than 1 must fail the rate check (default: %(default)s).",
    )
    parser.add_argument(
        "--inject-nan",
        action="store_true",
        help="Positive control: write a NaN into one measured value; the "
        "rate check must fail.",
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
        u_pe = gas["PESpecificEnergies"][:]
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
            if "H2I" in gas:
                inv_mu = inv_mu + 0.5 * gas["H2I"][:]
            mu = 1.0 / np.clip(inv_mu, 1e-10, None)
        else:
            mu = np.full(pos.shape[0], 1.0 / 0.82)

        # Species are mass fractions; "e" is in units of the proton mass.
        density_cgs_raw = gas["Densities"][:] * density_to_cgs
        if "e" in gas:
            electron_density_cgs = gas["e"][:] * density_cgs_raw / M_P_CGS
        else:
            electron_density_cgs = None
        if "HI" in gas:
            hydrogen_mass_fraction = gas["HI"][:] + gas["HII"][:]
            for name in ("HM", "H2I", "H2II"):
                if name in gas:
                    hydrogen_mass_fraction = hydrogen_mass_fraction + gas[name][:]
            hydrogen_density_cgs = hydrogen_mass_fraction * density_cgs_raw / M_P_CGS
        else:
            hydrogen_density_cgs = None

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
        electron_density_cgs=electron_density_cgs,
        hydrogen_density_cgs=hydrogen_density_cgs,
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


def pe_rate_over_rho_predicted(
    model: str,
    snapshot: Dict[str, np.ndarray],
    selected: np.ndarray,
    hydrogen_fraction: float,
    subtract_recombination: bool,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the reference photoelectric heating per unit mass of one model.

    Each model uses only its own published formula, with the gas state and
    the field G0 read from the heating-on snapshot.

    Parameters
    ----------
    model : str
        One of ``PE_MODELS``.
    snapshot : dict of str to numpy.ndarray
        Heating-on snapshot from :func:`load_snapshot`.
    selected : numpy.ndarray
        Boolean mask of the particles to evaluate.
    hydrogen_fraction : float
        GrackleCooling:HydrogenFractionByMass (Grackle models only).
    subtract_recombination : bool
        Whether to return the Wolfire et al. (1995) Eq. 9 recombination
        cooling (a reference run with the ISRF off).

    Returns
    -------
    tuple of numpy.ndarray
        ``(Gamma_PE / rho, Lambda_rec / rho)`` in erg g^-1 s^-1. The second
        entry is the Wolfire et al. (1995) Eq. 9 grain recombination cooling
        when ``subtract_recombination`` is set and the snapshot has
        electrons; otherwise it is zero.
    """
    rho = snapshot["density_cgs"][selected]
    z = snapshot["metallicity"][selected]
    temperature = snapshot["temperature"][selected]
    g0 = (
        C_LIGHT_CGS
        * rho
        * (snapshot["u_pe_cgs"][selected] + snapshot["u_lw_cgs"][selected])
        / RADIATION_HABING_FLUX_CGS
    )
    z_prime = z / RADIATION_GRACKLE_SOLAR_METAL_FRACTION
    if snapshot["hydrogen_density_cgs"] is not None:
        n_h = snapshot["hydrogen_density_cgs"][selected]
    else:
        n_h = hydrogen_fraction * rho * (1.0 - z) / M_P_CGS
    n_e = None
    if model in ("grackle-wolfire1995", "geajay-grain-charge"):
        if snapshot["electron_density_cgs"] is None:
            raise RuntimeError(
                f"{model} needs the electron density: the snapshot has no 'e' "
                "field (a cooling mode without species)"
            )
        n_e = snapshot["electron_density_cgs"][selected]

    recombination = np.zeros_like(rho)
    if model in ("grackle-constant", "geajay-constant-efficiency"):
        heating = GRACKLE_GAMMA_HA_CGS * GRACKLE_EPSILON * g0 * n_h * z_prime
    elif model == "geajay-fixed-rate":
        heating = GEAJAY_FIXED_RATE_CGS * n_h * z_prime
    elif model == "geajay-density-epsilon":
        epsilon = np.minimum(0.07, 0.0149 * n_h**0.235)
        heating = 1.3e-24 * epsilon * g0 * n_h * z_prime
    elif model in ("grackle-wolfire1995", "geajay-grain-charge"):
        with np.errstate(divide="ignore", invalid="ignore"):
            x = g0 * np.sqrt(temperature) / n_e
            epsilon = 4.9e-2 / (1.0 + (x / 1925.0) ** 0.73) + (
                3.7e-2 * (temperature / 1.0e4) ** 0.7 / (1.0 + x / 5000.0)
            )
        epsilon = np.where(n_e > 0.0, epsilon, 0.0)
        heating = GRACKLE_GAMMA_HA_CGS * epsilon * g0 * n_h * z_prime
    else:
        raise ValueError(f"unknown PE model {model!r}")
    # Without a species network (grackle_0) the electron density is internal
    # to the library, so the recombination cooling is not modelled.
    has_electrons = snapshot["electron_density_cgs"] is not None
    if subtract_recombination and has_electrons:
        with np.errstate(divide="ignore", invalid="ignore"):
            x = g0 * np.sqrt(temperature) / snapshot["electron_density_cgs"][selected]
            beta = 0.74 / temperature**0.068
            recombination = (
                4.65e-30
                * temperature**0.94
                * x**beta
                * snapshot["electron_density_cgs"][selected]
                * n_h
                * z_prime
            )
        n_e_all = snapshot["electron_density_cgs"][selected]
        recombination = np.where(n_e_all > 0.0, recombination, 0.0)
    return heating / rho, recombination / rho


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
    reference_field = opt.reference_field or (
        "on" if opt.pe_model.startswith("geajay") else "off"
    )

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

                heating, recombination = pe_rate_over_rho_predicted(
                    opt.pe_model,
                    on0,
                    selected,
                    opt.hydrogen_fraction,
                    reference_field == "off",
                )
                predicted_unscaled = heating - recombination
                predicted = predicted_unscaled * opt.reference_scale
                print(
                    f"Model {opt.pe_model}: heating median "
                    f"{np.median(heating):.4e}, recombination cooling median "
                    f"{np.median(recombination):.4e} erg/g/s"
                )

                # Hard precondition, not a precision knob: the on/off
                # cancellation below is only valid while the heating
                # perturbation is still small relative to u itself (see
                # docstring, "Why this needs its own...").
                perturbation = np.abs(heating * dt / on0["u"][selected])
                median_perturbation = float(np.median(perturbation))
                if not median_perturbation > 0.0:
                    rate_failures.append(
                        "the reference heating is zero in this window: "
                        "nothing to compare"
                    )
                    median_perturbation = np.nan
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

                    if opt.inject_nan:
                        measured = measured.copy()
                        measured[0] = np.nan
                    relative_error = np.abs(measured - predicted) / np.abs(predicted)
                    median_relative_error = float(np.median(relative_error))
                    # Bar from the mechanism: the forward difference of a
                    # rate that is constant over the window has a truncation
                    # error of at most half the perturbation; the snapshot
                    # floats limit each of the four stored energies to
                    # FLOAT32_EPS of u, against a signal of perturbation * u.
                    bar = (
                        0.5 * median_perturbation
                        + 4.0 * FLOAT32_EPS / median_perturbation
                        + opt.rate_extra_tol
                    )
                    message = check_gate(
                        "resolvability",
                        bar,
                        bar >= RESOLVABLE_BAR_MAX,
                        f"bar {bar:.4g} is not below {RESOLVABLE_BAR_MAX}: "
                        "the heating in this window is too small for the "
                        "snapshot floats to resolve, so no agreement can be "
                        "claimed",
                    )
                    if message is not None:
                        rate_failures.append(message)
                    print(
                        f"Measured Gamma_PE/rho, median "
                        f"{np.median(measured):.4e} erg/g/s; reference, "
                        f"median {np.median(predicted):.4e} erg/g/s"
                    )
                    print(
                        f"Relative error, median {median_relative_error:.4g}, "
                        "16th to 84th percentile "
                        f"{np.percentile(relative_error, 16):.4g} to "
                        f"{np.percentile(relative_error, 84):.4g} "
                        f"(bar {bar:.4g} = 0.5 x perturbation + "
                        "4 x 2^-24 / perturbation + "
                        f"{opt.rate_extra_tol:.4g})"
                    )

                    message = check_gate(
                        "rate relative error",
                        median_relative_error,
                        median_relative_error > bar,
                        f"relative error {median_relative_error:.4g} above "
                        f"{bar:.4g}",
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
