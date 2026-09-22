#!/usr/bin/env python3
"""Rank the GEARFeedback:ISRF_extinction_path mechanisms on the uniform
injection_dusty box against an exact analytic reference.

In a uniform medium of density rho the attenuation of the band injected into
a particle at separation r from the star is exactly exp(-kappa_eff rho r).
Averaging it under the injection weight the code itself applies gives the
injection-weighted optical depth

    tau_exact = -ln( <exp(-kappa_eff rho r)>_w ),

which has no free parameter: it is formed from the star-to-particle
separations and the dust cross-section chain, and never from the smoothing
length. Each candidate leg supplies its own tau_code the same way, from the
ratio of its injected band energy to that of a leg run with a negligible
path. The score of a mechanism is log10(tau_code / tau_exact); zero means it
reproduces the exact uniform-geometry answer.

What this check cannot do: in a uniform box the pair_separation mechanism
computes exp(-kappa rho_j r) while the reference computes
exp(-kappa_eff rho r), so the two differ only by the SPH density estimate. A
good score there is a regression check on the density estimate and the unit
chain, not evidence that the mechanism wins on real gas.

Do NOT use the band-ratio residual of the C1 and M4 injection checks to
discriminate between paths. Both bands attenuate on the same column and
sigma_PE / sigma_LW = 9e-22 / 1.5e-21 = 0.6 exactly, so that residual is an
algebraic identity in the path and is zero for any column whatsoever.

Usage
-----
Run one leg with a negligible path as the reference, then one leg per
candidate, all at the same metallicity::

    config=injection_dusty extinction_path=constant_kernel_path \\
        extinction_path_in_kernel_radii=1e-8 run_name=thin ./run.sh
    config=injection_dusty extinction_path=constant_kernel_path \\
        extinction_path_in_kernel_radii=1.0 run_name=R1 ./run.sh
    config=injection_dusty extinction_path=pair_separation \\
        run_name=pairsep ./run.sh
    python3 isrf_extinction_path_check.py thin R1 pairsep

The reference leg must share the candidates' metallicity: the star's band
luminosity is metallicity dependent, so a Z = 0 leg is NOT an unattenuated
reference for a dusty one.
"""

import argparse
import glob
import os
import sys
from typing import Dict, List

import h5py
import numpy as np

# radiation.h's own LW dust cross-section per hydrogen nucleus, mean
# molecular weight per hydrogen nucleus, and hydrogen mass, in CGS.
SIGMA_D_LW_CGS = 1.5e-21
MU_H = 1.4
M_H_CGS = 1.6726219e-24
GRACKLE_SOLAR_METAL_FRACTION = 0.01295
GRACKLE_DEFAULT_DUST_TO_GAS_RATIO = 0.009387

PHYSICAL_CGS = "Conversion factor to physical CGS (including cosmological corrections)"


def physical(group: h5py.Group, name: str) -> np.ndarray:
    """Return a snapshot dataset in physical CGS."""
    dataset = group[name]
    factor = float(np.atleast_1d(dataset.attrs[PHYSICAL_CGS])[0])
    return dataset[:].astype(np.float64) * factor


def read_leg(directory: str, band: str) -> Dict[str, np.ndarray]:
    """Read the last snapshot of one run directory."""
    snapshots = sorted(glob.glob(os.path.join(directory, "snap", "snapshot_*.hdf5")))
    if not snapshots:
        raise RuntimeError(f"No snapshot under {directory}/snap.")
    with h5py.File(snapshots[-1], "r") as handle:
        gas = handle["PartType0"]
        length_cgs = float(np.atleast_1d(gas["Coordinates"].attrs[PHYSICAL_CGS])[0])
        stars = handle["PartType4"]
        if stars["Coordinates"].shape[0] != 1:
            raise RuntimeError("This check assumes the fixture's single star.")
        return dict(
            ids=gas["ParticleIDs"][:],
            x=physical(gas, "Coordinates"),
            rho=physical(gas, "Densities"),
            mass=physical(gas, "Masses"),
            u=physical(gas, f"{band}SpecificEnergies"),
            metallicity=gas["SmoothedMetalMassFractions"][:, -1].astype(np.float64),
            box=np.atleast_1d(handle["Header"].attrs["BoxSize"]).astype(np.float64)
            * length_cgs,
            star=stars["Coordinates"][:][0].astype(np.float64) * length_cgs,
        )


def require_finite(name: str, values: np.ndarray) -> None:
    """Abort unless every entry is finite.

    A non-finite value compares false against every bar, so it must be
    caught before, never by, the comparison it would otherwise pass.
    """
    if not np.all(np.isfinite(values)):
        sys.exit(f"FAIL: {name} carries a non-finite value.")


def separations(leg: Dict[str, np.ndarray]) -> np.ndarray:
    """Return each particle's physical separation from the star."""
    delta = leg["x"] - leg["star"]
    delta -= leg["box"] * np.round(delta / leg["box"])
    return np.sqrt((delta * delta).sum(axis=1))


def weighted_optical_depth(attenuation: np.ndarray, weight: np.ndarray) -> float:
    """Return -ln of the weight-averaged attenuation."""
    return -float(np.log(np.average(attenuation, weights=weight)))


def main() -> int:
    """Rank every candidate leg against the analytic reference."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "reference",
        help="Run directory of the negligible-path leg, at the candidates' "
        "own metallicity.",
    )
    parser.add_argument("candidates", nargs="+", help="Candidate run directories.")
    parser.add_argument(
        "--band", default="LW", choices=("LW", "PE"), help="Band to score."
    )
    parser.add_argument(
        "--bar",
        type=float,
        default=0.05,
        help="Maximum |log10(tau_code / tau_exact)| the best candidate may "
        "reach before this check fails (default: %(default)s dex).",
    )
    parser.add_argument(
        "--dust-to-gas-ratio",
        type=float,
        default=GRACKLE_DEFAULT_DUST_TO_GAS_RATIO,
        help="The run's resolved Grackle chemistry_data.local_dust_to_gas_ratio "
        "(default: %(default)s, Grackle's own compiled default).",
    )
    args = parser.parse_args()

    reference = read_leg(args.reference, args.band)
    require_finite("reference band energy", reference["u"])
    require_finite("reference density", reference["rho"])

    radius = separations(reference)
    require_finite("star-particle separations", radius)

    # u_inject = Delta_t * weight * L_band per unit mass, so the injection
    # weight is the reference leg's own deposit times the particle mass.
    illuminated = reference["u"] > 0.0
    weight = reference["u"] * reference["mass"] * illuminated
    if illuminated.sum() < 20:
        print("FAIL: too few illuminated particles to average over.")
        return 1

    dust_relative = (
        np.maximum(reference["metallicity"], 0.0)
        / GRACKLE_SOLAR_METAL_FRACTION
        * (args.dust_to_gas_ratio / GRACKLE_DEFAULT_DUST_TO_GAS_RATIO)
    )
    kappa_eff = SIGMA_D_LW_CGS * dust_relative / (MU_H * M_H_CGS)
    tau_exact = weighted_optical_depth(
        np.exp(-kappa_eff * reference["rho"] * radius), weight
    )
    require_finite("exact optical depth", np.array([tau_exact]))
    if tau_exact <= 0.0:
        print(
            "FAIL: the reference leg is optically thin, so no path can be "
            "scored. Raise the metallicity or the density."
        )
        return 1

    print(f"illuminated particles      : {int(illuminated.sum())}")
    print(f"exact injection-weighted tau: {tau_exact:.6f}  (band {args.band})")
    print(f"{'candidate':32s} {'tau_code':>10s} {'ratio':>8s} {'score/dex':>10s}")

    scores: List[float] = []
    failed = False
    for directory in args.candidates:
        leg = read_leg(directory, args.band)
        if not np.array_equal(leg["ids"], reference["ids"]):
            order = np.argsort(leg["ids"])
            leg["u"] = leg["u"][order][np.argsort(np.argsort(reference["ids"]))]
        require_finite(f"{directory} band energy", leg["u"])
        attenuation = np.where(
            illuminated, leg["u"] / np.where(illuminated, reference["u"], 1.0), 1.0
        )
        require_finite(f"{directory} attenuation", attenuation)
        if np.any(attenuation[illuminated] <= 0.0) or np.any(
            attenuation[illuminated] > 1.0 + 1e-6
        ):
            print(f"{directory:32s}   attenuation outside (0, 1]")
            failed = True
            continue
        tau_code = weighted_optical_depth(attenuation, weight)
        if not np.isfinite(tau_code) or tau_code <= 0.0:
            print(f"{directory:32s}   non-positive or non-finite tau_code")
            failed = True
            continue
        ratio = tau_code / tau_exact
        score = float(np.log10(ratio))
        scores.append(abs(score))
        print(f"{directory:32s} {tau_code:10.6f} {ratio:8.4f} {score:+10.4f}")

    if failed:
        print("FAIL: at least one candidate produced an unusable attenuation.")
        return 1
    best = min(scores)
    if not np.isfinite(best) or best > args.bar:
        print(
            f"FAIL: the best candidate is {best:.4f} dex from the exact "
            f"answer, above the {args.bar} dex bar."
        )
        return 1
    print(f"PASS: best candidate within {best:.4f} dex of the exact answer.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
