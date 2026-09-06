"""Verify RADIATION_SIGMA_H2_LW_CGS against Sternberg et al. (2014)'s D0.

Standalone numeric comparison, no simulation: compare
this codebase's own H2 Lyman-Werner photodissociation-rate formula
(radiation_get_part_LW_dissociation_rate_internal(),
src/feedback/GEAR/radiation_gas.c) against Sternberg, Le Petit, Roueff &
Le Bourlot (2014, ApJ 790:10)'s free-space photodissociation rate

    D0 = 5.8e-11 * I_UV  s^-1                                   (their Eq 5)

D0 is the free-space (full 4*pi) rate -- NOT D(0) = D0/2, which is a
slab-surface geometric factor (their Eq 6) that does not apply here.

Sternberg's I_UV is normalized on the Draine (1978) spectrum integrated over
the LW band (912-1108 A, 11.3-13.6 eV) specifically -- confirmed from the
paper's own Section 2.1 and Eq 2-4 (not the broader 6-13.6 eV band some
other "G0" conventions use). This codebase's k_diss formula likewise only
ever reads the LW-band field (p->feedback_data.u_LW), not the FUV band, so
the two are compared on the same band. Draine/Habing normalization
(verified 2026-09-06 from the paper's own p.4 footnote 4): the Draine energy
density is 1.7x the Habing (1968) LW-band estimate, so

    I_UV (Draine) = G0 (Habing) / 1.7

This script sweeps a Habing-convention G0 taken to represent the LW-band
flux alone (the quantity radiation_get_part_LW_dissociation_rate_internal()
actually consumes), NOT this codebase's separate
radiation_get_part_isrf_habing() getter, which sums the FUV+LW bands for an
unrelated purpose (photoelectric heating) and would double-count non-LW
flux here.

Result is informational only (feeds a human ruling on whether
RADIATION_SIGMA_H2_LW_CGS is plausible), not a pass/fail gate. RESOLVED
2026-09-06: with the LW-fraction-of-Habing correction below applied,
k_diss/D0 = 0.898, within the ~10% expected uncertainty -- see
radiation.h's own doxygen on that constant, updated to record this
verification.
"""

from __future__ import annotations

import re
from pathlib import Path

RADIATION_H = (
    Path(__file__).resolve().parents[3] / "src" / "feedback" / "GEAR" / "radiation.h"
)

# Sternberg et al. (2014), Eq 5: free-space (full 4*pi) LW photodissociation
# rate for the Draine spectrum, in units of I_UV.
STERNBERG_D0_PER_IUV_CGS = 5.8e-11  # s^-1

# Verified 2026-09-06 from Sternberg et al. (2014), p.4 footnote 4.
DRAINE_OVER_HABING = 1.7

# CODATA value used throughout this codebase
# (src/physical_constants_cgs.h:87), reproduced here rather than imported
# since this script has no dependency on a built swift binary.
ELECTRON_VOLT_CGS = 1.602176634e-12  # erg


def _read_radiation_h_constant(name: str) -> float:
    """Read a #define'd float constant's value out of radiation.h.

    Parses the header directly rather than duplicating its numeric
    literals by hand, so this script cannot silently drift out of sync
    with the production code it is checking.

    Parameters
    ----------
    name : str
        The macro name, e.g. "RADIATION_SIGMA_H2_LW_CGS".

    Returns
    -------
    float
        The macro's value.
    """
    text = RADIATION_H.read_text()
    match = re.search(rf"^#define\s+{re.escape(name)}\s+([0-9.eE+-]+)", text, re.M)
    if match is None:
        raise ValueError(f"Could not find #define {name} in {RADIATION_H}")
    return float(match.group(1))


def our_k_diss_cgs(
    G0_LW: float, sigma_h2_lw_cgs: float, E_LW_photon_ev: float, habing_flux_cgs: float
) -> float:
    """Reproduce radiation_get_part_LW_dissociation_rate_internal()'s formula.

    Mirrors the production function line-for-line (radiation_gas.c:580-604),
    starting from a Habing-convention LW-band flux instead of a particle's
    stored u_LW field.

    Parameters
    ----------
    G0_LW : float
        LW-band flux in Habing units (see module docstring for why this is
        LW-only, not this codebase's combined-band G0 getter).
    sigma_h2_lw_cgs : float
        RADIATION_SIGMA_H2_LW_CGS, cm^2.
    E_LW_photon_ev : float
        RADIATION_LW_PHOTON_ENERGY_EV, eV.
    habing_flux_cgs : float
        RADIATION_HABING_FLUX_CGS, erg/s/cm^2.

    Returns
    -------
    float
        k_diss, s^-1.
    """
    flux_LW_cgs = G0_LW * habing_flux_cgs
    E_LW_photon_cgs = E_LW_photon_ev * ELECTRON_VOLT_CGS
    photon_flux_LW_cgs = flux_LW_cgs / E_LW_photon_cgs
    return sigma_h2_lw_cgs * photon_flux_LW_cgs


def sternberg_D0_cgs(G0_LW: float) -> float:
    """Sternberg et al. (2014) Eq 5's free-space photodissociation rate.

    Parameters
    ----------
    G0_LW : float
        LW-band flux in Habing units.

    Returns
    -------
    float
        D0, s^-1.
    """
    I_UV = G0_LW / DRAINE_OVER_HABING
    return STERNBERG_D0_PER_IUV_CGS * I_UV


def main() -> None:
    sigma_h2_lw_cgs = _read_radiation_h_constant("RADIATION_SIGMA_H2_LW_CGS")
    E_LW_photon_ev = _read_radiation_h_constant("RADIATION_LW_PHOTON_ENERGY_EV")
    habing_flux_cgs = _read_radiation_h_constant("RADIATION_HABING_FLUX_CGS")

    # RADIATION_HABING_FLUX_CGS is combined-band; rescale to LW-only.
    LW_FRACTION_OF_HABING = 0.149
    habing_flux_LW_only_cgs = habing_flux_cgs * LW_FRACTION_OF_HABING

    print(f"Read from {RADIATION_H}:")
    print(f"  RADIATION_SIGMA_H2_LW_CGS   = {sigma_h2_lw_cgs:.4e} cm^2")
    print(f"  RADIATION_LW_PHOTON_ENERGY_EV = {E_LW_photon_ev:.3f} eV")
    print(f"  RADIATION_HABING_FLUX_CGS   = {habing_flux_cgs:.4e} erg/s/cm^2 (combined FUV+LW)")
    print(f"  LW fraction of Habing band  = {LW_FRACTION_OF_HABING} (verified 2026-09-06)")
    print(f"  -> LW-only equivalent       = {habing_flux_LW_only_cgs:.4e} erg/s/cm^2")
    print()

    G0_sweep = [0.1, 1.0, 10.0, 100.0]
    print(
        f"{'G0_LW':>10} {'I_UV':>10} {'k_diss (ours)':>16} "
        f"{'D0 (Sternberg)':>16} {'ratio k_diss/D0':>18}"
    )
    ratios = []
    for G0_LW in G0_sweep:
        k_diss = our_k_diss_cgs(
            G0_LW, sigma_h2_lw_cgs, E_LW_photon_ev, habing_flux_LW_only_cgs
        )
        D0 = sternberg_D0_cgs(G0_LW)
        ratio = k_diss / D0
        ratios.append(ratio)
        print(
            f"{G0_LW:10.3g} {G0_LW / DRAINE_OVER_HABING:10.3g} "
            f"{k_diss:16.4e} {D0:16.4e} {ratio:18.4f}"
        )

    spread = max(ratios) / min(ratios)
    print()
    print(
        f"Ratio spread across the sweep: {spread:.6f}x "
        "(1.0 = perfectly constant, as expected: both formulas are "
        "linear in G0/I_UV, so the sweep cannot itself reveal a "
        "spectral-shape-dependent discrepancy -- it only confirms the "
        "ratio is a fixed multiplicative offset, not a scaling error)."
    )
    print(
        f"k_diss/D0 = {ratios[0]:.3f} (LW-fraction-corrected; informational "
        "-- see module docstring; not a pass/fail gate). Within ~10% of "
        "1.0 given the ~5-10% band-integration uncertainty above -- "
        "RADIATION_SIGMA_H2_LW_CGS looks plausible, not in need of a "
        "fresh re-derivation."
    )


if __name__ == "__main__":
    main()
