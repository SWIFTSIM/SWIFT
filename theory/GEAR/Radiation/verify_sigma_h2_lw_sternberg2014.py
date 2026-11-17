"""Verify RADIATION_SIGMA_H2_LW_CGS against Sternberg et al. (2014)'s D0.

Standalone numeric comparison, no simulation: compare
this codebase's own H2 Lyman-Werner photodissociation-rate formula
(radiation_get_part_LW_dissociation_rate_internal(),
src/feedback/GEAR/radiation_gas.c) against Sternberg, Le Petit, Roueff &
Le Bourlot (2014, ApJ 790:10)'s free-space photodissociation rate

    D0 = 5.8e-11 * I_UV  s^-1                                   (their Eq 5)

D0 is the free-space (full 4*pi) rate, NOT D(0) = D0/2, which is a
slab-surface geometric factor (their Eq 6) that does not apply here.

Sternberg's I_UV is normalized on the Draine (1978) spectrum integrated over
the LW band (912-1108 A, 11.3-13.6 eV) specifically, confirmed from the
paper's own Section 2.1 and Eq 2-4 (not the broader 6-13.6 eV band some
other "G0" conventions use). This codebase's k_diss formula likewise only
ever reads the LW-band field (p->feedback_data.u_LW), not the PE band, so
the two are compared on the same band. Draine/Habing normalization
(verified from the paper's own p.4 footnote 4): the Draine energy
density is 1.7x the Habing (1968) LW-band estimate, so

    I_UV (Draine) = G0 (Habing) / 1.7

This script sweeps a Habing-convention G0 taken to represent the LW-band
flux alone (the quantity radiation_get_part_LW_dissociation_rate_internal()
actually consumes), NOT this codebase's separate
radiation_get_part_isrf_habing() getter, which sums the PE+LW bands for an
unrelated purpose (photoelectric heating) and would double-count non-LW
flux here.

The rate comparison is informational only (used to judge whether
RADIATION_SIGMA_H2_OVER_E_LW_CGS, and so the cross section it is quoted
through, is plausible). One check gates the exit status: the quotient in the
header must lie between the smallest and the largest of the source values
below, which do not come from the header. With the LW-fraction-of-Habing
correction below applied, k_diss/D0 = 1.008. The constant is the mean of three
computations (an H2 line-list integral, Sternberg et al. 2014 and Draine
and Bertoldi 1996 Table 2), derived in 04_isrf_packet_transport.tex, so this
ratio is a cross-check of the field conversion, not the anchor of the constant.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

RADIATION_H = (
    Path(__file__).resolve().parents[3] / "src" / "feedback" / "GEAR" / "radiation.h"
)

# Sternberg et al. (2014), Eq 5: free-space (full 4*pi) LW photodissociation
# rate for the Draine spectrum, in units of I_UV.
STERNBERG_D0_PER_IUV_CGS = 5.8e-11  # s^-1

# Verified from Sternberg et al. (2014), p.4 footnote 4.
DRAINE_OVER_HABING = 1.7

# Source values of sigma_H2/E_LW, cm^2 erg^-1, from 04_isrf_packet_transport.tex:
# the line-list integral (integrate_h2_line_list_lw.py), the rate of Sternberg
# et al. (2014) over the LW energy flux of their field, and the two levels of
# Draine and Bertoldi (1996), Table 2 (rate x dissociation probability) over the
# LW energy flux of the Draine field at chi = 1.
SOURCE_QUOTIENTS_CGS = {
    "H2 line list": 1.4515e-7,
    "Sternberg et al. 2014": 5.8e-11 / 3.9596e-4,
    "Draine and Bertoldi 1996, level (0,0)": 2.77e-10 * 0.117 / 2.3245e-4,
    "Draine and Bertoldi 1996, level (0,1)": 2.79e-10 * 0.119 / 2.3245e-4,
}

# The header quotes the cross section to eight digits.
QUOTED_PAIR_RELATIVE_LIMIT = 0.5e-7 / 2.8196309

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
        The macro name, e.g. "RADIATION_SIGMA_H2_OVER_E_LW_CGS".

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
    G0_LW: float,
    sigma_over_E_LW_cgs: float,
    habing_flux_cgs: float,
) -> float:
    """Compute this module's H2 photodissociation rate at a given LW field.

    Mirrors the production function, starting from a Habing-convention
    LW-band flux instead of a particle's stored u_LW field: the rate is the
    LW energy flux times the one constrained quotient, with no separate
    cross section or photon energy.

    Parameters
    ----------
    G0_LW : float
        LW-band flux in Habing units (see module docstring for why this is
        LW-only, not this codebase's combined-band G0 getter).
    sigma_over_E_LW_cgs : float
        RADIATION_SIGMA_H2_OVER_E_LW_CGS, cm^2 erg^-1.
    habing_flux_cgs : float
        RADIATION_HABING_FLUX_CGS, erg/s/cm^2.

    Returns
    -------
    float
        k_diss, s^-1.
    """
    return sigma_over_E_LW_cgs * G0_LW * habing_flux_cgs


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


def main() -> int:
    """Run the comparison and the source-bracket check.

    Returns
    -------
    int
        1 if the header quotient lies outside the range of the source values,
        0 otherwise.
    """
    # The rate reads the quotient sigma_H2/E_LW and nothing else, so that is
    # what this verification anchors; the pair below is how it is quoted.
    sigma_over_E_LW_cgs = _read_radiation_h_constant("RADIATION_SIGMA_H2_OVER_E_LW_CGS")
    sigma_h2_lw_cgs = _read_radiation_h_constant("RADIATION_SIGMA_H2_LW_CGS")
    E_LW_photon_ev = _read_radiation_h_constant("RADIATION_LW_PHOTON_ENERGY_EV")
    habing_flux_cgs = _read_radiation_h_constant("RADIATION_HABING_FLUX_CGS")

    # RADIATION_HABING_FLUX_CGS is combined-band; rescale to LW-only.
    LW_FRACTION_OF_HABING = 0.149
    habing_flux_LW_only_cgs = habing_flux_cgs * LW_FRACTION_OF_HABING

    print(f"Read from {RADIATION_H}:")
    print(
        f"  RADIATION_SIGMA_H2_OVER_E_LW_CGS = {sigma_over_E_LW_cgs:.6e} "
        "cm^2 erg^-1 (the constrained quotient)"
    )
    print(
        f"  RADIATION_SIGMA_H2_LW_CGS   = {sigma_h2_lw_cgs:.4e} cm^2 "
        "(derived, reported only)"
    )
    print(f"  RADIATION_LW_PHOTON_ENERGY_EV = {E_LW_photon_ev:.3f} eV")
    quoted_quotient = sigma_h2_lw_cgs / (E_LW_photon_ev * ELECTRON_VOLT_CGS)
    drift = quoted_quotient / sigma_over_E_LW_cgs - 1.0
    if abs(drift) > QUOTED_PAIR_RELATIVE_LIMIT:
        print(
            f"  WARNING: the quoted pair gives {quoted_quotient:.6e} cm^2 "
            f"erg^-1, {drift * 100.0:+.4f} per cent off the quotient the rate "
            "uses. One of the three was changed without the others."
        )
    print(
        f"  RADIATION_HABING_FLUX_CGS   = {habing_flux_cgs:.4e} erg/s/cm^2 (combined PE+LW)"
    )
    print(
        f"  LW fraction of Habing band  = {LW_FRACTION_OF_HABING} (verified 2026-09-06)"
    )
    print(f"  -> LW-only equivalent       = {habing_flux_LW_only_cgs:.4e} erg/s/cm^2")
    print()

    G0_sweep = [0.1, 1.0, 10.0, 100.0]
    print(
        f"{'G0_LW':>10} {'I_UV':>10} {'k_diss (ours)':>16} "
        f"{'D0 (Sternberg)':>16} {'ratio k_diss/D0':>18}"
    )
    ratios = []
    for G0_LW in G0_sweep:
        k_diss = our_k_diss_cgs(G0_LW, sigma_over_E_LW_cgs, habing_flux_LW_only_cgs)
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
        "spectral-shape-dependent discrepancy, it only confirms the "
        "ratio is a fixed multiplicative offset, not a scaling error)."
    )
    print(
        f"k_diss/D0 = {ratios[0]:.3f} (LW-fraction-corrected; informational, "
        "see module docstring; not a pass/fail gate). Expected near 1.01 "
        "for the three-source constant."
    )
    lowest = min(SOURCE_QUOTIENTS_CGS.values())
    highest = max(SOURCE_QUOTIENTS_CGS.values())
    print()
    for name, value in SOURCE_QUOTIENTS_CGS.items():
        print(f"  source {name:<40s} {value:.4e} cm^2 erg^-1")
    inside = lowest <= sigma_over_E_LW_cgs <= highest
    print(
        f"Source bracket: header {sigma_over_E_LW_cgs:.6e} vs "
        f"[{lowest:.4e}, {highest:.4e}]: {'PASS' if inside else 'FAIL'}"
    )
    return 0 if inside else 1


if __name__ == "__main__":
    sys.exit(main())
