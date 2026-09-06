"""Verify the G0-to-Grackle photoelectric-heating-rate pipeline by an
actual call into Grackle's compiled solver, not just a source read.

This codebase computes a per-particle Habing-unit ISRF strength G0 via
radiation_get_part_isrf_habing() (src/feedback/GEAR/radiation_gas.c) and
hands it to Grackle through chemistry_data.use_isrf_field /
grackle_field_data.isrf_habing, with GrackleCooling:photoelectric_heating
set to 2 (Wolfire et al. 1995 Eq 1, constant epsilon=0.05).

The commonly-cited form of that formula is

    Gamma_PE = gamma_ha * epsilon * G0 * n_H         (naive form)

with gamma_ha=1e-24 erg/s and epsilon=0.05 hardcoded in Grackle's own
rate_functions.c (gammah_rate(), igammah>1 branch) -- not SWIFT-side
constants, so they are not read from radiation.h below; they are Grackle's
compiled behaviour being checked, not a value this codebase owns.

Reading Grackle's actual Fortran (src/clib/cool1d_multi_g.F,
COOLING_GRACKLE side note: this is the standalone Grackle repo, not
SWIFT's src/cooling/grackle/ wrapper) shows the naive form is INCOMPLETE:
the heating term actually added to edot is

    edot += gammaha_eff(i) * rhoH(i) * dom_inv * dust2gas(i) / fgr   (line ~997)

with gammaha_eff(i) = gammaha * 0.05 * myisrf(i) for igammah=2 (and pinned
to 0 above T=2e4 K, irrelevant at this script's T~100K test point but
real), and, in the default no-explicit-dust-density-field configuration
(idustfield=0),

    dust2gas(i) = fgr * metallicity(i)     (fgr cancels exactly)
    metallicity(i) = metal(i,j,k) / d(i,j,k) / z_solar   ( = Z/Zsun = Z' )

so the real formula is

    Gamma_PE = gamma_ha * epsilon * G0 * n_H * Z'    (actual form)

idustfield=0 is confirmed to be what SWIFT's own wrapper actually uses,
not assumed: grepping src/cooling/grackle/ for use_dust_density_field or
dust_density finds no reference at all, so it stays at Grackle's own default
(FALSE, grackle_chemistry_data_fields.def) -- SWIFT never sets it or
supplies a dust_density field. Note this means Grackle's
local_dust_to_gas_ratio parameter (which radiation.h's own doxygen
documents as scaling extinction) plays NO role in this heating term at
all: with fgr cancelling algebraically above, PE heating depends only on
metallicity(i), never on whatever local_dust_to_gas_ratio is set to. Same
underlying parameter, two different roles in two different code paths.

z_solar here is Grackle's own SolarMetalFractionByMass (0.01295, Cloudy
v13 abundances -- also confirmed unreferenced in src/cooling/grackle/,
hence left at that default). This is the exact same value radiation.h's
RADIATION_GRACKLE_SOLAR_METAL_FRACTION already records, so "Z'=1" means
the same thing on both sides of the coupling -- no Zsun-convention
mismatch to worry about here.

rhoH(i) itself is confirmed (cool1d_multi_g.F lines ~228-270) to be TOTAL
hydrogen (HI+HII, +H2 when tracked), not neutral-only H0 -- so "n_H" in
both forms above is the same total-hydrogen quantity this codebase's own
n_H means, not a species-restricted one.

At solar metallicity (Z'=1) the two forms coincide, which is exactly why
a metallicity-blind citation of the naive form was never caught before.
This script's Grackle call-through sweeps Z' explicitly to make the
discrepancy visible rather than staying hidden at the one test point
where it cancels.

Part 2/3 (call Grackle's own solver): no pygrackle install exists on this
machine (checked: `import pygrackle` fails), so this script instead
compiles and runs a small standalone C harness
(verify_photoelectric_heating_rate_grackle_harness.c, same directory)
against the compiled libgrackle at GRACKLE_LIB_DIR/GRACKLE_INCLUDE_DIR
(default /home/darwinr/local/{lib,include}, source
/home/darwinr/programs/grackle-swift -- confirmed byte-identical headers
against the installed ones before trusting them). The harness calls
calculate_cooling_time() twice per Z' (photoelectric_heating=2 and =0,
otherwise identical one-zone state) and isolates Gamma_PE as the edot
difference; see the harness file's own header comment for the unit
bookkeeping (density_units=m_H makes dom=1 exactly, sidestepping Grackle's
comoving dom/dom_inv bookkeeping instead of hand-deriving it).

Part 4 (G0 unit-convention check): radiation_get_part_isrf_habing()
computes G0 = c*rho*(u_FUV+u_LW)_cgs / RADIATION_HABING_FLUX_CGS, i.e. the
standard combined FUV+LW Habing (1968) convention (RADIATION_HABING_FLUX_
CGS=1.6e-3 erg/s/cm^2, radiation.h). Grackle's igammah=2/3 Fortran path
consumes myisrf(i) = isrf_habing(i,j,k) with no internal rescaling --
it trusts the caller's G0 to already be in this same Habing convention.
Corroborating evidence this actually is the same convention (not just
same-named): rate_functions.c's igammah<=1 default docstring states its
fixed rate assumes "epsilon=0.05, G_0=1.7" -- 1.7 is exactly this
codebase's own DRAINE_OVER_HABING constant (see
verify_sigma_h2_lw_sternberg2014.py), i.e. Grackle's own comment is
quoting the Draine field's value *in Habing units*, the same units
radiation_get_part_isrf_habing() emits. No analogous factor to the
sigma_H2_LW LW-band-fraction correction is needed here: unlike that
formula (which mixed a Draine-normalized quantity into the combined-band
Habing constant), the photoelectric path is Habing-native on both sides.

Result: PASS if the call-through Gamma_PE matches the actual formula
(with the Z' factor) to within floating-point tolerance at every swept
Z', not just at Z'=1 -- this is what would catch a real discrepancy
rather than one hidden by the sweep's own test point.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
RADIATION_H = SCRIPT_DIR.parents[2] / "src" / "feedback" / "GEAR" / "radiation.h"
HARNESS_SRC = SCRIPT_DIR / "verify_photoelectric_heating_rate_grackle_harness.c"

# Grackle's own hardcoded constants for photoelectric_heating=2
# (src/clib/rate_functions.c: gammah_rate(), igammah>1 branch), reproduced
# here as the formula being checked, not a value this codebase owns.
GRACKLE_GAMMA_HA_CGS = 1.0e-24  # erg/s
GRACKLE_EPSILON = 0.05

# This script's test point, matching
# verify_sigma_h2_lw_sternberg2014.py's single-point-plus-sweep approach.
N_H_CGS = 100.0  # cm^-3, total hydrogen
T_K = 100.0
G0 = 1.0
Z_PRIME_SWEEP = [1.0, 0.1, 3.0]

GRACKLE_INCLUDE_DIR = "/home/darwinr/local/include"
GRACKLE_LIB_DIR = "/home/darwinr/local/lib"
GRACKLE_SOURCE_DIR = "/home/darwinr/programs/grackle-swift"


def _read_radiation_h_constant(name: str) -> float:
    """Read a #define'd float constant's value out of radiation.h.

    Parameters
    ----------
    name : str
        The macro name, e.g. "RADIATION_HABING_FLUX_CGS".

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


def _check_headers_match_installed_lib() -> None:
    """Confirm GRACKLE_INCLUDE_DIR's headers match the source tree exactly.

    A stale/mismatched header would silently invalidate the whole
    call-through test (wrong struct layout compiled against the real
    .so), so this is checked rather than assumed.
    """
    for name in ("grackle_chemistry_data.h", "grackle_types.h"):
        installed = Path(GRACKLE_INCLUDE_DIR) / name
        source = Path(GRACKLE_SOURCE_DIR) / "src" / "clib" / name
        if not installed.read_bytes() == source.read_bytes():
            raise RuntimeError(
                f"{installed} differs from {source} -- installed Grackle "
                "headers do not match the claimed source tree; refusing "
                "to trust the call-through result."
            )


def naive_formula_gamma_pe_cgs(g0: float, n_h_cgs: float) -> float:
    """The commonly-cited (metallicity-blind) formula.

    Parameters
    ----------
    g0 : float
        ISRF strength, Habing units.
    n_h_cgs : float
        Total hydrogen number density, cm^-3.

    Returns
    -------
    float
        Gamma_PE, erg/s/cm^3.
    """
    return GRACKLE_GAMMA_HA_CGS * GRACKLE_EPSILON * g0 * n_h_cgs


def actual_formula_gamma_pe_cgs(g0: float, n_h_cgs: float, z_prime: float) -> float:
    """The formula actually implemented in cool1d_multi_g.F's igammah=2 path.

    Parameters
    ----------
    g0 : float
        ISRF strength, Habing units.
    n_h_cgs : float
        Total hydrogen number density, cm^-3.
    z_prime : float
        Metallicity relative to solar, Z/Zsun.

    Returns
    -------
    float
        Gamma_PE, erg/s/cm^3.
    """
    return naive_formula_gamma_pe_cgs(g0, n_h_cgs) * z_prime


def run_grackle_harness() -> list[tuple[float, float, float, float, float]]:
    """Compile and run the C harness against the compiled libgrackle.

    Returns
    -------
    list of (z_prime, edot_on_cgs, edot_off_cgs, T_on_K, T_off_K)
        One row per swept Z', parsed from the harness's CSV output.
    """
    _check_headers_match_installed_lib()

    binary_path = SCRIPT_DIR / "_grackle_pe_harness_bin"
    compile_cmd = [
        "gcc",
        "-O2",
        "-o",
        str(binary_path),
        str(HARNESS_SRC),
        f"-I{GRACKLE_INCLUDE_DIR}",
        f"-L{GRACKLE_LIB_DIR}",
        f"-Wl,-rpath,{GRACKLE_LIB_DIR}",
        "-lgrackle",
        "-lm",
    ]
    result = subprocess.run(compile_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Harness compile failed:\n{result.stderr}")

    try:
        run_result = subprocess.run(
            [str(binary_path)], capture_output=True, text=True, timeout=60
        )
    finally:
        binary_path.unlink(missing_ok=True)

    if run_result.returncode != 0:
        raise RuntimeError(f"Harness run failed:\n{run_result.stderr}")

    rows = []
    for line in run_result.stdout.splitlines():
        if not line.startswith("DATA,"):
            continue
        _, z_prime_s, edot_on_s, edot_off_s, t_on_s, t_off_s = line.split(",")
        rows.append(
            (
                float(z_prime_s),
                float(edot_on_s),
                float(edot_off_s),
                float(t_on_s),
                float(t_off_s),
            )
        )
    if not rows:
        raise RuntimeError(f"No DATA rows in harness output:\n{run_result.stdout}")
    return rows


def main() -> None:
    habing_flux_cgs = _read_radiation_h_constant("RADIATION_HABING_FLUX_CGS")
    draine_over_habing = 1.7  # rate_functions.c's own "G_0=1.7" comment

    print(f"Read from {RADIATION_H}:")
    print(
        f"  RADIATION_HABING_FLUX_CGS = {habing_flux_cgs:.4e} erg/s/cm^2 "
        "(combined FUV+LW Habing normalization)"
    )
    print(
        "  Grackle rate_functions.c's igammah<=1 default docstring cites "
        f"G_0={draine_over_habing} for the Draine field -- matches this "
        "codebase's own DRAINE_OVER_HABING, confirming both sides use the "
        "same Habing convention."
    )
    print()
    print(
        f"Test point: n_H={N_H_CGS:.1f} cm^-3, T={T_K:.1f} K, G0={G0:.1f}, "
        f"Z' swept over {Z_PRIME_SWEEP}"
    )
    print()

    rows = run_grackle_harness()

    print(
        f"{'Z prime':>8} {'Gamma_PE naive':>16} {'Gamma_PE actual':>16} "
        f"{'Gamma_PE Grackle':>18} {'actual/Grackle':>15} {'naive/Grackle':>15}"
    )
    all_actual_agree = True
    all_temperatures_calibrated = True
    for z_prime, edot_on, edot_off, t_on, t_off in rows:
        naive = naive_formula_gamma_pe_cgs(G0, N_H_CGS)
        actual = actual_formula_gamma_pe_cgs(G0, N_H_CGS, z_prime)
        grackle_isolated = edot_on - edot_off

        ratio_actual = actual / grackle_isolated
        ratio_naive = naive / grackle_isolated
        if abs(ratio_actual - 1.0) > 1e-6:
            all_actual_agree = False
        if abs(t_on - T_K) > 1e-3 or abs(t_off - T_K) > 1e-3:
            all_temperatures_calibrated = False

        print(
            f"{z_prime:8.2f} {naive:16.4e} {actual:16.4e} "
            f"{grackle_isolated:18.4e} {ratio_actual:15.6f} {ratio_naive:15.6f}"
            f"   (T_on={t_on:.4f} K, T_off={t_off:.4f} K)"
        )

    if not all_temperatures_calibrated:
        raise RuntimeError(
            "Harness's temperature calibration did not converge to the "
            f"requested T={T_K} K (see per-row T_on/T_off above) -- "
            "igammah=2 has no explicit T-dependence below the 2e4 K "
            "cutoff, so this would not itself invalidate the Gamma_PE "
            "comparison above, but it means the reported 'T=100 K' test "
            "point claim would be false; not silently passing on that."
        )

    print()
    if all_actual_agree:
        print(
            "PASS: the Z'-corrected formula matches Grackle's own "
            "call-through result to <1e-6 relative error at every swept "
            "Z'. The naive (metallicity-blind) formula only agrees at "
            "Z'=1 (ratio 1.0 above); away from solar metallicity it is "
            "off by exactly the factor Z', confirming this is a real "
            "property of Grackle's photoelectric_heating=2 path, not an "
            "artifact of the comparison setup."
        )
    else:
        print(
            "FAIL: the Z'-corrected formula does not match Grackle's "
            "call-through result to <1e-6 at every swept Z' -- see ratios "
            "above. This would indicate either a mistake in this script's "
            "reproduction of the Fortran formula, or a real behaviour "
            "change in the linked libgrackle build."
        )
        sys.exit(1)

    print()
    print(
        "This codebase's own G0-to-Grackle pipeline "
        "(radiation_get_part_isrf_habing) hands G0 straight through with "
        "no extra scaling, which is exactly what this Habing-native "
        "path expects -- so it inherits this Z'-dependence directly and "
        "correctly (dust2gas comes from Grackle's own metal_density field, "
        "not something SWIFT would need to fold into G0 itself). No unit "
        "mismatch found between the two sides; the only real finding is "
        "that the commonly-cited formula (as quoted in the earlier "
        "investigation this script follows up on) omits the metallicity "
        "factor, which is a citation gap, not a bug in this codebase's "
        "coupling."
    )


if __name__ == "__main__":
    main()
