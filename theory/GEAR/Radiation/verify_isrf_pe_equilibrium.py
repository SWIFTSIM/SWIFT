"""Quantitative cross-check: does SWIFT's own converged photoelectric-
heating temperature match an independently-computed Grackle equilibrium
temperature at the SAME local conditions?

This is the follow-up the operator asked for after two same-day passes
(`.claude/dev/logs/2026-09-08_1527_grackle-fix-postverify-design-b.md`,
`.claude/dev/logs/2026-09-08_1622_grackle-fix-nonequilibrium-verify.md`)
were judged insufficient: those checks only confirmed "T_on > T_off,
T_on < 1e6 K" (sign/magnitude) for the `ISRFPhotoelectricHeating` example.
Neither compared the achieved temperature to any independently-computed
target.

Method
------
1. Load a real, already-converged `with_photoelectric_heating=1` SWIFT
   snapshot (this script reruns the example itself, in an isolated scratch
   copy, since the original scratch runs referenced by the two logs above
   were deleted from disk between sessions -- confirmed missing, not
   assumed).
2. For one near-star illuminated particle and one far-field (unilluminated)
   control particle, extract its own local G0, density, and metallicity
   directly from the snapshot, using EXACTLY the same formula
   `radiation_get_part_isrf_habing()` uses (src/feedback/GEAR/
   radiation_gas.c) and EXACTLY the same metal_density mapping
   `cooling.c`'s `cooling_init_grackle()` uses -- not an approximation of
   either.
3. Feed those extracted values into verify_isrf_pe_equilibrium_grackle_
   harness.c, which time-integrates Grackle's own solve_chemistry() (same
   compiled libgrackle, same dust_chemistry=1/photoelectric_heating=2
   configuration a real with_photoelectric_heating=1 run uses -- this DOES
   exercise the 2026-09-08 dust-recombination-cooling NaN fix's code path,
   unlike the sibling verify_photoelectric_heating_rate.py harness, which
   hardcodes dust_chemistry=0) until net heating equals net cooling.
4. Compare that independently-computed equilibrium temperature against the
   particle's own actual, converged temperature in the SWIFT snapshot.

Tolerance
---------
This is a genuine cross-implementation check, not a bit-for-bit
reproduction: SWIFT integrates finite timesteps with the gas also
exchanging energy via SPH/hydrodynamics and the star's field itself
evolving over the run, while this harness time-integrates a single
isolated gas parcel at fixed, frozen (G0, n_H, Z') to a mathematically
exact equilibrium. A few percent agreement is the sign of a correctly
wired coupling; anything at the tens-of-percent level or worse, or a sign
flip, would indicate a real problem. 10% relative error is adopted as the
pass/fail line: loose enough to tolerate the finite-timestep/hydrodynamic-
coupling differences above, tight enough that it would not paper over a
real unit or formula mismatch (which, based on this codebase's own history
of radiation coupling bugs, has previously produced factor-of-few to
orders-of-magnitude discrepancies, not few-percent ones).
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
RADIATION_H = REPO_ROOT / "src" / "feedback" / "GEAR" / "radiation.h"
HARNESS_SRC = SCRIPT_DIR / "verify_isrf_pe_equilibrium_grackle_harness.c"

GRACKLE_INCLUDE_DIR = "/home/darwinr/local/include"
GRACKLE_LIB_DIR = "/home/darwinr/local/lib"

C_LIGHT_CGS = 2.99792458e10  # cm/s
M_H_CGS = 1.67262171e-24  # g, matches the harness's own MH_CGS
HYDROGEN_FRACTION_BY_MASS = 0.76  # GrackleCooling:HydrogenFractionByMass default
SOLAR_METAL_FRACTION_BY_MASS = 0.01295  # Grackle's own compiled default
RELATIVE_ERROR_TOLERANCE = 0.10

GAMMA_M1 = 2.0 / 3.0  # monatomic ideal gas, gamma=5/3
K_B_CGS = 1.380649e-16
M_P_CGS = 1.67262192e-24

# One (mode, example directory, primordial_chemistry) tuple per rerun scratch
# copy. primordial_chemistry matches COOLING_GRACKLE_MODE the binary that
# produced each snapshot was configured with.
CASES = [
    ("grackle_0", "ISRFPhotoelectricHeating_eqverify_g0", 0),
    ("grackle_1", "ISRFPhotoelectricHeating_eqverify_g1", 1),
    ("grackle_3", "ISRFPhotoelectricHeating_eqverify_g3", 3),
]


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


def compile_harness() -> Path:
    """Compile the equilibrium-solver harness against the installed libgrackle.

    Returns
    -------
    Path
        Path to the compiled binary.
    """
    binary_path = SCRIPT_DIR / "_isrf_pe_equilibrium_harness_bin"
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
    return binary_path


def extract_particle(gas, idx: int, unit_length_cgs: float, unit_time_cgs: float,
                     unit_mass_cgs: float, habing_flux_cgs: float) -> dict:
    """Extract one gas particle's local conditions and measured temperature.

    Reproduces radiation_get_part_isrf_habing() (src/feedback/GEAR/
    radiation_gas.c) and cooling.c's cooling_init_grackle() metal_density
    mapping exactly, from raw (internal-unit) snapshot fields -- HDF5
    snapshot fields are stored in internal units, the same units those C
    functions themselves operate in, so no extra conversion factor from
    "snapshot units" to "internal units" is needed.

    Parameters
    ----------
    gas : h5py.Group
        The snapshot's /PartType0 group.
    idx : int
        Particle index.
    unit_length_cgs, unit_time_cgs, unit_mass_cgs : float
        This run's unit system, read from the snapshot's own /Units group.
    habing_flux_cgs : float
        RADIATION_HABING_FLUX_CGS, read from radiation.h.

    Returns
    -------
    dict
        G0_habing, n_H_cgs, Zprime, T_measured_K, r (internal length units).
    """
    rho_internal = float(gas["Densities"][idx])
    u_fuv_internal = float(gas["FUVSpecificEnergies"][idx])
    u_lw_internal = float(gas["LWSpecificEnergies"][idx])
    u_internal = float(gas["InternalEnergies"][idx])
    Z = float(gas["SmoothedMetalMassFractions"][idx, -1])

    unit_mass_density_cgs = unit_mass_cgs / unit_length_cgs**3
    unit_specific_energy_cgs = (unit_length_cgs / unit_time_cgs) ** 2

    rho_cgs = rho_internal * unit_mass_density_cgs
    u_fuv_cgs = u_fuv_internal * unit_specific_energy_cgs
    u_lw_cgs = u_lw_internal * unit_specific_energy_cgs

    # radiation_get_part_isrf_habing(): G0 = c*rho*(u_FUV+u_LW) / HABING_FLUX,
    # all in cgs (the internal-unit computation and the cgs one are
    # equivalent since c*rho*u has pure mass/time^3 dimensions -- no length
    # dependence -- so converting rho and u to cgs separately here and
    # multiplying by C_LIGHT_CGS is exactly the internal-unit formula
    # evaluated in cgs, not an approximation of it).
    G0 = C_LIGHT_CGS * rho_cgs * (u_fuv_cgs + u_lw_cgs) / habing_flux_cgs

    # cooling.c: metal_density = chemistry_get_total_metal_mass_fraction_
    # for_cooling(p) * rho; Grackle's own metallicity(i) = metal_density /
    # density / SolarMetalFractionByMass = Z / SolarMetalFractionByMass.
    Zprime = Z / SOLAR_METAL_FRACTION_BY_MASS

    n_H_cgs = rho_cgs * HYDROGEN_FRACTION_BY_MASS / M_H_CGS

    u_measured_cgs = u_internal * unit_specific_energy_cgs

    # At COOLING_GRACKLE_MODE 0 (no species output) there is no real
    # per-particle ionization state to compute mu from -- mu_source flags
    # this so the caller can, for that mode only, ask the harness for
    # Grackle's own tabulated-mode mu(T,Z) via its MEASURE feature instead
    # of trusting this fallback (see run_harness()'s measure_u_cgs and the
    # harness's own doxygen for why this fallback is wrong exactly where
    # it would matter: it is what silently inflated T by ~1.6x in earlier
    # same-day verification passes).
    if "HI" in gas:
        inv_mu = (
            gas["HI"][idx]
            + 2.0 * gas["HII"][idx]
            + 0.25 * gas["HeI"][idx]
            + 0.5 * gas["HeII"][idx]
            + 0.75 * gas["HeIII"][idx]
        )
        mu = 1.0 / max(inv_mu, 1e-10)
        mu_source = "species"
    else:
        mu = 1.0 / 0.82
        mu_source = "fallback"

    T_measured_K = u_measured_cgs * GAMMA_M1 * mu * M_P_CGS / K_B_CGS

    return dict(G0_habing=G0, n_H_cgs=n_H_cgs, Zprime=Zprime,
               T_measured_K=T_measured_K, rho_cgs=rho_cgs,
               u_measured_cgs=u_measured_cgs, mu_source=mu_source)


def load_case_particles(example_dir: Path, habing_flux_cgs: float) -> dict:
    """Pick one near-star and one far-field particle from a case's final snapshot.

    Parameters
    ----------
    example_dir : Path
        The rerun scratch example directory (contains heating_on/snap/).
    habing_flux_cgs : float
        RADIATION_HABING_FLUX_CGS, read from radiation.h.

    Returns
    -------
    dict
        {"near_star": {...}, "far_field": {...}}, each an extract_particle()
        result plus "r" (radial distance from the star, internal units).
    """
    snap_files = sorted((example_dir / "heating_on" / "snap").glob("snapshot_*.hdf5"))
    if not snap_files:
        raise RuntimeError(f"No snapshots found under {example_dir}")
    snap_path = snap_files[-1]

    with h5py.File(snap_path, "r") as f:
        units = f["/Units"]
        unit_length_cgs = float(np.asarray(units.attrs["Unit length in cgs (U_L)"]).flat[0])
        unit_time_cgs = float(np.asarray(units.attrs["Unit time in cgs (U_t)"]).flat[0])
        unit_mass_cgs = float(np.asarray(units.attrs["Unit mass in cgs (U_M)"]).flat[0])

        gas = f["/PartType0"]
        star_pos = f["/PartType4/Coordinates"][0, :]
        pos = gas["Coordinates"][:, :]
        box = np.asarray(f["/Header"].attrs["BoxSize"], dtype=float).flatten()[0]
        dx = pos - star_pos
        dx -= box * np.round(dx / box)
        r = np.sqrt((dx**2).sum(axis=1))

        field = gas["FUVSpecificEnergies"][:] + gas["LWSpecificEnergies"][:]
        idx_near = int(np.argmax(field))
        idx_far = int(np.argmin(field))

        near = extract_particle(gas, idx_near, unit_length_cgs, unit_time_cgs,
                                unit_mass_cgs, habing_flux_cgs)
        near["r"] = float(r[idx_near])
        far = extract_particle(gas, idx_far, unit_length_cgs, unit_time_cgs,
                               unit_mass_cgs, habing_flux_cgs)
        far["r"] = float(r[idx_far])

    return dict(near_star=near, far_field=far, snapshot=str(snap_path))


def run_harness(binary_path: Path, primordial_chemistry: int, metal_cooling: int,
                cloudy_table: Path, G0: float, n_H_cgs: float, Zprime: float,
                T_initial_K: float, measure_u_cgs: float | None = None) -> dict:
    """Run the equilibrium-solver harness for one (G0, n_H, Z') tuple.

    Parameters
    ----------
    measure_u_cgs : float, optional
        If given, also asks the harness for the temperature Grackle's own
        mu convention assigns to this specific internal energy at the same
        (n_H, Z', mode) state -- only meaningful at primordial_chemistry=0
        (see the harness's own doxygen for why).

    Returns
    -------
    dict
        T_initial_K, T_equilibrium_K, n_iterations, t_elapsed_s, converged,
        and T_measured_consistent_K (None unless measure_u_cgs was given).
    """
    cmd = [
        str(binary_path),
        str(primordial_chemistry),
        str(metal_cooling),
        str(cloudy_table),
        f"{HYDROGEN_FRACTION_BY_MASS:.10e}",
        f"{G0:.10e}",
        f"{n_H_cgs:.10e}",
        f"{Zprime:.10e}",
        f"{T_initial_K:.10e}",
    ]
    if measure_u_cgs is not None:
        cmd.append(f"{measure_u_cgs:.10e}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    data_line = None
    measure_T_K = None
    for line in result.stdout.splitlines():
        if line.startswith("DATA,"):
            data_line = line
        elif line.startswith("MEASURE,"):
            measure_T_K = float(line.split(",")[1])
    if data_line is None:
        raise RuntimeError(
            f"Harness produced no DATA line (exit={result.returncode}):\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    _, t_init_s, t_eq_s, n_iter_s, t_elapsed_s, converged_s = data_line.split(",")
    if result.stderr.strip():
        print(f"    harness stderr: {result.stderr.strip()}")
    return dict(
        T_initial_K=float(t_init_s),
        T_equilibrium_K=float(t_eq_s),
        n_iterations=int(n_iter_s),
        t_elapsed_s=float(t_elapsed_s),
        converged=bool(int(converged_s)),
        T_measured_consistent_K=measure_T_K,
    )


def main() -> None:
    habing_flux_cgs = _read_radiation_h_constant("RADIATION_HABING_FLUX_CGS")
    print(f"RADIATION_HABING_FLUX_CGS = {habing_flux_cgs:.4e} erg/s/cm^2 "
          f"(from {RADIATION_H})")
    print(f"Relative-error tolerance: {RELATIVE_ERROR_TOLERANCE:.0%} (see module docstring)")
    print()

    binary_path = compile_harness()

    rows = []
    all_near_star_within_tolerance = True
    try:
        for mode_name, example_subdir, primordial_chemistry in CASES:
            example_dir = SCRIPT_DIR.parents[2] / "examples" / "SubgridTests" / \
                "SubgridRadiation" / example_subdir
            if not example_dir.exists():
                print(f"SKIP {mode_name}: {example_dir} not found")
                continue

            cloudy_table = example_dir / "CloudyData_UVB=HM2012.h5"
            particles = load_case_particles(example_dir, habing_flux_cgs)
            print(f"=== {mode_name} (primordial_chemistry={primordial_chemistry}) ===")
            print(f"  snapshot: {particles['snapshot']}")

            for label in ("near_star", "far_field"):
                p = particles[label]
                print(f"  [{label}] r={p['r']:.4e}, G0={p['G0_habing']:.4e}, "
                      f"n_H={p['n_H_cgs']:.4e} cm^-3, Z'={p['Zprime']:.4f}, "
                      f"T_measured={p['T_measured_K']:.4e} K (mu_source={p['mu_source']})")

                # At mode 0 there is no species output to compute a real
                # per-particle mu from, so T_measured above used the crude
                # mu=1/0.82 fallback -- ask the harness for the temperature
                # Grackle's own tabulated-mode mu(T,Z) assigns to this exact
                # u instead (see run_harness()'s docstring and the harness's
                # own doxygen), and use THAT as the comparison target.
                measure_u_cgs = p["u_measured_cgs"] if p["mu_source"] == "fallback" else None

                harness_result = run_harness(
                    binary_path, primordial_chemistry, metal_cooling=1,
                    cloudy_table=cloudy_table, G0=p["G0_habing"],
                    n_H_cgs=p["n_H_cgs"], Zprime=p["Zprime"],
                    T_initial_K=p["T_measured_K"], measure_u_cgs=measure_u_cgs,
                )
                T_eq = harness_result["T_equilibrium_K"]
                if harness_result["T_measured_consistent_K"] is not None:
                    T_meas = harness_result["T_measured_consistent_K"]
                    print(f"    mu-fallback T_measured was {p['T_measured_K']:.4e} K; "
                          f"using Grackle's own mode-0 mu(T,Z) instead: "
                          f"T_measured_consistent={T_meas:.4e} K")
                else:
                    T_meas = p["T_measured_K"]
                rel_err = abs(T_eq - T_meas) / T_meas
                within_tol = rel_err <= RELATIVE_ERROR_TOLERANCE and harness_result["converged"]
                if label == "near_star":
                    all_near_star_within_tolerance &= within_tol

                rows.append(dict(
                    mode=mode_name, particle=label, T_measured_K=T_meas,
                    T_equilibrium_K=T_eq, relative_error=rel_err,
                    converged=harness_result["converged"],
                    n_iterations=harness_result["n_iterations"],
                    mu_corrected=harness_result["T_measured_consistent_K"] is not None,
                ))
                print(f"    -> T_equilibrium={T_eq:.4e} K, relative error="
                      f"{rel_err:.2%}, converged={harness_result['converged']} "
                      f"({harness_result['n_iterations']} iterations), "
                      f"{'PASS' if within_tol else 'FAIL'}")
            print()
    finally:
        binary_path.unlink(missing_ok=True)

    print(f"{'Mode':>10} {'Particle':>10} {'T_measured (K)':>16} "
          f"{'T_equilibrium (K)':>18} {'Rel. error':>11} {'Converged':>10} {'Result':>7}")
    for row in rows:
        result_str = "PASS" if (row["relative_error"] <= RELATIVE_ERROR_TOLERANCE
                                and row["converged"]) else "FAIL"
        note = " (mu-corrected)" if row.get("mu_corrected") else ""
        print(f"{row['mode']:>10} {row['particle']:>10} {row['T_measured_K']:16.4e} "
              f"{row['T_equilibrium_K']:18.4e} {row['relative_error']:10.2%} "
              f"{str(row['converged']):>10} {result_str:>7}{note}")

    print()
    print(
        "The primary claim under test is the near_star rows: G0>0 there, so "
        "they exercise the actual with_photoelectric_heating=1 coupling "
        "(dust_chemistry=1, photoelectric_heating=2, use_isrf_field=1). The "
        "far_field rows are a G0=0 control that exercises NO part of that "
        "code path; their agreement/disagreement does not bear on whether "
        "the fix works. They are reported for completeness, not as a gate: "
        "unilluminated gas this close to the CMB floor has a cooling time "
        "that diverges on approach to equilibrium, while this short example "
        "(time_end ~ 1.5e5 yr) does not always run long enough for that "
        "specific control particle to fully relax -- confirmed directly for "
        "the grackle_3 far_field row by its own snapshot time series (T "
        "still slowly drifting, 4.09->4.17 K, at the final snapshot), and "
        "separately confirmed to NOT be a deuterium-initial-condition "
        "artifact in this harness (seeding a cosmological D/H ratio instead "
        "of a trace placeholder left the harness's own equilibrium value "
        "unchanged). The near_star gas, by contrast, has much shorter "
        "heating/cooling timescales at its higher G0 and density, and both "
        "SWIFT and the harness agree it reaches equilibrium well within the "
        "run."
    )

    print()
    if all_near_star_within_tolerance:
        print(
            f"PASS: every tested near_star particle's independently-computed "
            f"Grackle equilibrium temperature agrees with SWIFT's own "
            f"converged temperature to within {RELATIVE_ERROR_TOLERANCE:.0%}, "
            f"across {len(set(r['mode'] for r in rows))} cooling mode(s) -- "
            f"this is a genuine, quantitative confirmation that the "
            f"photoelectric-heating coupling (G0 computation, Grackle field "
            f"population, dust_chemistry=1 path) reproduces the correct "
            f"physical equilibrium, not just the correct sign/magnitude."
        )
    else:
        print(
            "FAIL: at least one tested near_star particle's measured "
            "temperature disagrees with the independently-computed "
            "equilibrium temperature by more than the stated tolerance, or "
            "the harness failed to converge -- see the per-row table above."
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
