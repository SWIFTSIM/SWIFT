#!/usr/bin/env python3
"""Check the H2 formation rate of the geajay chemistry in a dark, uniform box.

With no star there is no light. H2 forms on dust and through H-, and nothing
else acts, so the rate at which the H2 abundance x2 = n_H2 / n_H grows has an
exact answer from the gas state alone. Per hydrogen nucleus,

    dx2/dt = R_d(T, T_d) D n_HI                              (dust)
           + k8 n_H- n_HI / n_H                              (H-)
           + k10 n_H2+ n_HI / n_H + k22 n_HI^3 / n_H         (H2+, three body)

with n_H- in steady state,

    n_H- = k7 n_e n_HI / (k8 n_HI + k14 n_e + k15 n_HI + (k16 + k17) n_HII + k19 n_H2+).

The coefficients are the fits that the library names in rates.c: k7 radiative
attachment (Stancil, Lepp and Dalgarno 1998), k8 associative detachment (Kreckel et
al. 2010), k16 mutual neutralization (Croft et al. 1999) and the dust term R_d
(Omukai 2000). k10, k14, k15, k17, k19 and k22 are the standard fits of Abel et al.
(1997) and Glover and Abel (2008); rates.c gives them no more precise citation. In the
H- arm k8 n_HI is 99.93 per cent of the H- loss rate, so these six fits move the result
by less than the bar. The dust temperature T_d is the root of the grain energy balance
(heating by the field and by gas-grain collisions against emission to the CMB), solved
here with the library's grain opacity and gas-grain coupling. D is the dust-to-gas
ratio of the linear model, D = D_ref Z / Z_sun.

The reference never calls SWIFT or the library. It reads the species of each
snapshot and writes the rate of that state from the formulas above. A mistake that
the library and this file share in a published fit is not seen.

Metric. For each snapshot interval and each particle the measured mean rate is
F_meas = (x2_b - x2_a) / dt. The reference is the trapezoid mean of the rate
of the two snapshot states. The statistic is the median over the particles of
ln(F_meas / F_ref) per interval.

Bar. Per interval it is the sum of terms, each computed from the data:

* the trapezoid error, dt^2 |F''| / (12 F), with F'' from the divided differences of the
  reference rate over the snapshot times that touch the interval (they are not evenly
  spaced);
* the fast species: the relaxation time of H- and of H2+ over the interval length, because
  the reference takes them in balance;
* float32 rounding of x2: 2^-24 (x2_a + x2_b) / |x2_b - x2_a|;
* the temperature and density conversion: the check uses CODATA k_B and m_H, the library uses
  its own (phys_constants.h). That shifts T by 1.1e-6 and n by 1.3e-7. T comes from u and the
  species, which are float32: three roundings of 2^-24. The rate follows T with its local
  slope |d ln F / d ln T|. T uses the library's particle count (no deuterium, metals at mass
  16) and heat capacity (3/2 k_B per particle, 5/2 k_B per H2);
* the solver: the accurate solver keeps the local error of each step below
  atol + rtol |y| per species (error weights, bdf_network.h). The errors of the steps add, so
  the term is n_steps (atol + rtol x2) / |x2_b - x2_a|, with n_steps the SWIFT steps of the
  interval. The solver may take several internal steps per SWIFT step, which the check
  cannot see: the term is an allowance of this form, not a proven bound;
* the dust temperature: |d ln R_d / d ln T_d| times the shift of T_d for a CMB
  temperature that is 0.005 K higher.

The species of a snapshot are not drifted: they hold the state of the last step that
ended at or before the snapshot time. The elapsed time between two snapshots is the
difference of those step ends, read from timesteps.txt.

The first snapshot holds the start values, not a state that the library wrote. The
reference renormalises its species to the sum of the second snapshot and takes its H2+ from
the balance of the second snapshot, scaled with n_HII. The first interval is gated like the
others.

A particle with F_meas <= 0 fails the interval: a flipped sign or a missing channel has no
ln. So does a particle whose reference rate is 0. A value that is not finite fails before
any comparison.

What a PASS proves: in every interval the median measured formation rate
equals the rate of the published channel formulas at the gas state, within the
bar. What it does not prove: that the published fit is the right physics, any
error that appears only in the light (H- photodetachment, the dust temperature
in a transported field), a wrong k8 (H- turns into H2 with a share of 0.999 in
the H- arm, so k8 is not tested), and three-body formation (needs n_H above
1e8 cm^-3).

Exit status: 0 pass, 1 a failed bar, 2 a missing file or field, 3 a value that
is not finite, 4 a void run (a star, light, a solver other than 1, or a
temperature that moved).

The --plant-* options change the data as a known error would, to show that the
check fails when it should. They are offline: they act on the arrays read from
the snapshots, not on SWIFT.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
from numpy.typing import NDArray

Arr = NDArray[np.float64]

KB = 1.380649e-16
MH = 1.67262192369e-24
SIGMA_SB = 5.670374419e-5
T_EV_K = 11604.5
T_CMB = 2.725
T_CMB_SHIFT = 0.005
# Constants of the dust model of the library (rates.c, dust.c, dust_coupling.h)
FGR = 0.009387
D_REF_DEFAULT = 0.009387
Z_SOLAR = 0.01295
HABING_FLUX = 1.6e-3
C_LIGHT = 2.99792458e10
SOLVER_REL_TOL = 1e-6
SOLVER_ABS_TOL = 1e-18
# Constants of the library (phys_constants.h, kernel_base.h)
KB_LIB = 1.3806504e-16
MH_LIB = 1.67262171e-24
MU_METAL = 16.0
FLOAT32_EPS = 2.0**-24

SPECIES = ["HI", "HII", "HM", "H2I", "H2II", "HeI", "HeII", "HeIII", "DI", "DII", "HDI"]

EXIT_PASS, EXIT_FAIL, EXIT_MISSING, EXIT_NONFINITE, EXIT_VOID = 0, 1, 2, 3, 4


def k7_rate(t: Arr) -> Arr:
    """Return k7, HI + e- -> H- + photon [cm^3 s^-1], Stancil et al. (1998)."""
    return 3.0e-16 * (t / 3.0e2) ** 0.95 * np.exp(-t / 9.32e3)


def k8_rate(t: Arr) -> Arr:
    """Return k8, H- + HI -> H2 + e- [cm^3 s^-1], Kreckel et al. (2010)."""
    return (
        1.35e-9
        * (t**9.8493e-2 + 3.2852e-1 * t**5.5610e-1 + 2.771e-7 * t**2.1826)
        / (1.0 + 6.191e-3 * t**1.0461 + 8.9712e-11 * t**3.0424 + 3.2576e-14 * t**3.7741)
    )


def k14_rate(t: Arr) -> Arr:
    """Return k14, H- + e- -> HI + 2e- [cm^3 s^-1], Abel et al. (1997); zero below 0.04 eV."""
    te = t / T_EV_K
    ln = np.log(np.maximum(te, 1e-30))
    c = [
        -18.01849334273,
        2.360852208681,
        -0.2827443061704,
        0.01623316639567,
        -0.03365012031362999,
        0.01178329782711,
        -0.001656194699504,
        0.0001068275202678,
        -2.631285809207e-6,
    ]
    return np.where(te <= 0.04, 0.0, np.exp(sum(ci * ln**i for i, ci in enumerate(c))))


def k15_rate(t: Arr) -> Arr:
    """Return k15, H- + HI -> 2HI + e- [cm^3 s^-1], Abel et al. (1997)."""
    te = t / T_EV_K
    ln = np.log(np.maximum(te, 1e-30))
    c = [
        -20.37260896533324,
        1.139449335841631,
        -0.1421013521554148,
        0.00846445538663,
        -0.0014327641212992,
        0.0002012250284791,
        0.0000866396324309,
        -0.00002585009680264,
        2.4555011970392e-6,
        -8.06838246118e-8,
    ]
    return np.where(
        te <= 0.1,
        2.56e-9 * te**1.78186,
        np.exp(sum(ci * ln**i for i, ci in enumerate(c))),
    )


def k16_rate(t: Arr) -> Arr:
    """Return k16, H- + HII -> 2HI [cm^3 s^-1], Croft et al. (1999)."""
    return 2.4e-6 * (1.0 + t / 2.0e4) / np.sqrt(t)


def k17_rate(t: Arr) -> Arr:
    """Return k17, H- + HII -> H2+ + e- [cm^3 s^-1], Glover and Abel (2008)."""
    return np.where(
        t > 1.0e4,
        4.0e-4 * t**-1.4 * np.exp(-15100.0 / t),
        1.0e-8 * t**-0.4,
    )


def k19_rate(t: Arr) -> Arr:
    """Return k19, H- + H2+ -> H2 + HI [cm^3 s^-1], Glover and Abel (2008)."""
    return 5.0e-7 * np.sqrt(100.0 / t)


def k22_rate(t: Arr) -> Arr:
    """Return k22, 3HI -> H2 + HI [cm^6 s^-1], Glover and Abel (2008)."""
    return np.where(
        t <= 300.0, 1.3e-32 * (t / 300.0) ** -0.38, 1.3e-32 * (t / 300.0) ** -1.0
    )


def h2dust_rate(t: Arr, t_dust: Arr) -> Arr:
    """Return the H2-on-dust coefficient per unit dust-to-gas ratio [cm^3 s^-1], Omukai (2000)."""
    step = 1.0 / (1.0 + np.exp(7.5e2 * (1.0 / 75.0 - 1.0 / t_dust)))
    sticking = 1.0 / (1.0 + 4.0e-2 * np.sqrt(t + t_dust) + 2.0e-3 * t + 8.0e-6 * t**2)
    return 6.0e-17 / FGR * np.sqrt(t / 300.0) * step * sticking


def kappa_grain(t_dust: Arr) -> Arr:
    """Return the grain opacity [cm^2 g^-1]."""
    return np.where(
        t_dust < 200.0,
        4.0e-4 * t_dust**2,
        np.where(t_dust <= 1500.0, 16.0, 16.0 * (t_dust / 1500.0) ** -12),
    )


def dust_temperature(
    t: Arr, n_h: Arr, g0: Arr, d_ref: float, t_cmb: float = T_CMB
) -> Arr:
    """Return the root of the grain energy balance [K], by bisection between T_CMB and T.

    G_0 3.9e-24 / (m_H f_gr) + 4 sigma kappa (T_cmb^4 - T_d^4)
    + gasgr n_H (T - T_d) = 0, with the gas-grain coupling of Hollenbach and McKee
    (1989), 1.2e-31 (T / 1000)^0.5 (1 - 0.8 exp(-75 / T)) n^2 (T - T_d) erg cm^-3 s^-1.
    """
    gasgr = (
        d_ref * 1.2e-31 * 1e3**-0.5 / FGR * np.sqrt(t) * (1.0 - 0.8 * np.exp(-75.0 / t))
    ) / MH
    gamma = g0 * 3.9e-24 / MH / FGR
    lo = np.full_like(t, t_cmb)
    hi = np.maximum(t, t_cmb)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        bal = (
            gamma
            + 4.0 * SIGMA_SB * kappa_grain(mid) * (t_cmb**4 - mid**4)
            + gasgr * n_h * (t - mid)
        )
        lo = np.where(bal > 0.0, mid, lo)
        hi = np.where(bal > 0.0, hi, mid)
    return 0.5 * (lo + hi)


def read_snapshot(path: Path) -> dict[str, Arr]:
    """Read the gas fields of one snapshot in cgs, sorted by particle ID."""
    with h5py.File(path, "r") as f:
        um = float(f["Units"].attrs["Unit mass in cgs (U_M)"][0])
        ul = float(f["Units"].attrs["Unit length in cgs (U_L)"][0])
        ut = float(f["Units"].attrs["Unit time in cgs (U_t)"][0])
        g = f["PartType0"]
        order = np.argsort(g["ParticleIDs"][:])
        out: dict[str, Arr] = {"id": g["ParticleIDs"][:][order].astype(float)}
        out["rho"] = g["Densities"][:][order].astype(float) * um / ul**3
        out["u"] = g["InternalEnergies"][:][order].astype(float) * (ul / ut) ** 2
        for k in SPECIES:
            out[k] = g[k][:][order].astype(float)
        out["u_lw"] = g["LWSpecificEnergies"][:][order].astype(float) * (ul / ut) ** 2
        out["u_pe"] = g["PESpecificEnergies"][:][order].astype(float) * (ul / ut) ** 2
        out["z"] = g["SmoothedMetalMassFractions"][:][order][:, -1].astype(float)
        out["time"] = np.array([float(f["Header"].attrs["Time"][0]) * ut])
        out["time_code"] = np.array([float(f["Header"].attrs["Time"][0])])
        n_star = int(f["Header"].attrs["NumPart_Total"][4])
        out["n_star"] = np.array([float(n_star)])
        par = f["Parameters"].attrs
        solver = par.get("GeajayCooling:solver", None)
        out["solver"] = np.array([float(solver) if solver is not None else np.nan])
        for key, name, default in (
            ("GeajayCooling:accurate_relative_tolerance", "rtol", SOLVER_REL_TOL),
            ("GeajayCooling:accurate_absolute_tolerance", "atol", SOLVER_ABS_TOL),
        ):
            raw = par.get(key, None)
            value = float(raw) if raw is not None else -1.0
            out[name] = np.array([value if value > 0.0 else default])
    return out


def hydrogen(s: dict[str, Arr]) -> dict[str, Arr]:
    """Return the number densities and the temperature of a snapshot state.

    The temperature follows the library: the particles that carry the energy are H, H+, H-,
    H2, H2+, He and the electrons, plus the metals at mass 16; deuterium species do not
    count. H2 has the heat capacity 5/2 k_B per molecule, everything else 3/2 k_B.
    """
    rho = s["rho"]
    n_elec = (
        s["HII"]
        - s["HM"]
        + (s["HeII"] + 2.0 * s["HeIII"]) / 4.0
        + s["H2II"] / 2.0
        + s["DII"] / 2.0
    )
    h2 = 0.5 * s["H2I"]
    n_tot = (
        (s["HeI"] + s["HeII"] + s["HeIII"]) / 4.0
        + s["HI"]
        + s["HII"]
        + s["HM"]
        + h2
        + 0.5 * s["H2II"]
        + n_elec
        + s["z"] / MU_METAL
    )
    heat_capacity = (2.5 * h2 + 1.5 * (n_tot - h2)) / n_tot
    temp = s["u"] * MH / (n_tot * KB * heat_capacity)
    h_sum = s["HI"] + s["HII"] + s["HM"] + s["H2I"] + s["H2II"]
    return {
        "n_hi": s["HI"] * rho / MH,
        "n_hii": s["HII"] * rho / MH,
        "n_h2ii": s["H2II"] * rho / (2.0 * MH),
        "n_e": np.maximum(n_elec * rho / MH, 0.0),
        "n_h": h_sum * rho / MH,
        "x2": 0.5 * s["H2I"] / h_sum,
        "temp": temp,
    }


def hminus_loss(h: dict[str, Arr], t: Arr) -> Arr:
    """Return the destruction rate of H- [s^-1]."""
    return (
        k8_rate(t) * h["n_hi"]
        + k14_rate(t) * h["n_e"]
        + k15_rate(t) * h["n_hi"]
        + (k16_rate(t) + k17_rate(t)) * h["n_hii"]
        + k19_rate(t) * h["n_h2ii"]
    )


def formation_rate(
    s: dict[str, Arr], h: dict[str, Arr], t: Arr, d_ref: float, t_cmb: float = T_CMB
) -> Arr:
    """Return the H2 formation rate per hydrogen nucleus [s^-1] of the state, all channels."""
    g0 = (
        C_LIGHT
        * s["rho"]
        * (np.maximum(s["u_pe"], 0.0) + np.maximum(s["u_lw"], 0.0))
        / HABING_FLUX
    )
    t_dust = dust_temperature(t, h["n_h"], g0, d_ref, t_cmb)
    dust_to_gas = d_ref * s["z"] / Z_SOLAR
    n_hm = k7_rate(t) * h["n_e"] * h["n_hi"] / hminus_loss(h, t)
    volume = (
        h2dust_rate(t, t_dust) * dust_to_gas * h["n_hi"] * h["n_h"]
        + k8_rate(t) * n_hm * h["n_hi"]
        + 6.0e-10 * h["n_h2ii"] * h["n_hi"]
        + k22_rate(t) * h["n_hi"] ** 3
    )
    return volume / h["n_h"]


def step_end_times(timesteps: Path) -> Arr:
    """Return the end time of every step of timesteps.txt (code time)."""
    rows = [
        line.split()
        for line in timesteps.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    return np.array([float(r[1]) for r in rows])


def last_update_time(steps: Arr, t_snap: float) -> float:
    """Return the end of the last step at or before a snapshot time (code time).

    A snapshot drifts positions and energies to its time, but the species stay at the
    value of the last step end: they are not drifted. The elapsed chemistry time between
    two snapshots is therefore the difference of these step ends, not of the snapshot times.
    """
    tol = 1e-6 * float(steps[-1])
    return float(steps[steps <= t_snap + tol].max())


def plant_formation_scale(snaps: list[dict[str, Arr]], scale: float) -> None:
    """Scale the H2 increments of every snapshot by scale and lower HI to keep the mass."""
    h2_0 = snaps[0]["H2I"].copy()
    for s in snaps[1:]:
        new = h2_0 + scale * (s["H2I"] - h2_0)
        s["HI"] = s["HI"] - (new - s["H2I"])
        s["H2I"] = new


def plant_hminus_sign(snaps: list[dict[str, Arr]]) -> None:
    """Mirror the H2 increments of every snapshot about the first, as a flipped H- term would."""
    plant_formation_scale(snaps, -1.0)


def start_up_state(first: dict[str, Arr], second: dict[str, Arr]) -> dict[str, Arr]:
    """Return the first snapshot with the species the library holds after its first solve.

    The first snapshot holds the start values. The library renormalises the species to
    the sum it keeps afterwards (hydrogen, helium and deuterium plus the metals), and H2+
    starts at 1e-20 although its balance is reached within a step. The reference takes the
    renormalised species and H2+ from the balance of the second snapshot, scaled with n_HII,
    which is what its production follows.
    """
    state = {k: v.copy() for k, v in first.items()}
    total_first = sum(first[k] for k in SPECIES) + first["z"]
    total_second = sum(second[k] for k in SPECIES) + second["z"]
    for k in SPECIES:
        state[k] = first[k] * total_second / total_first
    state["H2II"] = second["H2II"] * state["HII"] / second["HII"]
    return state


def median_ln(ratio: Arr) -> float:
    """Return the median of ln(ratio) over the positive finite entries."""
    good = ratio[np.isfinite(ratio) & (ratio > 0.0)]
    return float(np.median(np.log(good))) if good.size else float("nan")


def term_trapezoid(f_ref: list[Arr], t_s: list[float], j: int, f_mean: Arr) -> float:
    """Return the relative error of the trapezoid rule over interval j.

    It is dt^2 |F''| / (12 F), with F'' from the divided differences of the reference
    rate over the three snapshot times that touch the interval, which are not evenly spaced.
    """
    dt = t_s[j + 1] - t_s[j]
    curvature = []
    for i in (j - 1, j):
        if i < 0 or i + 2 >= len(f_ref):
            continue
        slope_up = (f_ref[i + 2] - f_ref[i + 1]) / (t_s[i + 2] - t_s[i + 1])
        slope_down = (f_ref[i + 1] - f_ref[i]) / (t_s[i + 1] - t_s[i])
        curvature.append(np.abs(2.0 * (slope_up - slope_down) / (t_s[i + 2] - t_s[i])))
    if not curvature:
        return 0.0
    return float(np.median(dt**2 * np.max(curvature, axis=0) / (12.0 * f_mean)))


def term_fast_species(a: dict[str, Arr], dt: float) -> float:
    """Return the relaxation time of the slowest fast species, H- or H2+, over the interval."""
    t_hminus = 1.0 / hminus_loss(a, a["temp"])
    t_h2ii = 1.0 / (6.0e-10 * a["n_hi"])
    return float(np.median(np.maximum(t_hminus, t_h2ii) / dt))


def term_float32(a: dict[str, Arr], b: dict[str, Arr]) -> float:
    """Return the float32 rounding of the two x2 values over their change."""
    change = np.maximum(np.abs(b["x2"] - a["x2"]), 1e-300)
    return float(np.median(FLOAT32_EPS * (a["x2"] + b["x2"]) / change))


def term_temperature(s: dict[str, Arr], a: dict[str, Arr], d_ref: float) -> float:
    """Return the error of the rate from the temperature and density conversion.

    The check uses CODATA k_B and m_H. The library uses its own values (phys_constants.h),
    which shifts T by (m_H,lib / m_H) (k_B / k_B,lib) - 1 and n by m_H / m_H,lib - 1.
    The temperature comes from u and the species, float32 in the snapshot: three roundings.
    The rate follows T with the local slope d ln F / d ln T.
    """
    slope = np.log(
        formation_rate(s, a, a["temp"] * 1.001, d_ref)
        / formation_rate(s, a, a["temp"] / 1.001, d_ref)
    ) / np.log(1.001**2)
    shift_t = abs((MH_LIB / MH) * (KB / KB_LIB) - 1.0) + 3.0 * FLOAT32_EPS
    shift_n = abs(MH_LIB / MH - 1.0) + FLOAT32_EPS
    return float(np.median(np.abs(slope) * shift_t + shift_n))


def term_solver(
    a: dict[str, Arr], b: dict[str, Arr], s: dict[str, Arr], n_steps: float
) -> float:
    """Return the accumulated error weight of the solver over the increment of x2.

    The accurate solver keeps the local error of each step below atol + rtol |y| per
    species (bdf_network.h, error weights). The errors of the steps add.
    """
    weight = float(s["atol"][0]) + float(s["rtol"][0]) * 0.5 * (
        np.abs(a["x2"]) + np.abs(b["x2"])
    )
    change = np.maximum(np.abs(b["x2"] - a["x2"]), 1e-300)
    return float(np.median(n_steps * weight / change))


def term_dust_temperature(a: dict[str, Arr], t_mid: Arr, d_ref: float) -> float:
    """Return the error of the dust term from a CMB temperature 0.005 K higher."""
    zero = np.zeros_like(a["temp"])
    t_shift = dust_temperature(a["temp"], a["n_h"], zero, d_ref, T_CMB + T_CMB_SHIFT)
    t_base = dust_temperature(a["temp"], a["n_h"], zero, d_ref)
    return float(
        np.median(
            np.abs(np.log(h2dust_rate(t_mid, t_shift) / h2dust_rate(t_mid, t_base)))
        )
    )


def interval_row(
    j: int,
    snaps: list[dict[str, Arr]],
    states: list[dict[str, Arr]],
    f_ref: list[Arr],
    t_upd: list[float],
    steps: Arr,
    t_unit: float,
    d_ref: float,
    arm: str,
) -> dict:
    """Return the statistic, the bar terms and the counts of interval j."""
    a, b = states[j], states[j + 1]
    dt = (t_upd[j + 1] - t_upd[j]) * t_unit
    f_meas = (b["x2"] - a["x2"]) / dt
    f_mean = 0.5 * (f_ref[j] + f_ref[j + 1])
    n_steps = float(
        np.sum(
            (steps > t_upd[j] + 1e-6 * steps[-1])
            & (steps <= t_upd[j + 1] + 1e-6 * steps[-1])
        )
    )
    terms = {
        "trapezoid": term_trapezoid(f_ref, [t * t_unit for t in t_upd], j, f_mean),
        "fast_species": term_fast_species(a, dt),
        "float32": term_float32(a, b),
        "temperature": term_temperature(snaps[j], a, d_ref),
        "solver": term_solver(a, b, snaps[j], n_steps),
        "dust_temperature": (
            term_dust_temperature(a, 0.5 * (a["temp"] + b["temp"]), d_ref)
            if arm == "dust"
            else 0.0
        ),
    }
    return {
        "interval": j,
        "t_mid_Myr": float(0.5 * (t_upd[j] + t_upd[j + 1]) * t_unit / 3.15576e13),
        "dt_s": dt,
        "statistic": median_ln(f_meas / f_mean),
        "n_nonpositive": int(np.sum(~(f_meas > 0.0))),
        "n_zero_reference": int(np.sum(~(f_mean > 0.0))),
        "n_particles": int(f_meas.size),
        "bar": float(sum(terms.values())),
        "terms": terms,
        "f_ref_median": float(np.median(f_mean)),
        "f_meas_median": float(np.median(f_meas)),
    }


def evaluate(snaps: list[dict[str, Arr]], steps: Arr, d_ref: float, arm: str) -> dict:
    """Return the per-interval statistic, bar terms and verdict of the arm."""
    snaps = [start_up_state(snaps[0], snaps[1])] + snaps[1:]
    states = [hydrogen(s) for s in snaps]
    t_unit = float(snaps[1]["time"][0] / snaps[1]["time_code"][0])
    t_upd = [last_update_time(steps, float(s["time_code"][0])) for s in snaps]
    f_ref = [formation_rate(s, h, h["temp"], d_ref) for s, h in zip(snaps, states)]
    rows = [
        interval_row(j, snaps, states, f_ref, t_upd, steps, t_unit, d_ref, arm)
        for j in range(len(snaps) - 1)
    ]
    return {"intervals": rows}


def nonfinite_report(snaps: list[dict[str, Arr]]) -> list[str]:
    """Return one message per non-finite field of any snapshot."""
    out = []
    for i, s in enumerate(snaps):
        for k, v in s.items():
            if v.dtype.kind == "f" and k != "solver" and not np.all(np.isfinite(v)):
                out.append(
                    f"snapshot {i}: {int(np.sum(~np.isfinite(v)))} non-finite values in {k}"
                )
    return out


def void_reasons(snaps: list[dict[str, Arr]], arm: str) -> list[str]:
    """Return the reasons why the run is void for this check."""
    out = []
    for i, s in enumerate(snaps):
        if s["n_star"][0] > 0:
            out.append(f"snapshot {i}: a star is present")
        if np.any(s["u_lw"] != 0.0) or np.any(s["u_pe"] != 0.0):
            out.append(f"snapshot {i}: nonzero LW or PE energy")
        if s["solver"][0] != 1.0:
            out.append(f"snapshot {i}: GeajayCooling:solver is {s['solver'][0]}, not 1")
    temps = [hydrogen(s)["temp"] for s in snaps]
    drift = max(float(np.max(np.abs(t / temps[0] - 1.0))) for t in temps)
    if drift > 0.01:
        out.append(
            f"the temperature moved by {drift:.3g} of its start value (limit 0.01)"
        )
    z = float(np.median(snaps[0]["z"]))
    if arm == "hminus" and z != 0.0:
        out.append(f"the H- arm needs zero metallicity, got {z:.3g}")
    if arm == "dust" and z <= 0.0:
        out.append("the dust arm needs a nonzero metallicity")
    return out


def interval_verdict(r: dict) -> bool:
    """Return True if the interval meets its bar with a positive rate for every particle."""
    return (
        r["n_nonpositive"] == 0
        and r["n_zero_reference"] == 0
        and abs(r["statistic"]) <= r["bar"]
    )


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("snapdir", type=Path, help="directory of the snapshots")
    p.add_argument("--arm", required=True, choices=["hminus", "dust"])
    p.add_argument(
        "--timesteps",
        type=Path,
        default=None,
        help="timesteps.txt (default: next to snap/)",
    )
    p.add_argument(
        "--expected-dust-to-gas",
        type=float,
        default=D_REF_DEFAULT,
        help="reference dust-to-gas ratio D_ref of the model, the library value by default",
    )
    p.add_argument(
        "--plant-formation-scale",
        type=float,
        default=None,
        help="plant: scale the H2 increments",
    )
    p.add_argument(
        "--plant-hm-sign", action="store_true", help="plant: mirror the H2 increments"
    )
    p.add_argument(
        "--plant-nan",
        action="store_true",
        help="plant: one NaN in the density of the last snapshot",
    )
    p.add_argument("--json", type=Path, default=None, help="write the result as JSON")
    o = p.parse_args()

    files = sorted(o.snapdir.glob("snapshot_*.hdf5"))
    if len(files) < 4:
        print(f"MISSING: {len(files)} snapshots in {o.snapdir}, need at least 4")
        return EXIT_MISSING
    timesteps = o.timesteps or (o.snapdir.parent / "timesteps.txt")
    if not timesteps.exists():
        print(f"MISSING: {timesteps}")
        return EXIT_MISSING
    try:
        snaps = [read_snapshot(f) for f in files]
    except (KeyError, OSError) as err:
        print(f"MISSING: {err}")
        return EXIT_MISSING

    if o.plant_nan:
        snaps[-1]["rho"][0] = np.nan
    bad = nonfinite_report(snaps)
    if bad:
        print("NON-FINITE:\n  " + "\n  ".join(bad))
        return EXIT_NONFINITE
    void = void_reasons(snaps, o.arm)
    if void:
        print("VOID:\n  " + "\n  ".join(void))
        return EXIT_VOID

    if o.plant_formation_scale is not None:
        plant_formation_scale(snaps, o.plant_formation_scale)
        print(f"PLANT: H2 increments x {o.plant_formation_scale}")
    if o.plant_hm_sign:
        plant_hminus_sign(snaps)
        print("PLANT: H2 increments mirrored (flipped H- sign)")

    res = evaluate(snaps, step_end_times(timesteps), o.expected_dust_to_gas, o.arm)
    z_over_zsun = float(np.median(snaps[0]["z"]) / Z_SOLAR)
    print(
        f"arm {o.arm}: Z / Z_sun = {z_over_zsun:.4g}, D_ref = {o.expected_dust_to_gas:.6g}"
    )
    print(
        f"{'int':>3} {'t[Myr]':>7} {'median ln(F_meas/F_ref)':>24} {'bar':>10} {'<=0':>5}  verdict"
    )
    status = EXIT_PASS
    for r in res["intervals"]:
        inputs = [r["bar"], r["f_ref_median"], *r["terms"].values()]
        if not all(np.isfinite(inputs)):
            print(f"NON-FINITE bar or reference in interval {r['interval']}")
            return EXIT_NONFINITE
        ok = interval_verdict(r)
        if not ok:
            status = EXIT_FAIL
        print(
            f"{r['interval']:>3} {r['t_mid_Myr']:>7.4f} {r['statistic']:>24.3e} {r['bar']:>10.2e} "
            f"{r['n_nonpositive'] + r['n_zero_reference']:>5}  {'ok' if ok else 'FAIL'}"
        )
    finite = [r for r in res["intervals"] if np.isfinite(r["statistic"])]
    if finite:
        worst = max(finite, key=lambda r: abs(r["statistic"]) / r["bar"])
        ratio = abs(worst["statistic"]) / worst["bar"]
        print(f"largest |statistic| / bar = {ratio:.3g} (interval {worst['interval']})")
    else:
        print("no interval has a positive measured rate for any particle")
    print("RESULT:", "PASS" if status == EXIT_PASS else "FAIL")
    if o.json:
        o.json.write_text(json.dumps({"arm": o.arm, "exit": status, **res}, indent=1))
    return status


if __name__ == "__main__":
    sys.exit(main())
