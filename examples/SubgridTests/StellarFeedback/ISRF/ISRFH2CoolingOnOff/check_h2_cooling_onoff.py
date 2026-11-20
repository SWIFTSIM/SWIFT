#!/usr/bin/env python3
"""Compare the H2 abundance of two legs, line absorption on and off.

The on leg writes the line transmission T_line of every gas particle. In
photodissociation equilibrium, with x = n_H2 / n_H, the formation rate per H
atom F over the star rate k_* is x_on T_line with the switch on and
x_off S(n_H2,off L_J) with the switch off. The two legs hold different gas, so
each leg enters with its own state:
x_off S(x_off n_H,off L_off) = x_on T_line (F_off / F_on) (k_on / k_off).
k_* is proportional to rho u_LW (radiation_gas.c). F is the sum of the library
formation channels of the leg: on dust, h2dust(T, T_dust) d2g n_HI n_H
(rates.c h2dust_rate, element_blocks.h HI_to_H2I); through H-,
k8 n_H- n_HI with n_H- solved from the leg's own electrons (k7, k8, k14, k15,
k16, k17, k19); through H2+, k10 n_H2+ n_HI; three-body, k22 n_HI^3. T_dust is
the root of the library grain balance (dust.c gr_balance) of the leg. At low
metallicity H2 forms mostly through H-, at high metallicity on dust. The
solve calls no SWIFT or library code, but it uses the library rate
coefficients and shielding fit.

What a PASS proves: in each decided bin, the median change of x between the
switch-on and switch-off runs, particle by particle, equals the change
predicted by photodissociation equilibrium from the T_line written by the on
leg, the library shielding fit at the off-leg Jeans column, and the ratio of the
two legs' own library formation rates and rho u_LW, within the bar. What it
does not prove: that T_line, the library formation rates, the shielding fit or
the LW energy are right. It cannot see a fault that lowers u_LW together with
the dissociation rate (the two cancel in the ratio), an error in the absolute
size or sign of a formation channel (only F_off / F_on enters), or an error in
T_line or in the column smaller than the bar, which is set by the run-to-run
noise and by the local slope s of the shielding fit.

Statistic and bar of a density bin. The discrepancy is the median over the bin
of the per-particle residual, the measured minus the derived ln(x_on / x_off).
The legs share their particles, so the error of this median is the bootstrap
standard error of the residual itself (se, particles resampled together). The
bar is 2 se + noise + lag. The noise is |median ln(x_on / x_repeat)| plus twice
its bootstrap standard error, from --repeat, a second run of the on leg. The lag
is the one-step lead of the library temperature: the cooling step integrates
from the state at the start of the step, so the temperature that sets x can
differ from the snapshot temperature by up to the change of one step (library
cooling.c: the energy starts at u + dt du/dt). The gate takes the change of
temperature per step of each particle in each leg from the previous snapshot
(change divided by the number of steps in timesteps.txt, an average) and uses
the median shift of the derived value when both legs use T minus that change.
Without the previous snapshots or timesteps.txt the lag is unknown and the bin
is INDETERMINATE.

The reference legs MUST use the accurate (BDF) solver, GeajayCooling:solver 1.
The cheap solver holds the library self-shielding factor f_loc of the
switch-off leg fixed over each sub-step (geajay primordial_chemistry.c:1033-1035
and :1297-1300). While the H2 falls, that leg lags its equilibrium by
s/(1-s) times the change of ln x over its last sub-step, which is at most
dt/2 (primordial_chemistry.c:1174-1176); the gate has no per-step rate of the
off leg, so it cannot evaluate this one-sided bound. Cheap-solver legs are
therefore reported INDETERMINATE with this reason, never FAIL against the bar.
The solver is read from the Parameters group of every snapshot, printed, and
the gate refuses to decide if the legs differ.

The root search of each particle stops at n_H2 = n_H / 2. A particle is masked
only if it has no light, has no formation (F = 0 in a leg), has no root, has
more than one root in that range, or is not in photodissociation equilibrium;
the counts are printed by reason. Any non-finite input of any particle, and any
negative or non-finite formation rate, is counted and fails the snapshot.

A bin is decided only with at least N_MIN particles and a bar of at most
MAX_BAR. Otherwise it is INDETERMINATE and decides nothing.

Exit status per snapshot: 0 every decided bin agrees and at least one bin is
decided, 1 failure (non-finite value, bad formation rate, or a decided bin
beyond its bar), 2 inconclusive or indeterminate (noise floor unknown, lag
unknown, cheap solver, or no decided bin), 3 nothing left after the masks. The
run exits 1 on any 1 or when every snapshot is 3, 0 only if every snapshot is 0,
else 2.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
from numpy.typing import NDArray

KB = 1.380649e-16
MH = 1.67262192369e-24
GRAV = 6.6743e-8
GAMMA = 5.0 / 3.0
KDISS_DRAINE = 5.85e-11
FLUX_DRAINE = 4.053e-4
C_LIGHT = 2.99792458e10
BINS = {"B1": (1e1, 1e2), "B2": (1e2, 1e3), "B3": (1e3, 1e4)}
X_H = 0.76
# Smallest bin the gate decides. At 30 particles the bootstrap standard error of
# the median is unbiased (0.00514 against 0.00517 measured over 400 subsamples of
# the second fixture, snapshot 6) and scatters by 28 percent; at 10 it scatters by 68 percent.
N_MIN = 30
# A bin whose bar exceeds ln 2 cannot see a factor 2 in the H2 abundance.
MAX_BAR = float(np.log(2.0))
SPECIES_A = {
    "HI": 1.0,
    "HII": 1.0,
    "HM": 1.0,
    "H2I": 2.0,
    "H2II": 2.0,
    "HeI": 4.0,
    "HeII": 4.0,
    "HeIII": 4.0,
    "DI": 2.0,
    "DII": 2.0,
    "HDI": 3.0,
}
# Constants of the geajay build (rates.c, dust.c, dust_coupling.h, radiation.h).
FGR = 0.009387  # GEAJAY_FGR, rates.c
LOCAL_D2G = 0.009387  # local_dust_to_gas_ratio default (the run log prints 0.00939)
Z_SOLAR = 0.01295  # SolarMetalFractionByMass default
HABING_FLUX = 1.6e-3  # RADIATION_HABING_FLUX_CGS [erg s^-1 cm^-2]
SIGMA_SB = 5.670374419e-5
T_CMB = 2.725  # dust T_rad floor: geajay_tcmb at redshift 0
T_EV_K = 11604.5
Y_GRID = 10.0 ** np.arange(-16.0, 3.0001, 0.01)
Arr = NDArray[np.float64]


def wg_alpha(temp: Arr, n_tot: Arr) -> Arr:
    """Return the exponent of the shielding fit, Eq. 7 of Wolcott-Green 2019."""
    tc = np.clip(temp, 1e2, 8e3)
    nc = np.clip(n_tot, 1e-8, 1e7)
    return (0.8711 * np.log10(tc) - 1.928) * np.exp(-0.2856 * np.log10(nc)) + (
        -0.9639 * np.log10(tc) + 3.892
    )


def wg_factor(col: Arr, temp: Arr, n_tot: Arr) -> Arr:
    """Return the shielding factor S at the column col in cm^-2."""
    x = 2e-15 * col
    b5 = 1e-5 * np.sqrt(KB * temp / MH)
    s = 0.965 / (1.0 + x / b5) ** wg_alpha(temp, n_tot) + 0.035 * np.exp(
        -8.5e-4 * np.sqrt(1.0 + x)
    ) / np.sqrt(1.0 + x)
    return np.minimum(s, 1.0)


def h2dust_rate(temp: Arr, t_dust: Arr) -> Arr:
    """Return the H2-on-dust coefficient [cm^3 s^-1], Omukai (2000), rates.c h2dust_rate.

    Excludes the dust-to-gas ratio.
    """
    term1 = 1.0 / (1.0 + np.exp(7.5e2 * (1.0 / 75.0 - 1.0 / t_dust)))
    term2 = 1.0 / (
        1.0 + 4.0e-2 * np.sqrt(temp + t_dust) + 2.0e-3 * temp + 8.0e-6 * temp**2
    )
    return 6.0e-17 / FGR * np.sqrt(temp / 300.0) * term1 * term2


def kappa_grain(t_dust: Arr) -> Arr:
    """Return the grain opacity [cm^2 g^-1] of dust.c geajay_kappa_gr."""
    return np.where(
        t_dust < 200.0,
        4.0e-4 * t_dust**2,
        np.where(t_dust <= 1500.0, 16.0, 16.0 * (t_dust / 1500.0) ** -12),
    )


def dust_temperature(temp: Arr, n_h_lib: Arr, g0: Arr) -> Arr:
    """Return the root of the grain balance of dust.c gr_balance [K].

    Gamma_ISRF + 4 sigma kappa (T_rad^4 - T_d^4) + gasgr n_H (T - T_d) = 0, with
    Gamma_ISRF = G0 3.9e-24 / (m_H f_gr) [erg g^-1 s^-1] and the gas-grain
    coupling gasgr = d2g_local gasGrain_rate(T) / m_H [erg cm^3 g^-1 s^-1 K^-1]
    with gasGrain_rate = 1.2e-31 1000^-0.5 / f_gr sqrt(T) (1 - 0.8 exp(-75/T))
    (dust_coupling.h:186-187, rates.c:1466-1470). Bisection between T_rad and T.
    """
    gasgr = (
        LOCAL_D2G
        * 1.2e-31
        * 1e3**-0.5
        / FGR
        * np.sqrt(temp)
        * (1.0 - 0.8 * np.exp(-75.0 / temp))
    ) / MH
    gamma = g0 * 3.9e-24 / MH / FGR
    lo = np.full_like(temp, T_CMB)
    hi = np.maximum(temp, T_CMB)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        bal = (
            gamma
            + 4.0 * SIGMA_SB * kappa_grain(mid) * (T_CMB**4 - mid**4)
            + gasgr * n_h_lib * (temp - mid)
        )
        lo = np.where(bal > 0.0, mid, lo)
        hi = np.where(bal > 0.0, hi, mid)
    return 0.5 * (lo + hi)


def k8_rate(t: Arr) -> Arr:
    """Return k8, HM + HI -> H2 + e- [cm^3 s^-1] (rates.c k8_rate)."""
    return (
        1.35e-9
        * (t**9.8493e-2 + 3.2852e-1 * t**5.5610e-1 + 2.771e-7 * t**2.1826)
        / (
            1.0
            + 6.191e-3 * t**1.0461
            + 8.9712e-11 * t**3.0424
            + 3.2576e-14 * t**3.7741
        )
    )


def k14_rate(t: Arr) -> Arr:
    """Return k14, HM + e- -> HI + 2e- [cm^3 s^-1] (zero below 0.04 eV, the floor)."""
    ln = np.log(np.maximum(t / T_EV_K, 1e-30))
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
    return np.where(
        t / T_EV_K <= 0.04, 0.0, np.exp(sum(ci * ln**i for i, ci in enumerate(c)))
    )


def k15_rate(t: Arr) -> Arr:
    """Return k15, HM + HI -> 2HI + e- [cm^3 s^-1] (rates.c k15_rate)."""
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


HM_SIGN = 1.0  # harness only: --plant-hm-sign sets -1


def formation_rate(leg: dict[str, Arr], st: dict[str, Arr], temp: Arr) -> Arr:
    """Return the H2 formation rate per H atom F [s^-1] of one leg.

    The sum of the formation channels of the library from the leg's own species,
    density and the temperature temp. Per volume,
    dn_H2/dt = h2dust d2g n_HI n_H + k8 n_H- n_HI + k10 n_H2+ n_HI + k22 n_HI^3,
    with n_H the hydrogen of geajay_rhoH (HI, HII, H2, H2+) and n_H- at its
    steady state from the leg's electrons. F is that rate over n_H = X_H rho / m_H.
    """
    rho = leg["rho"]
    n_hi = leg["HI"] * rho / MH
    n_hii = leg["HII"] * rho / MH
    n_h2p = leg["H2II"] * rho / (2.0 * MH)
    n_e = np.maximum(
        (
            leg["HII"]
            - leg["HM"]
            + (leg["HeII"] + 2.0 * leg["HeIII"]) / 4.0
            + leg["H2II"] / 2.0
            + leg["DII"] / 2.0
        )
        * rho
        / MH,
        0.0,
    )
    n_h_lib = (leg["HI"] + leg["HII"] + leg["H2I"] + leg["H2II"]) * rho / MH
    g0 = (
        C_LIGHT
        * rho
        * (np.maximum(leg["u_pe"], 0.0) + np.maximum(leg["u_lw"], 0.0))
        / HABING_FLUX
    )
    t_dust = dust_temperature(temp, n_h_lib, g0)
    d2g = LOCAL_D2G * leg["z_smooth"] / Z_SOLAR
    dust = h2dust_rate(temp, t_dust) * d2g * n_hi * n_h_lib
    k16 = 2.4e-6 * (1.0 + temp / 2e4) / np.sqrt(temp)
    k17 = np.where(
        temp > 1e4,
        4.0e-4 * temp**-1.4 * np.exp(-15100.0 / temp),
        1.0e-8 * temp**-0.4,
    )
    k19 = 5.0e-7 * np.sqrt(100.0 / temp)
    k7 = 3.0e-16 * (temp / 300.0) ** 0.95 * np.exp(-temp / 9.32e3)
    n_hm = (
        k7
        * n_e
        * n_hi
        / (
            k8_rate(temp) * n_hi
            + k14_rate(temp) * n_e
            + k15_rate(temp) * n_hi
            + (k16 + k17) * n_hii
            + k19 * n_h2p
        )
    )
    k22 = np.where(
        temp <= 300.0,
        1.3e-32 * (temp / 300.0) ** -0.38,
        1.3e-32 * (temp / 300.0) ** -1.0,
    )
    volume_rate = (
        dust
        + HM_SIGN * k8_rate(temp) * n_hm * n_hi
        + 6.0e-10 * n_h2p * n_hi
        + k22 * n_hi**3
    )
    return volume_rate / st["n_h"]


def jeans_length(temp: Arr, rho: Arr, mu: Arr) -> Arr:
    """Return the Jeans length in cm for gas of density rho in g cm^-3."""
    return np.sqrt(GAMMA * np.pi * KB * temp / (GRAV * mu * MH * rho))


def solve_off(
    y_on: Arr, t_line: Arr, temp: Arr, n_tot: Arr, length: Arr, y_max: Arr
) -> tuple[Arr, NDArray[np.int64]]:
    """Solve y f(y L) = y_on T_line for y, the H2 number density in cm^-3.

    The search stops at y_max, the physical bound n_H / 2 of each particle.
    Returns the root, NaN unless the root is unique, and the root count.
    """
    target = np.log(y_on * t_line)
    y_root = np.full(y_on.shape, np.nan)
    count = np.zeros(y_on.shape, dtype=np.int64)
    ly = np.log(Y_GRID)
    for i in range(y_on.size):
        n = int(np.searchsorted(Y_GRID, y_max[i], side="right"))
        s = wg_factor(Y_GRID[:n] * length[i], temp[i], n_tot[i])
        d = ly[:n] + np.log(s) - target[i]
        idx = np.nonzero(np.sign(d[:-1]) * np.sign(d[1:]) <= 0.0)[0]
        count[i] = idx.size
        if idx.size == 1:
            j = idx[0]
            w = d[j] / (d[j] - d[j + 1]) if d[j] != d[j + 1] else 0.0
            y_root[i] = np.exp(ly[j] + w * (ly[j + 1] - ly[j]))
    return y_root, count


def slope(col: Arr, temp: Arr, n_tot: Arr, eps: float = 1e-3) -> Arr:
    """Return s = -d ln S / d ln N of the fit."""
    hi = np.log(wg_factor(col * (1 + eps), temp, n_tot))
    lo = np.log(wg_factor(col * (1 - eps), temp, n_tot))
    return -(hi - lo) / (np.log1p(eps) - np.log1p(-eps))


# GeajayCooling:solver as printed by cooling_print_backend of the geajay build used
# (leg_on/output.log:202 prints "solver = 0" for the parameter value -1).
SOLVER_NAMES = {-1: "cheap", 0: "cheap", 1: "accurate (BDF)"}


def read_leg(path: Path, line_fields: bool) -> dict[str, Arr]:
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
        for k in SPECIES_A:
            out[k] = g[k][:][order].astype(float)
        out["u_lw"] = g["LWSpecificEnergies"][:][order].astype(float) * (ul / ut) ** 2
        out["u_pe"] = g["PESpecificEnergies"][:][order].astype(float) * (ul / ut) ** 2
        out["z_smooth"] = g["SmoothedMetalMassFractions"][:][order][:, -1].astype(float)
        if line_fields:
            out["t_line"] = g["LWLineTransmissions"][:][order].astype(float)
        t0 = float(f["Header"].attrs["Time"][0]) * ut
        births = f["PartType4/BirthTimes"][:] * ut if "PartType4" in f else np.array([])
        out["since_star"] = np.array([t0 - births.min() if births.size else np.nan])
        out["time"] = np.array([t0])
        out["time_raw"] = np.array([float(f["Header"].attrs["Time"][0])])
        raw = f["Parameters"].attrs.get("GeajayCooling:solver", None)
        out["solver"] = np.array([float(raw) if raw is not None else np.nan])
    return out


def pair(leg: dict[str, Arr], ids: Arr) -> dict[str, Arr]:
    """Keep the particles whose ID is in ids; scalars pass through."""
    keep = np.isin(leg["id"], ids)
    return {k: (v[keep] if v.size == keep.size else v) for k, v in leg.items()}


def state(leg: dict[str, Arr]) -> dict[str, Arr]:
    """Return n_tot, mu, temperature and n_H2 from the species fractions."""
    n_part = sum(leg[k] / a for k, a in SPECIES_A.items())
    n_elec = (
        leg["HII"]
        - leg["HM"]
        + (leg["HeII"] + 2 * leg["HeIII"]) / 4.0
        + leg["H2II"] / 2.0
        + leg["DII"] / 2.0
    )
    n_part = n_part + n_elec
    mu = 1.0 / n_part
    n_tot = leg["rho"] * n_part / MH
    temp = (GAMMA - 1.0) * leg["u"] * mu * MH / KB
    n_h2 = 0.5 * leg["H2I"] * leg["rho"] / MH
    return {
        "mu": mu,
        "n_tot": n_tot,
        "temp": temp,
        "n_h2": n_h2,
        "n_h": X_H * leg["rho"] / MH,
    }


def boot_se(v: Arr, rng: np.random.Generator, n: int = 400) -> float:
    """Return the bootstrap standard error of the median."""
    if v.size < 2:
        return float("nan")
    m = [np.median(rng.choice(v, v.size)) for _ in range(n)]
    return float(np.std(m))


def finite_counts(legs: dict[str, dict[str, Arr] | None], n: int) -> dict[str, int]:
    """Return the non-finite count of every per-particle field of every leg."""
    out: dict[str, int] = {}
    for name, leg in legs.items():
        if leg is None:
            continue
        for k, v in leg.items():
            if v.size == n and np.issubdtype(v.dtype, np.floating):
                bad = int((~np.isfinite(v)).sum())
                if bad:
                    out[f"{name}.{k}"] = bad
    return out


def evaluate(
    on: dict[str, Arr],
    off: dict[str, Arr],
    rep: dict[str, Arr] | None,
    synthetic: bool,
    x_on_scale: float,
    t_line_power: float,
    column_scale: float = 1.0,
    lead_on: Arr | None = None,
    lead_off: Arr | None = None,
    noise_floor: float | None = None,
) -> dict:
    """Return the per-bin discrepancy, bar terms and the mask counts.

    Parameters
    ----------
    on, off, rep
        Gas fields of the on leg, the off leg and the optional repeat of the on leg.
    synthetic
        Harness check: take the off-leg H2 from the solver.
    x_on_scale
        Planted scale of the on-leg H2 in the measured ratio only.
    t_line_power
        Planted power of the written T_line in the solve only (1 for none).
    column_scale
        Planted scale of the Jeans column in the solve only (1 for none).
    lead_on, lead_off
        Change of temperature per step [K] of every particle in each leg, or None.
    noise_floor
        Harness only: global noise in ln units used instead of the repeat run.
    """
    son, soff = state(on), state(off)
    y_on = son["n_h2"]
    x_on = y_on / son["n_h"]
    t_line = on["t_line"] ** t_line_power
    res: dict = {"n_all": y_on.size}
    bad = finite_counts({"on": on, "off": off, "repeat": rep}, y_on.size)
    res["nonfinite_fields"] = bad
    res["n_nonfinite"] = int(sum(bad.values()))
    res["finite"] = res["n_nonfinite"] == 0
    if not res["finite"]:
        return res
    y_max = 0.5 * soff["n_h"]
    lit = (t_line > 0.0) & (on["u_lw"] > 0.0) & (off["u_lw"] > 0.0) & (y_on > 0.0)
    res["n_unlit"] = int((~lit).sum())
    res["n_lit"] = int(lit.sum())
    k_unit = KDISS_DRAINE * C_LIGHT / FLUX_DRAINE
    k_on = k_unit * on["rho"] * on["u_lw"]
    k_off = k_unit * off["rho"] * off["u_lw"]

    def sub(leg: dict[str, Arr], mask: NDArray[np.bool_]) -> dict[str, Arr]:
        return {k: v[mask] for k, v in leg.items() if v.size == mask.size}

    def sub_state(st: dict[str, Arr], mask: NDArray[np.bool_]) -> dict[str, Arr]:
        return {k: v[mask] for k, v in st.items()}

    f_on = np.full(y_on.shape, np.nan)
    f_off = np.full(y_on.shape, np.nan)
    f_on[lit] = formation_rate(sub(on, lit), sub_state(son, lit), son["temp"][lit])
    f_off[lit] = formation_rate(sub(off, lit), sub_state(soff, lit), soff["temp"][lit])
    bad_f = lit & ~(
        (np.isfinite(f_on) & (f_on >= 0.0)) & (np.isfinite(f_off) & (f_off >= 0.0))
    )
    res["n_bad_formation"] = int(bad_f.sum())
    res["n_no_formation"] = int((lit & ~bad_f & ((f_on == 0.0) | (f_off == 0.0))).sum())
    lit = lit & ~bad_f & (f_on > 0.0) & (f_off > 0.0)

    def derive(
        mask: NDArray[np.bool_], t_on: Arr, t_off: Arr
    ) -> tuple[Arr, NDArray[np.int64]]:
        """Solve x_off S(x_off n_H,off L_off) = x_on T_line F_off/F_on k_on/k_off."""
        length = jeans_length(t_off, off["rho"][mask], soff["mu"][mask]) * column_scale
        f_ratio = formation_rate(sub(off, mask), sub_state(soff, mask), t_off) / (
            formation_rate(sub(on, mask), sub_state(son, mask), t_on)
        )
        k_ratio = k_on[mask] / k_off[mask]
        y_eff = x_on[mask] * f_ratio * k_ratio * soff["n_h"][mask]
        return solve_off(
            y_eff,
            t_line[mask],
            t_off,
            soff["n_tot"][mask],
            length,
            y_max[mask],
        )

    y_pred = np.full(y_on.shape, np.nan)
    count = np.zeros(y_on.shape, dtype=np.int64)
    y_pred[lit], count[lit] = derive(lit, son["temp"][lit], soff["temp"][lit])
    unique = lit & (count == 1)
    res["n_multi"] = int((lit & (count >= 2)).sum())
    res["n_noroot"] = int((lit & (count == 0)).sum())
    length_off = jeans_length(soff["temp"], off["rho"], soff["mu"]) * column_scale
    col_pred = np.where(unique, y_pred * length_off, 1.0)
    s_root = slope(col_pred, soff["temp"], soff["n_tot"])
    res["n_slope_ge1_kept"] = int((unique & (s_root >= 1.0)).sum())
    f_loc = wg_factor(col_pred, soff["temp"], soff["n_tot"])
    window = on["since_star"][0]
    with np.errstate(divide="ignore", invalid="ignore"):
        t_rel_on = 1.0 / (k_on * t_line)
        t_rel_off = 1.0 / (k_off * f_loc)
    eq = unique & (t_rel_on <= window) & (t_rel_off <= window)
    res["n_unique"] = int(unique.sum())
    res["n_not_equilibrated"] = int((unique & ~eq).sum())
    x_pred = y_pred / soff["n_h"]
    x_meas_off = x_pred if synthetic else soff["n_h2"] / soff["n_h"]
    x_rep = (
        state(rep)["n_h2"] / state(rep)["n_h"]
        if rep is not None
        else np.ones_like(y_on)
    )
    ok_h2 = (x_meas_off > 0.0) & (y_on > 0.0) & (x_rep > 0.0)
    res["n_zero_h2"] = int((eq & ~ok_h2).sum())
    keep = eq & ok_h2
    res["n_kept"] = int(keep.sum())
    res["dlnn"] = (
        float(np.median(np.abs(np.log(soff["n_tot"][keep] / son["n_tot"][keep]))))
        if keep.any()
        else float("nan")
    )
    res["dlnT"] = (
        float(np.median(np.abs(np.log(soff["temp"][keep] / son["temp"][keep]))))
        if keep.any()
        else float("nan")
    )
    meas = np.log(x_on_scale * x_on / x_meas_off)
    rep_ln = np.log(x_on / x_rep)
    pred = np.log(x_on / x_pred)
    resid = meas - pred
    # One-step lead of the library temperature: both legs at T minus the change of one step.
    res["lag_known"] = lead_on is not None and lead_off is not None
    lag_shift = np.full(y_on.shape, np.nan)
    res["n_lag_dropped"] = 0
    if res["lag_known"]:
        t_on_lib = np.maximum(son["temp"] - lead_on, 0.5 * son["temp"])
        t_off_lib = np.maximum(soff["temp"] - lead_off, 0.5 * soff["temp"])
        ok_t = keep & np.isfinite(t_on_lib) & np.isfinite(t_off_lib)
        yp = np.full(y_on.shape, np.nan)
        yp[ok_t], _ = derive(ok_t, t_on_lib[ok_t], t_off_lib[ok_t])
        lag_shift = np.log(x_on / (yp / soff["n_h"])) - pred
        res["n_lag_dropped"] = int((keep & ~np.isfinite(lag_shift)).sum())
    rng = np.random.default_rng(1)
    bins = {}
    for b, (lo, hi) in BINS.items():
        m = keep & (son["n_h"] >= lo) & (son["n_h"] < hi)
        if m.sum() == 0:
            bins[b] = {"n": 0}
            continue
        se = boot_se(resid[m], rng)
        noise = float("nan")
        if rep is not None:
            r = rep_ln[m]
            noise = float(abs(np.median(r)) + 2.0 * boot_se(r, rng))
        elif noise_floor is not None:
            noise = noise_floor
        lag = float("nan")
        if res["lag_known"] and np.isfinite(lag_shift[m]).any():
            lag = float(abs(np.median(lag_shift[m][np.isfinite(lag_shift[m])])))
        bins[b] = {
            "noise": noise,
            "n": int(m.sum()),
            "meas": float(np.median(meas[m])),
            "pred": float(np.median(pred[m])),
            "disc": float(np.median(resid[m])),
            "se": float(se),
            "lag": lag,
            "bar": float(2.0 * se + noise + lag),
        }
    res["bins"] = bins
    return res


def verdict(
    res: dict, noise: float | None, repeat: bool, accurate: bool = True
) -> tuple[int, str]:
    """Return the exit code and the reason.

    noise is a global harness-only floor in ln units; with a repeat run each
    bin carries its own measured noise, which is part of its bar. A bin is
    decided only with at least N_MIN particles and a bar of at most MAX_BAR.
    With accurate False (cheap-solver legs) no bin is decided: the cheap solver
    lags its equilibrium while the H2 falls, which the bar does not hold.
    """
    if not res["finite"]:
        fields = ", ".join(f"{k} {v}" for k, v in res["nonfinite_fields"].items())
        return (
            1,
            f"FAIL: {res['n_nonfinite']} non-finite values in the inputs ({fields})",
        )
    if res["n_bad_formation"]:
        return (
            1,
            f"FAIL: {res['n_bad_formation']} particles with a negative or non-finite formation rate",
        )
    good = {b: v for b, v in res["bins"].items() if "meas" in v}
    if res["n_kept"] == 0 or not good:
        return 3, "NOT EVALUATED: no particle or no bin left after the masks"
    if not accurate:
        return (
            2,
            "INDETERMINATE: cheap-solver legs lag their equilibrium while the H2 falls "
            "(s/(1-s) times the last sub-step change of ln x, at most dt/2); "
            "rerun both legs with GeajayCooling:solver 1",
        )
    if noise is None and not repeat:
        return 2, "INCONCLUSIVE: noise floor between two runs is UNKNOWN"
    if not res["lag_known"]:
        return (
            2,
            "INDETERMINATE: the one-step temperature lag is unknown (no previous snapshot or timesteps.txt)",
        )
    for b, v in good.items():
        if not np.isfinite([v["disc"], v["bar"]]).all():
            return 1, f"FAIL: non-finite result in {b}"
    decided = {b: v for b, v in good.items() if v["n"] >= N_MIN and v["bar"] <= MAX_BAR}
    if not decided:
        return (
            2,
            "INDETERMINATE: no bin with enough particles and a bar below the maximum",
        )
    worst = max(abs(v["disc"]) / v["bar"] for v in decided.values())
    skipped = len(good) - len(decided)
    if worst > 1.0:
        return 1, f"FAIL: discrepancy {worst:.2f} times the bar in a decided bin"
    return (
        0,
        f"PASS: largest discrepancy {worst:.2f} of the bar in {len(decided)} decided bins ({skipped} indeterminate)",
    )


def plot(rows: list[tuple[int, dict]], path: Path) -> None:
    """Plot measured and derived median ln ratio per density bin, all snapshots."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    names = list(BINS)
    seen = False
    for i, (snap, r) in enumerate(rows):
        for j, b in enumerate(names):
            v = r["bins"].get(b, {})
            if "meas" not in v:
                continue
            x = j + (i - len(rows) / 2) * 0.08
            if np.isfinite(v["noise"]):
                ax.errorbar(
                    x + 0.03,
                    v["pred"],
                    yerr=v["noise"],
                    fmt="none",
                    ecolor="0.6",
                    elinewidth=3,
                    label="noise floor (repeat run), half-width" if not seen else None,
                )
            ax.errorbar(
                x,
                v["meas"],
                yerr=v["bar"],
                fmt="o",
                color="C0",
                ms=4,
                label="measured median, bar half-width" if not seen else None,
            )
            ax.plot(
                x,
                v["pred"],
                "_",
                color="C3",
                ms=12,
                mew=2,
                label="derived (equilibrium solve)" if not seen else None,
            )
            seen = True
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([f"{lo:g} to {hi:g}" for lo, hi in BINS.values()])
    ax.set_xlabel("n_H of the on-leg gas [cm^-3]")
    ax.set_ylabel("median ln(x on / x off), x = n_H2 / n_H")
    ax.axhline(0.0, color="k", lw=0.5)
    ax.legend(fontsize=8)
    ax.set_title(
        "Snapshots "
        + ", ".join(str(s) for s, _ in rows)
        + " (left to right in each bin)",
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=110, pil_kwargs={"quality": 70})


def temperature_lead(snap_dir: Path, s: int, cur: dict[str, Arr]) -> Arr | None:
    """Return the change of T per step [K] of each particle of cur, or None.

    The change is that between snapshot s - 1 and s of the same leg divided by the
    number of steps of timesteps.txt (next to the snapshot directory) in between:
    an average over the interval.
    """
    prev_path = snap_dir / f"snapshot_{s - 1:04d}.hdf5"
    steps = snap_dir.parent / "timesteps.txt"
    if not prev_path.exists() or not steps.exists():
        return None
    prev = read_leg(prev_path, False)
    t_now, t_prev = cur["time_raw"][0], prev["time_raw"][0]
    times = np.loadtxt(steps, comments="#", usecols=(1,), ndmin=1)
    n_steps = int(((times > t_prev * (1 + 1e-9)) & (times <= t_now * (1 + 1e-9))).sum())
    if n_steps < 1:
        return None
    pos = np.searchsorted(prev["id"], cur["id"])
    pos = np.clip(pos, 0, prev["id"].size - 1)
    found = prev["id"][pos] == cur["id"]
    temp_prev = np.where(found, state(prev)["temp"][pos], np.nan)
    return (state(cur)["temp"] - temp_prev) / n_steps


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--on", type=Path, required=True, help="snapshot directory of the on leg"
    )
    ap.add_argument(
        "--off", type=Path, required=True, help="snapshot directory of the off leg"
    )
    ap.add_argument(
        "--repeat",
        type=Path,
        default=None,
        help="snapshot directory of a second run of the on leg, for the noise floor",
    )
    ap.add_argument("--first", type=int, default=3)
    ap.add_argument("--last", type=int, default=6)
    ap.add_argument(
        "--noise",
        type=float,
        default=None,
        help="harness only: global noise floor in ln units, instead of --repeat",
    )
    ap.add_argument(
        "--synthetic",
        action="store_true",
        help="harness check: take the off value from the solver",
    )
    ap.add_argument(
        "--x-on-scale",
        type=float,
        default=1.0,
        help="planted scale of the on-leg H2 in the measured ratio only",
    )
    ap.add_argument(
        "--t-line-power",
        type=float,
        default=1.0,
        help="planted power of the written T_line in the solve only (harness)",
    )
    ap.add_argument(
        "--column-scale",
        type=float,
        default=1.0,
        help="planted scale of the Jeans column in the solve only (harness)",
    )
    ap.add_argument(
        "--plant-nan-field",
        default=None,
        metavar="LEG:FIELD",
        help="harness: plant NaN in 5 particles of a field, e.g. on:H2I, off:z_smooth, off:u_pe",
    )
    ap.add_argument(
        "--plant-hm-sign",
        action="store_true",
        help="harness: flip the sign of the H- channel in the formation rate",
    )
    ap.add_argument("--plot", type=Path, default=None)
    a = ap.parse_args()
    global HM_SIGN
    if a.plant_hm_sign:
        HM_SIGN = -1.0
    rows = []
    codes = []
    for s in range(a.first, a.last + 1):
        on = read_leg(a.on / f"snapshot_{s:04d}.hdf5", True)
        off = read_leg(a.off / f"snapshot_{s:04d}.hdf5", False)
        rep = (
            read_leg(a.repeat / f"snapshot_{s:04d}.hdf5", False)
            if a.repeat is not None
            else None
        )
        legs = {"on": on, "off": off, "repeat": rep}
        names = {
            k: SOLVER_NAMES.get(int(v["solver"][0]), "unknown")
            for k, v in legs.items()
            if v is not None and np.isfinite(v["solver"][0])
        }
        print(
            "snapshot %d solver: " % s
            + ", ".join(
                f"{k} {names.get(k, 'not recorded')}"
                for k, v in legs.items()
                if v is not None
            )
        )
        if len(set(names.values())) != 1 or len(names) != sum(
            v is not None for v in legs.values()
        ):
            codes.append(2)
            print(
                f"snapshot {s}: exit 2 INDETERMINATE: legs differ or do not record the solver"
            )
            continue
        accurate = next(iter(names.values())) == "accurate (BDF)"
        common = np.intersect1d(on["id"], off["id"])
        if rep is not None:
            common = np.intersect1d(common, rep["id"])
        on, off = (pair(d, common) for d in (on, off))
        rep = pair(rep, common) if rep is not None else None
        lead_on = temperature_lead(a.on, s, on)
        lead_off = temperature_lead(a.off, s, off)
        if a.plant_nan_field is not None:
            leg_name, field = a.plant_nan_field.split(":")
            target = {"on": on, "off": off, "repeat": rep}[leg_name]
            target[field] = target[field].copy()
            target[field][:5] = np.nan
        r = evaluate(
            on,
            off,
            rep,
            a.synthetic,
            a.x_on_scale,
            a.t_line_power,
            a.column_scale,
            lead_on,
            lead_off,
            a.noise,
        )
        code, why = verdict(r, a.noise, rep is not None, accurate)
        codes.append(code)
        rows.append((s, r))
        print(f"snapshot {s}: exit {code} {why}")
        if not r["finite"]:
            print(
                f"  non-finite values: {r['n_nonfinite']} in fields {r['nonfinite_fields']}"
            )
        if r["finite"]:
            print(
                f"  of {r['n_all']}: no light {r['n_unlit']}, no formation {r['n_no_formation']}, "
                f"bad formation rate {r['n_bad_formation']}, no root {r['n_noroot']}, "
                f"more than one root up to n_H/2 {r['n_multi']}, not equilibrated {r['n_not_equilibrated']}, "
                f"zero H2 {r['n_zero_h2']}, kept {r['n_kept']} "
                f"(lag dropped {r['n_lag_dropped']}), "
                f"kept with s>=1 {r['n_slope_ge1_kept']}, median |dln n| {r['dlnn']:.3g}, median |dln T| {r['dlnT']:.3g}"
            )
            for b, v in r["bins"].items():
                if "meas" in v:
                    tag = (
                        "informational (cheap solver)"
                        if not accurate
                        else "decided"
                        if (v["n"] >= N_MIN and v["bar"] <= MAX_BAR)
                        else "INDETERMINATE"
                    )
                    print(
                        f"  {b} {tag}: n {v['n']}, measured {v['meas']:+.3f}, derived {v['pred']:+.3f}, "
                        f"median residual {v['disc']:+.2e}, discrepancy/bar {abs(v['disc']) / v['bar']:.2f}, "
                        f"bar {v['bar']:.3g} = 2 se {2 * v['se']:.3g} + noise {v['noise']:.3g} + lag {v['lag']:.3g}"
                    )
                else:
                    print(f"  {b}: no particle")
    if a.plot is not None:
        plot([(s, r) for s, r in rows if r["finite"] and r.get("bins")], a.plot)
    if 1 in codes or all(c == 3 for c in codes):
        return 1
    return 0 if all(c == 0 for c in codes) else 2


if __name__ == "__main__":
    sys.exit(main())
