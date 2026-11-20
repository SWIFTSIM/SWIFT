#!/usr/bin/env python3
"""Compare the H2 abundance of two legs, line absorption on and off.

The on leg writes the line transmission T_line of every gas particle. In
photodissociation equilibrium, with x = n_H2 / n_H, the formation rate per H
atom over the star rate is x_on T_line with the switch on and x_off S(n_H2,off L_J)
with the switch off. The two legs hold different gas: the formation rate per H
atom is F = h2dust(T) n_H (geajay rates.c h2dust_rate, element_blocks.h:191, rhoH)
and the star rate k_* is proportional to rho u_LW, so
x_off S(x_off n_H,off L_off) = x_on T_line (F_off / F_on) (k_on / k_off)
with each leg's own density, temperature, mean molecular weight and Jeans length
L_off. The script solves this for x_off of every gas particle and compares the
median of ln(x_on / x_off) per density bin with the measured one. The old check
compared n_H2 and solved with the on-leg density, temperature and column for
the off leg, which is exact only for two legs with equal gas. The solve is independent of the SWIFT
code: it uses the Wolcott-Green and Haiman (2019) fit written out below.

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

The noise floor is measured with --repeat, a second run of the on leg: per
density bin it is |median ln(n_H2,on / n_H2,repeat)| plus twice the bootstrap
standard error of that median, on the particles kept for the comparison. A bin
is decided (PASS or FAIL) only where its bar exceeds that noise; otherwise it
is INCONCLUSIVE. The bar is never widened.

The root search of each particle stops at n_H2 = n_H / 2. A particle is masked
only if it has no light, no root, more than one root in that range, or is not in
photodissociation equilibrium; the counts are printed by reason. Any non-finite
input is counted and fails the snapshot.

A bin is decided only with at least N_MIN particles, a bar of at most MAX_BAR
and a bar above its noise. Otherwise it is INDETERMINATE and decides nothing.

Exit status per snapshot: 0 every decided bin agrees and at least one bin is
decided, 1 failure (non-finite value, or a decided bin beyond its bar),
2 inconclusive or indeterminate (noise floor unknown, or no decided bin),
3 nothing left after the masks. The run exits 1 on any 1, 0 only if every
snapshot is 0, else 2.
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
# Dust temperature of the H2 formation rate. It is not in the snapshots. For T_dust
# below about 40 K the Omukai (2000) term in T_dust is 1; the bar carries the change
# between T_DUST_LOW and T_DUST_HIGH.
T_DUST = 10.0
T_DUST_LOW, T_DUST_HIGH = 3.0, 50.0
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


def h2dust_shape(temp: Arr, t_dust: float) -> Arr:
    """Return the H2-on-dust rate coefficient of Omukai (2000) up to a constant.

    geajay rates.c h2dust_rate for h2_dust_rate 1 (the library default).
    """
    term1 = 1.0 / (1.0 + np.exp(7.5e2 * (1.0 / 75.0 - 1.0 / t_dust)))
    term2 = 1.0 / (
        1.0 + 4.0e-2 * np.sqrt(temp + t_dust) + 2.0e-3 * temp + 8.0e-6 * temp**2
    )
    return np.sqrt(temp / 300.0) * term1 * term2


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
        if line_fields:
            out["t_line"] = g["LWLineTransmissions"][:][order].astype(float)
        t0 = float(f["Header"].attrs["Time"][0]) * ut
        births = f["PartType4/BirthTimes"][:] * ut if "PartType4" in f else np.array([])
        out["since_star"] = np.array([t0 - births.min() if births.size else np.nan])
        out["time"] = np.array([t0])
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


def evaluate(
    on: dict[str, Arr],
    off: dict[str, Arr],
    rep: dict[str, Arr] | None,
    synthetic: bool,
    x_on_scale: float,
    nan_plant: bool,
    t_line_power: float,
    column_scale: float = 1.0,
) -> dict:
    """Return the per-bin measured and derived medians and the mask counts.

    Parameters
    ----------
    on, off, rep
        Gas fields of the on leg, the off leg and the optional repeat of the on leg.
    synthetic
        Harness check: take the off-leg H2 from the solver.
    x_on_scale
        Planted scale of the on-leg H2 in the measured ratio only.
    nan_plant
        Plant NaN in the on-leg H2.
    t_line_power
        Planted power of the written T_line in the solve only (1 for none).
    column_scale
        Planted scale of the Jeans column in the solve only (1 for none).
    """
    son, soff = state(on), state(off)
    y_on = son["n_h2"]
    if nan_plant:
        y_on = y_on.copy()
        y_on[:5] = np.nan
    x_on = y_on / son["n_h"]
    t_line = on["t_line"] ** t_line_power
    res: dict = {"n_all": y_on.size}
    fields = [y_on, t_line, on["u_lw"], off["u_lw"], son["temp"], son["n_tot"]]
    fields += [soff["n_h2"], soff["temp"], soff["n_tot"]]
    if rep is not None:
        fields.append(state(rep)["n_h2"])
    res["n_nonfinite"] = int(sum((~np.isfinite(f)).sum() for f in fields))
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

    def derive(
        mask: NDArray[np.bool_], t_fac: float, t_dust: float
    ) -> tuple[Arr, NDArray[np.int64]]:
        """Solve x_off S(x_off n_H,off L_off) = x_on T_line F_off/F_on k_on/k_off."""
        t_on, t_off = son["temp"][mask] * t_fac, soff["temp"][mask] * t_fac
        length = (
            jeans_length(t_off, off["rho"][mask], soff["mu"][mask]) * column_scale
        )
        f_ratio = (soff["n_h"][mask] / son["n_h"][mask]) * (
            h2dust_shape(t_off, t_dust) / h2dust_shape(t_on, t_dust)
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
    y_pred[lit], count[lit] = derive(lit, 1.0, T_DUST)
    unique = lit & (count == 1)
    res["n_multi"] = int((lit & (count >= 2)).sum())
    res["n_noroot"] = int((lit & (count == 0)).sum())
    length_off = jeans_length(soff["temp"], off["rho"], soff["mu"]) * column_scale
    col_pred = np.where(unique, y_pred * length_off, 1.0)
    s_root = slope(col_pred, soff["temp"], soff["n_tot"])
    res["n_slope_ge1_kept"] = int((unique & (s_root >= 1.0)).sum())
    f_off = wg_factor(col_pred, soff["temp"], soff["n_tot"])
    window = on["since_star"][0]
    with np.errstate(divide="ignore", invalid="ignore"):
        t_rel_on = 1.0 / (k_on * t_line)
        t_rel_off = 1.0 / (k_off * f_off)
    eq = unique & (t_rel_on <= window) & (t_rel_off <= window)
    res["n_unique"] = int(unique.sum())
    res["n_not_equilibrated"] = int((unique & ~eq).sum())
    x_pred = y_pred / soff["n_h"]
    x_meas_off = x_pred if synthetic else soff["n_h2"] / soff["n_h"]
    x_rep = (
        state(rep)["n_h2"] / state(rep)["n_h"] if rep is not None else np.ones_like(y_on)
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
    # Sensitivity of the derived value to a 10 per cent error in T (gamma, mu) of
    # both legs, and to the unknown dust temperature of the formation rate.
    sens = {}
    res["n_sens_dropped"] = 0
    for key, t_fac, t_dust in (
        ("hot", 1.1, T_DUST),
        ("cold", 1 / 1.1, T_DUST),
        ("td_low", 1.0, T_DUST_LOW),
        ("td_high", 1.0, T_DUST_HIGH),
    ):
        yp = np.full(y_on.shape, np.nan)
        yp[keep], _ = derive(keep, t_fac, t_dust)
        sens[key] = np.log(x_on / (yp / soff["n_h"]))[keep]
        res["n_sens_dropped"] += int((~np.isfinite(sens[key])).sum())
    rng = np.random.default_rng(1)
    bins = {}
    for b, (lo, hi) in BINS.items():
        m = keep & (son["n_h"] >= lo) & (son["n_h"] < hi)
        if m.sum() == 0:
            bins[b] = {"n": 0}
            continue
        idx_in_keep = np.nonzero(m[keep])[0]


        def half_range(a: str, b: str) -> float:
            ok = np.isfinite(sens[a][idx_in_keep]) & np.isfinite(sens[b][idx_in_keep])
            if not ok.any():
                return float("nan")
            return 0.5 * abs(
                np.median(sens[a][idx_in_keep][ok]) - np.median(sens[b][idx_in_keep][ok])
            )

        sys_t = half_range("hot", "cold") + half_range("td_low", "td_high")
        se = np.hypot(boot_se(meas[m], rng), boot_se(pred[m], rng))
        noise = float("nan")
        if rep is not None:
            r = rep_ln[m]
            noise = float(abs(np.median(r)) + 2.0 * boot_se(r, rng))
        bins[b] = {
            "noise": noise,
            "n": int(m.sum()),
            "meas": float(np.median(meas[m])),
            "pred": float(np.median(pred[m])),
            "se": float(se),
            "sys_T": float(sys_t),
            "bar": float(2.0 * se + sys_t),
        }
    res["bins"] = bins
    return res


def verdict(
    res: dict, noise: float | None, repeat: bool, accurate: bool = True
) -> tuple[int, str]:
    """Return the exit code and the reason.

    noise is a global harness-only floor in ln units; with a repeat run each
    bin carries its own measured noise. A bin is decided only with at least
    N_MIN particles, a bar of at most MAX_BAR and a bar above its noise.
    With accurate False (cheap-solver legs) no bin is decided: the cheap solver
    lags its equilibrium while the H2 falls, which the bar does not hold.
    """
    if not res["finite"]:
        return 1, f"FAIL: {res['n_nonfinite']} non-finite values in the inputs"
    good = {b: v for b, v in res["bins"].items() if "meas" in v}
    if res["n_kept"] == 0 or not good:
        return 3, "NOT EVALUATED: no particle or no bin left after the masks"
    for b, v in good.items():
        if not np.isfinite([v["meas"], v["pred"], v["bar"]]).all():
            return 1, f"FAIL: non-finite result in {b}"
    if not accurate:
        return (
            2,
            "INDETERMINATE: cheap-solver legs lag their equilibrium while the H2 falls "
            "(s/(1-s) times the last sub-step change of ln x, at most dt/2); "
            "rerun both legs with GeajayCooling:solver 1",
        )
    if noise is None and not repeat:
        return 2, "INCONCLUSIVE: noise floor between two runs is UNKNOWN"
    floor = {b: (v["noise"] if repeat else noise) for b, v in good.items()}
    if not all(np.isfinite(list(floor.values()))):
        return 1, "FAIL: non-finite noise floor"
    decided = {
        b: v
        for b, v in good.items()
        if v["n"] >= N_MIN and v["bar"] <= MAX_BAR and v["bar"] > floor[b]
    }
    if not decided:
        return (
            2,
            "INDETERMINATE: no bin with enough particles, a bar below the maximum and above the noise",
        )
    worst = max(abs(v["meas"] - v["pred"]) / v["bar"] for v in decided.values())
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
        "--plant-nan", action="store_true", help="plant NaN in the on-leg H2"
    )
    ap.add_argument("--plot", type=Path, default=None)
    a = ap.parse_args()
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
        r = evaluate(
            on,
            off,
            rep,
            a.synthetic,
            a.x_on_scale,
            a.plant_nan,
            a.t_line_power,
            a.column_scale,
        )
        code, why = verdict(r, a.noise, rep is not None, accurate)
        codes.append(code)
        rows.append((s, r))
        print(f"snapshot {s}: exit {code} {why}")
        if not r["finite"]:
            print(
                f"  non-finite values: {r['n_nonfinite']} of {r['n_all']} particles' inputs"
            )
        if r["finite"]:
            print(
                f"  of {r['n_all']}: no light {r['n_unlit']}, no root {r['n_noroot']}, "
                f"more than one root up to n_H/2 {r['n_multi']}, not equilibrated {r['n_not_equilibrated']}, "
                f"zero H2 {r['n_zero_h2']}, kept {r['n_kept']} "
                f"(T-variation dropped {r['n_sens_dropped']}), "
                f"kept with s>=1 {r['n_slope_ge1_kept']}, median |dln n| {r['dlnn']:.3g}, median |dln T| {r['dlnT']:.3g}"
            )
            for b, v in r["bins"].items():
                if "meas" in v:
                    tag = (
                        "informational (cheap solver)"
                        if not accurate
                        else "decided"
                        if (
                            v["n"] >= N_MIN
                            and v["bar"] <= MAX_BAR
                            and v["bar"]
                            > (
                                v["noise"]
                                if rep is not None
                                else (a.noise or float("inf"))
                            )
                        )
                        else "INDETERMINATE"
                    )
                    print(
                        f"  {b} {tag}: n {v['n']}, measured {v['meas']:+.3f}, derived {v['pred']:+.3f}, "
                        f"discrepancy/bar {abs(v['meas'] - v['pred']) / v['bar']:.2f}, bar {v['bar']:.2g} (se {v['se']:.3f}, T-model {v['sys_T']:.3f}), noise {v['noise']:.2g}"
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
