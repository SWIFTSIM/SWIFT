#!/usr/bin/env python3
"""Compare the H2 abundance of two legs, line absorption on and off.

The on leg writes the line transmission T_line of every gas particle. In
photodissociation equilibrium the formation rate over the star rate is
x_on T_line with the switch on and x_off S(n_H2,off L_J) with the switch off.
The script solves the second equation for the off-leg H2 number density from
the on-leg state and compares the median of ln(n_H2,on / n_H2,off) per
density bin with the measured one. The solve is independent of the SWIFT
code: it uses the Wolcott-Green and Haiman (2019) fit written out below.

The noise floor is measured with --repeat, a second run of the on leg: per
density bin it is |median ln(n_H2,on / n_H2,repeat)| plus twice the bootstrap
standard error of that median, on the particles kept for the comparison. A bin
is decided (PASS or FAIL) only where its bar exceeds that noise; otherwise it
is INCONCLUSIVE. The bar is never widened.

Exit status per snapshot: 0 every decided bin agrees and at least one bin is
decided, 1 failure (non-finite value, or a decided bin beyond its bar),
2 inconclusive (noise floor unknown, or no bin with a bar above the noise),
3 nothing left after the masks. The run exits 1 on any 1, else 0 if any
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


def jeans_length(temp: Arr, rho: Arr, mu: Arr) -> Arr:
    """Return the Jeans length in cm for gas of density rho in g cm^-3."""
    return np.sqrt(GAMMA * np.pi * KB * temp / (GRAV * mu * MH * rho))


def solve_off(
    y_on: Arr, t_line: Arr, temp: Arr, n_tot: Arr, length: Arr
) -> tuple[Arr, NDArray[np.int64]]:
    """Solve y f(y L) = y_on T_line for y, the H2 number density in cm^-3.

    Returns the root, NaN unless the root is unique, and the root count.
    """
    target = np.log(y_on * t_line)
    y_root = np.full(y_on.shape, np.nan)
    count = np.zeros(y_on.shape, dtype=np.int64)
    ly = np.log(Y_GRID)
    for i in range(y_on.size):
        s = wg_factor(Y_GRID * length[i], temp[i], n_tot[i])
        d = ly + np.log(s) - target[i]
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
        if line_fields:
            out["t_line"] = g["LWLineTransmissions"][:][order].astype(float)
            out["u_lw"] = (
                g["LWSpecificEnergies"][:][order].astype(float) * (ul / ut) ** 2
            )
        t0 = float(f["Header"].attrs["Time"][0]) * ut
        births = f["PartType4/BirthTimes"][:] * ut if "PartType4" in f else np.array([])
        out["since_star"] = np.array([t0 - births.min() if births.size else np.nan])
        out["time"] = np.array([t0])
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
    snap_dt: float,
) -> dict:
    """Return the per-bin measured and derived medians and the mask counts."""
    son, soff = state(on), state(off)
    y_on = son["n_h2"]
    if nan_plant:
        y_on = y_on.copy()
        y_on[:5] = np.nan
    t_line = on["t_line"]
    res: dict = {"n_all": y_on.size}
    need = np.concatenate(
        [
            y_on,
            t_line,
            on["u_lw"],
            son["temp"],
            son["n_tot"],
            soff["n_h2"],
            soff["temp"],
        ]
    )
    res["finite"] = bool(np.isfinite(need).all())
    if not res["finite"]:
        return res
    length = jeans_length(son["temp"], on["rho"], son["mu"])
    lit = (t_line > 0.0) & (on["u_lw"] > 0.0) & (y_on > 0.0)
    res["n_lit"] = int(lit.sum())
    y_pred = np.full(y_on.shape, np.nan)
    count = np.zeros(y_on.shape, dtype=np.int64)
    y_pred[lit], count[lit] = solve_off(
        y_on[lit], t_line[lit], son["temp"][lit], son["n_tot"][lit], length[lit]
    )
    unique = lit & (count == 1)
    res["n_multi"] = int((lit & (count >= 2)).sum())
    res["n_noroot"] = int((lit & (count == 0)).sum())
    col_pred = np.where(unique, y_pred * length, 1.0)
    s_root = slope(col_pred, son["temp"], son["n_tot"])
    res["n_slope_ge1_kept"] = int((unique & (s_root >= 1.0)).sum())
    k_star = (KDISS_DRAINE * C_LIGHT / FLUX_DRAINE) * on["rho"] * on["u_lw"]
    f_off = wg_factor(col_pred, son["temp"], son["n_tot"])
    window = on["since_star"][0]
    with np.errstate(divide="ignore", invalid="ignore"):
        t_on = 1.0 / (k_star * t_line)
        t_off = 1.0 / (k_star * f_off)
    eq = unique & (t_on <= window) & (t_off <= window)
    res["n_unique"] = int(unique.sum())
    res["n_not_equilibrated"] = int((unique & ~eq).sum())
    y_meas_off = y_pred if synthetic else soff["n_h2"]
    keep = eq & (y_meas_off > 0.0) & (y_on > 0.0)
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
    meas = np.log(x_on_scale * y_on / y_meas_off)
    if rep is not None:
        y_rep = state(rep)["n_h2"]
        res["finite"] = res["finite"] and bool(np.isfinite(y_rep).all())
        if not res["finite"]:
            return res
        with np.errstate(divide="ignore", invalid="ignore"):
            rep_ln = np.log(y_on / y_rep)
    pred = np.log(y_on / y_pred)
    # Sensitivity of the derived value to a 10 per cent error in T (gamma, mu).
    sens = {}
    for fac in (1.1, 1 / 1.1):
        tt = son["temp"] * fac
        ll = jeans_length(tt, on["rho"], son["mu"])
        yp, ct = solve_off(
            y_on[keep], t_line[keep], tt[keep], son["n_tot"][keep], ll[keep]
        )
        sens[fac] = np.log(y_on[keep] / yp)
    rng = np.random.default_rng(1)
    bins = {}
    for b, (lo, hi) in BINS.items():
        m = keep & (son["n_h"] >= lo) & (son["n_h"] < hi)
        if m.sum() < 10:
            bins[b] = {"n": int(m.sum())}
            continue
        idx_in_keep = np.nonzero(m[keep])[0]
        sys_t = 0.5 * abs(
            np.nanmedian(sens[1.1][idx_in_keep])
            - np.nanmedian(sens[1 / 1.1][idx_in_keep])
        )
        se = np.hypot(boot_se(meas[m], rng), boot_se(pred[m], rng))
        noise = float("nan")
        if rep is not None:
            r = rep_ln[m]
            r = r[np.isfinite(r)]
            noise = (
                float(abs(np.median(r)) + 2.0 * boot_se(r, rng))
                if r.size > 1
                else float("nan")
            )
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


def verdict(res: dict, noise: float | None, repeat: bool) -> tuple[int, str]:
    """Return the exit code and the reason.

    noise is a global harness-only floor in ln units; with a repeat run each
    bin carries its own measured noise.
    """
    if not res["finite"]:
        return 1, "FAIL: a non-finite value in the inputs"
    good = {b: v for b, v in res["bins"].items() if "meas" in v}
    if res["n_kept"] == 0 or not good:
        return 3, "NOT EVALUATED: no particle or no bin left after the masks"
    for b, v in good.items():
        if not np.isfinite([v["meas"], v["pred"], v["bar"]]).all():
            return 1, f"FAIL: non-finite result in {b}"
    if noise is None and not repeat:
        return 2, "INCONCLUSIVE: noise floor between two runs is UNKNOWN"
    floor = {b: (v["noise"] if repeat else noise) for b, v in good.items()}
    if not all(np.isfinite(list(floor.values()))):
        return 1, "FAIL: non-finite noise floor"
    decided = {b: v for b, v in good.items() if v["bar"] > floor[b]}
    if not decided:
        return 2, "INCONCLUSIVE: no bin has a bar above the noise floor"
    worst = max(abs(v["meas"] - v["pred"]) / v["bar"] for v in decided.values())
    skipped = len(good) - len(decided)
    if worst > 1.0:
        return 1, f"FAIL: discrepancy {worst:.2f} times the bar in a decided bin"
    return (
        0,
        f"PASS: largest discrepancy {worst:.2f} of the bar in {len(decided)} decided bins ({skipped} inconclusive)",
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
    ax.set_ylabel("median ln(n_H2 on / n_H2 off)")
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
    ap.add_argument("--first", type=int, default=43)
    ap.add_argument("--last", type=int, default=50)
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
        common = np.intersect1d(on["id"], off["id"])
        if rep is not None:
            common = np.intersect1d(common, rep["id"])
        on, off = (pair(d, common) for d in (on, off))
        rep = pair(rep, common) if rep is not None else None
        r = evaluate(on, off, rep, a.synthetic, a.x_on_scale, a.plant_nan, 0.0)
        code, why = verdict(r, a.noise, rep is not None)
        codes.append(code)
        rows.append((s, r))
        print(f"snapshot {s}: exit {code} {why}")
        if r["finite"]:
            print(
                f"  lit {r['n_lit']}, multi-root masked {r['n_multi']}, no-root masked {r['n_noroot']}, "
                f"not equilibrated {r['n_not_equilibrated']}, kept {r['n_kept']} of {r['n_all']}, "
                f"kept with s>=1 {r['n_slope_ge1_kept']}, median |dln n| {r['dlnn']:.3g}, median |dln T| {r['dlnT']:.3g}"
            )
            for b, v in r["bins"].items():
                if "meas" in v:
                    print(
                        f"  {b}: n {v['n']}, measured {v['meas']:+.3f}, derived {v['pred']:+.3f}, "
                        f"discrepancy/bar {abs(v['meas'] - v['pred']) / v['bar']:.2f}, bar {v['bar']:.2g} (se {v['se']:.3f}, T-model {v['sys_T']:.3f}), noise {v['noise']:.2g}"
                    )
                else:
                    print(f"  {b}: n {v['n']} (fewer than 10, not evaluated)")
    if a.plot is not None:
        plot([(s, r) for s, r in rows if r["finite"] and r.get("bins")], a.plot)
    if 1 in codes or all(c == 3 for c in codes):
        return 1
    return 0 if 0 in codes else 2


if __name__ == "__main__":
    sys.exit(main())
