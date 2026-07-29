#!/usr/bin/env python3
"""Multi-metric comparison of GEAR FVPM hyperbolic-diffusion runs against the
exact Cattaneo/telegrapher solution, resolving the actual discrete seed.

The exact reference is built per seed particle with closed-form kernels:
  - smooth interior: shell-redistribution kernel k(r;u,d) = r/(2ud) on
    |u-d| < r < u+d, integrated over the interior Green's function;
  - ballistic front: monopole shell of mass fraction e^{-at}(1+at+(at)^2/2)
    PLUS the delta' dipole (-ct e^{-at} d'(r-ct)), whose seeded form is a
    band -r/(2(ct)^2 d) with signed edge point-masses at r = ct +/- d.
    The dipole is genuine 3D wave physics (Huygens dipole layer); it makes
    the exact binned reference locally signed near the inner front edge.
The dipole's second moment 2(ct)^2 e^{-at} is required for <r^2> to satisfy
the exact moment ODE  <r^2> = 6 kappa (t - tau(1 - e^{-t/tau})), which is
also what `exact_rms` uses below (closed form, dimension-aware, plus the
actual seed's own second moment).

Line-density (x-profile) references use the identity that the y,z-marginal
of the 3D Green's function IS the 1D Green's function (verified to machine
precision), so any 3D run's x-profile reference is a sum of 1D solutions
over seed particles.

Stages:  --compute  scan runs, write metrics CSV + per-run profile cache
         --plot     read cache, write figures
"""
import os
import sys
import glob
import argparse
import pickle
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import swiftsimio as sw
from multiprocessing import Pool
from profile_vs_exact import u_3d, u_1d, metal_mass

ROOT = "hyperbolic.limiter_dev_20260725"
OUT = os.path.join(ROOT, "analysis_v2")
KERNEL_GAMMA = 1.936492  # Wendland C2, 3D

# ---------------------------------------------------------------- exact refs


def front_terms(t, tau, kappa):
    c = np.sqrt(kappa / tau)
    a = 0.5 / tau
    at = a * t
    e = np.exp(-at)
    return c, a, at, e, e * (1.0 + at + 0.5 * at**2)


def exact_rms(t, tau, kappa, ndim, seed_m2):
    """Exact RMS radius from the moment ODE (holds for the PDE regardless of
    the front decomposition), plus the actual seed's own second moment."""
    m2 = 2 * ndim * kappa * (t - tau * (1.0 - np.exp(-t / tau)))
    return np.sqrt(m2 + seed_m2)


def exact_radial_binned(edges, t, seed_d, seed_m, tau, kappa, include_dipole=True):
    """Exact Fe mass per radial bin for a discrete seed (3D). Returns
    (mass_per_bin, fine_r, fine_smooth_dMdr) — the fine curve excludes the
    front point-masses, which live only in the binned output."""
    c, a, at, e_at, ff = front_terms(t, tau, kappa)
    ct = c * t
    M = seed_m.sum()
    rmax = edges[-1]
    rg = np.linspace(0.0, rmax, 6001)
    mu = np.zeros_like(rg)  # smooth dM/dr
    points = []  # (r, signed mass) front point masses

    # interior Green's function radial mass density and its /u cumulative
    nu = 6001
    ug = np.linspace(0.0, ct, nu)
    muG = 4 * np.pi * ug**2 * u_3d(ug, t, tau, kappa)
    with np.errstate(divide="ignore", invalid="ignore"):
        psi = np.where(ug > 0, muG / ug, 0.0)
    Psi = np.concatenate([[0.0], np.cumsum(0.5 * (psi[1:] + psi[:-1]) * np.diff(ug))])
    d_floor = 0.75 * (rg[1] - rg[0])
    for d, m in zip(seed_d, seed_m):
        w = m / M
        if d < d_floor:
            # effectively a point source at the centre
            mu += w * np.interp(rg, ug, muG, left=0.0, right=0.0)
            points.append((ct, w * ff))
            # point-seed dipole: integral of delta' over any bin whose
            # interior contains ct is exactly zero -- no binned contribution
        else:
            # interior: (r/2d) * [Psi(min(r+d,ct)) - Psi(|r-d|)]
            hi = np.interp(np.minimum(rg + d, ct), ug, Psi)
            lo = np.interp(np.abs(rg - d), ug, Psi)
            mu += w * rg / (2 * d) * np.clip(hi - lo, 0.0, None)
            # shell-kernel support is |ct-d| < r < ct+d (also right for ct<d)
            band = (rg > abs(ct - d)) & (rg < ct + d)
            # front monopole shell redistributed over the band
            mu[band] += w * ff * rg[band] / (2 * ct * d)
            if include_dipole:
                # dipole: band part + signed edge point-masses; the lower
                # edge weight -(ct-d)/(2d) has the correct sign for ct<d too
                mu[band] += -w * e_at * rg[band] / (2 * ct * d)
                points.append((ct + d, w * e_at * (ct + d) / (2 * d)))
                points.append((abs(ct - d), -w * e_at * (ct - d) / (2 * d)))

    mass = np.zeros(len(edges) - 1)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (mu[1:] + mu[:-1]) * np.diff(rg))])
    mass += np.diff(np.interp(edges, rg, cum))
    for r0, pm in points:
        if r0 < edges[-1]:
            mass[min(np.searchsorted(edges, r0, side="right") - 1, len(mass) - 1)] += pm
    return M * mass, rg, M * mu


# ---------------------------------------------------------------- run scan


def load_seed(snapdir):
    f = os.path.join(snapdir, "snapshot_0000.hdf5")
    d = sw.load(f)
    m = metal_mass(d)
    box = d.metadata.boxsize.value
    ndim = 3 if (box[1] > 1e-3 and box[2] > 1e-3) else 1
    pos = d.gas.coordinates.value - box / 2.0
    sel = m > 0
    r = np.linalg.norm(pos[sel], axis=1) if ndim == 3 else np.abs(pos[sel, 0])
    spacing = box[0] / round(len(m) ** (1.0 / ndim))
    return dict(
        ndim=ndim,
        box=box,
        spacing=spacing,
        d=r,
        x=pos[sel, 0],
        m=m[sel],
        m2=(m[sel] * r**2).sum() / m[sel].sum(),
        rmax_seed=float(r.max()) if r.size else 0.0,
    )


def scan_run(rundir):
    try:
        return scan_run_inner(rundir)
    except Exception as exc:  # a corrupt run must not kill the whole sweep
        print(f"SCAN FAILED {rundir}: {type(exc).__name__}: {exc}")
        return None


def scan_run_inner(rundir):
    snapdir = os.path.join(rundir, "snap")
    files = sorted(glob.glob(os.path.join(snapdir, "snapshot_*.hdf5")))
    if len(files) < 3:
        return None
    try:
        p = sw.load(files[0]).metadata.parameters
        kappa = float(p["GEARChemistry:diffusion_coefficient"])
        tau = float(p["GEARChemistry:tau"])
        psi = float(p.get("GEARChemistry:hll_riemann_solver_psi", np.nan))
    except Exception:
        return None
    if tau <= 0 or kappa <= 0:
        return None
    seed = load_seed(snapdir)
    c = np.sqrt(kappa / tau)
    half = seed["box"][0] / 2.0
    rows, profiles = [], {}
    # profile times: three snapshots inside the clean (unsaturated) window
    t_clean = (half - seed["rmax_seed"]) / c
    want_t = [f * min(t_clean, 5.0) for f in (0.3, 0.6, 0.95)]

    for f in files[1:]:
        d = sw.load(f)
        t = float(d.metadata.time.value)
        if t <= 0:
            continue
        m = metal_mass(d)
        box = d.metadata.boxsize.value
        pos = d.gas.coordinates.value - box / 2.0
        r = np.linalg.norm(pos, axis=1) if seed["ndim"] == 3 else np.abs(pos[:, 0])
        M = m.sum()
        if M <= 0:
            continue
        ct = c * t
        front_hi = ct + seed["rmax_seed"] + 0.5 * seed["spacing"]
        saturated = front_hi > half
        h_med = float(np.median(d.gas.smoothing_lengths.value))
        reach = KERNEL_GAMMA * h_med

        rms_ex = exact_rms(t, tau, kappa, seed["ndim"], seed["m2"])
        row = dict(
            t=t,
            saturated=int(saturated),
            excess=m[r > front_hi].sum() / M if not saturated else np.nan,
            rms_ratio=np.sqrt((m * r**2).sum() / M) / rms_ex,
            # largest rms a POSITIVE causal solution can reach: the
            # exact signed dipole inflates <r^2> beyond the support
            rms_ceiling=front_hi / rms_ex,
            n_neg=int((m < 0).sum()),
            worst_neg=float(min(m.min(), 0)),
            front_capture=np.nan,
            L1=np.nan,
            W1=np.nan,
        )

        if seed["ndim"] == 3 and not saturated:
            # fine-grid exact cumulative: W1 + smooth window metrics
            # dipole-free (positive) reference: same mass, same causal
            # support; the fair transport target for a positive scheme
            fine = np.linspace(0.0, half, 2001)
            ex_fine, _, _ = exact_radial_binned(
                fine, t, seed["d"], seed["m"], tau, kappa, include_dipole=False
            )
            cum_ex = np.concatenate([[0.0], np.cumsum(ex_fine)])
            sim_fine, _ = np.histogram(r, bins=fine, weights=m)
            cum_sim = np.concatenate([[0.0], np.cumsum(sim_fine)])
            row["W1"] = np.trapezoid(np.abs(cum_sim - cum_ex), fine) / M
            wlo = max(ct - seed["rmax_seed"] - reach, 0.0)
            whi = front_hi + reach
            ex_w = np.interp(whi, fine, cum_ex) - np.interp(wlo, fine, cum_ex)
            sim_w = m[(r > wlo) & (r < whi)].sum()
            row["front_capture"] = sim_w / ex_w if ex_w > 0.05 * M else np.nan

            # end exactly at half so near-saturation front points stay in
            edges = np.append(np.arange(0.0, half, seed["spacing"]), half)
            sim_mass, _ = np.histogram(r, bins=edges, weights=m)
            ex_mass, rg, mu = exact_radial_binned(
                edges, t, seed["d"], seed["m"], tau, kappa
            )
            row["L1"] = 0.5 * np.abs(sim_mass - ex_mass).sum() / M
            tgt = min(want_t, key=lambda w: abs(w - t))
            if abs(t - tgt) < 0.051 and tgt not in profiles:
                profiles[tgt] = dict(
                    t=t,
                    edges=edges,
                    sim=sim_mass,
                    exact=ex_mass,
                    fine_r=rg[::10],
                    fine_mu=mu[::10],
                    ct=ct,
                    rseed=seed["rmax_seed"],
                    M=M,
                )
        rows.append(row)
    if not rows:
        return None
    rel = os.path.relpath(rundir, ROOT)
    return dict(
        run=rel,
        group=rel.split(os.sep)[0],
        name=os.path.basename(rundir),
        tau=tau,
        kappa=kappa,
        psi=psi,
        ndim=seed["ndim"],
        nseed=len(seed["m"]),
        spacing=seed["spacing"],
        rows=rows,
        profiles=profiles,
    )


def compute(nproc):
    rundirs = sorted(
        os.path.dirname(s) for s in glob.glob(os.path.join(ROOT, "*", "*", "snap"))
    )
    os.makedirs(OUT, exist_ok=True)
    with Pool(nproc) as pool:
        results = [r for r in pool.map(scan_run, rundirs) if r]
    with open(os.path.join(OUT, "cache.pkl"), "wb") as fh:
        pickle.dump(results, fh)
    with open(os.path.join(OUT, "metrics.csv"), "w") as fh:
        fh.write(
            "group,run,ndim,tau,psi,nseed,t,saturated,excess,"
            "rms_ratio,rms_ceiling,L1,W1,front_capture,n_neg,"
            "worst_neg\n"
        )
        for r in results:
            for row in r["rows"]:
                fh.write(
                    f"{r['group']},{r['name']},{r['ndim']},{r['tau']:g},"
                    f"{r['psi']:g},{r['nseed']},{row['t']:.4f},"
                    f"{row['saturated']},{row['excess']:.6f},"
                    f"{row['rms_ratio']:.6f},{row['rms_ceiling']:.6f},"
                    f"{row['L1']:.6f},{row['W1']:.6f},"
                    f"{row['front_capture']:.6f},{row['n_neg']},"
                    f"{row['worst_neg']:.3e}\n"
                )
    print(f"cached {len(results)} runs -> {OUT}")
    return results


# ---------------------------------------------------------------- plotting

# family -> (label, color); eras separate panels in the profile galleries
FAMILY = {
    "03_limiter_scope_compare/l5_tau_*_limiter_all": ("flux-limiter all", "#d62728"),
    "03_limiter_scope_compare/l5_tau_*_limiter_density": (
        "flux-limiter density",
        "#ff7f0e",
    ),
    "04_minmod/l5_tau_*_minmod": ("minmod psi=0.1", "#1f77b4"),
    "05_diagnostics/l5_tau_*_nodiss": ("F_diss=0", "#7f7f7f"),
    "05_diagnostics/l5_tau_*_noblend": ("no alpha-blend", "#bcbd22"),
    "09_muscl_fix/l5_tau_*_musclfix*": ("MUSCL-fix+minmod", "#2ca02c"),
    "13_psi_retune_postfix/l5_tau_*_psi*": ("psi sweep", "#9467bd"),
    "13_psi_retune_postfix/l5_tau_*_fdiss0_postfix": ("F_diss=0 postfix", "#8c564b"),
    "14_donorcell_diagnostic/*": ("donor-cell", "#e377c2"),
    "16_reduced_neighbours/*": ("eta=1.0", "#17becf"),
    "11_epsilon_fixed/l5*": ("l5 eps=0.04", "#2ca02c"),
    "11_epsilon_fixed/l6*": ("l6 eps=0.04", "#00441b"),
    "12_final_defaults/*": ("final defaults", "#000080"),
    "17_1d_check/*": ("1D check", "#654321"),
    "20_causal_diss/base_tau_*": ("baseline (controlled)", "#808080"),
    "20_causal_diss/causal_tau_*": ("causal F_diss bound", "#e6194b"),
    "21_more_tests/cfl0p9_*": ("CFL=0.9", "#3cb44b"),
    "21_more_tests/cfl0p1_*": ("CFL=0.1", "#4363d8"),
    "21_more_tests/psi20_base_*": ("psi=20 base", "#f58231"),
    "21_more_tests/psi20_causal_*": ("psi=20 + causal", "#911eb4"),
    "21_more_tests/l6_base_*": ("l6 baseline", "#469990"),
}


def family_of(run):
    import fnmatch

    for pat, (lab, colr) in FAMILY.items():
        if fnmatch.fnmatch(run["run"], pat):
            if lab == "psi sweep":
                lab = f"psi={run['psi']:g}"
            return lab, colr
    return run["name"], "#333333"


ERAS = [
    ("pre-fix era", ("03_", "04_", "05_")),
    ("MUSCL-fix era", ("09_", "13_", "14_", "16_")),
    ("eps-fixed defaults", ("11_", "12_")),
    ("causal-F_diss experiment", ("20_",)),
    ("CFL / psi20 / l6 tests", ("21_",)),
]


def plot_profiles(results):
    taus = sorted({r["tau"] for r in results if r["ndim"] == 3})
    for tau in taus:
        sel = [
            r for r in results if r["tau"] == tau and r["ndim"] == 3 and r["profiles"]
        ]
        if not sel:
            continue
        tkeys = sorted({round(k, 3) for r in sel for k in r["profiles"]})
        if not tkeys:
            continue
        tkeys = tkeys[-3:]
        eras = [
            (era, pref)
            for era, pref in ERAS
            if any(r["group"].startswith(pref) for r in sel)
        ]
        fig, axes = plt.subplots(
            len(tkeys),
            len(eras),
            figsize=(4.6 * len(eras), 3.4 * len(tkeys)),
            squeeze=False,
            sharex="row",
        )
        for i, tk in enumerate(tkeys):
            for j, (era, prefixes) in enumerate(eras):
                ax = axes[i][j]
                # one exact curve per distinct seed in the panel (l5 vs l6
                # seeds differ; a single reference would misattribute that)
                exact_seeds = {}
                seen_labels = set()
                for r in sel:
                    if not r["group"].startswith(prefixes):
                        continue
                    pk = min(r["profiles"], key=lambda k: abs(k - tk))
                    if abs(pk - tk) > 0.15:
                        continue
                    p = r["profiles"][pk]
                    wid = np.diff(p["edges"])
                    lab, colr = family_of(r)
                    if r["nseed"] not in exact_seeds:
                        ls = ["-", "--", ":"][len(exact_seeds) % 3]
                        exact_seeds[r["nseed"]] = True
                        ax.stairs(
                            p["exact"] / p["M"] / wid,
                            p["edges"],
                            color="k",
                            lw=2.2,
                            ls=ls,
                            label=f"exact ({r['nseed']}-part seed)",
                        )
                        ax.axvspan(
                            p["ct"] - p["rseed"],
                            p["ct"] + p["rseed"],
                            color="gold",
                            alpha=0.25,
                            lw=0,
                        )
                        ax.set_title(f"{era}   t={p['t']:.2f}", fontsize=10)
                    ax.stairs(
                        p["sim"] / p["M"] / wid,
                        p["edges"],
                        color=colr,
                        lw=1.4,
                        label=None if lab in seen_labels else lab,
                    )
                    seen_labels.add(lab)
                if exact_seeds:
                    ax.legend(fontsize=6.5)
                ax.grid(alpha=0.25)
                if j == 0:
                    ax.set_ylabel("dM/dr / M$_{tot}$ [kpc$^{-1}$]")
                if i == len(tkeys) - 1:
                    ax.set_xlabel("r [kpc]")
        fig.suptitle(
            f"Radial Fe mass profiles vs exact, tau={tau:g}\n"
            "gold band: exact ballistic front; negative dip inside "
            "the front = genuine Huygens dipole, not an error",
            fontsize=11,
        )
        plt.tight_layout()
        out = os.path.join(OUT, f"profiles_tau{tau:g}.png")
        plt.savefig(out, dpi=130)
        plt.close(fig)
        print("wrote", out)


def plot_metrics(results):
    taus = sorted({r["tau"] for r in results if r["ndim"] == 3})
    specs = [
        ("excess", "mass fraction beyond causal front", None),
        ("rms_ratio", "RMS radius / exact", 1.0),
        ("W1", "W1 distance to exact (dipole-free) [kpc]", None),
        ("front_capture", "front-window mass / exact", 1.0),
    ]
    for tau in taus:
        sel = [r for r in results if r["tau"] == tau and r["ndim"] == 3]
        if not sel:
            continue
        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        seen = set()
        for r in sel:
            lab, colr = family_of(r)
            show = lab not in seen
            seen.add(lab)
            t = np.array([w["t"] for w in r["rows"]])
            for ax, (key, ylab, ref) in zip(axes.flat, specs):
                v = np.array([w[key] for w in r["rows"]], dtype=float)
                ok = np.isfinite(v)
                ax.plot(t[ok], v[ok], color=colr, lw=1.4, label=lab if show else None)
        # positivity ceiling on rms: a positive causal solution cannot match
        # the exact signed <r^2>; plot the best it could do (front_hi/rms_ex)
        r0 = max(sel, key=lambda r: r["nseed"])
        tt = np.array([w["t"] for w in r0["rows"]])
        cc = np.array([w["rms_ceiling"] for w in r0["rows"]])
        axes.flat[1].plot(tt, cc, color="k", ls=":", lw=1.2, label="positivity ceiling")
        axes.flat[1].legend(fontsize=7)
        for ax, (key, ylab, ref) in zip(axes.flat, specs):
            if ref is not None:
                ax.axhline(ref, color="k", ls="--", lw=0.9)
            ax.set_xlabel("t")
            ax.set_ylabel(ylab)
            ax.grid(alpha=0.25)
        axes.flat[0].legend(fontsize=7, ncol=2)
        fig.suptitle(
            f"Metrics vs exact solution, tau={tau:g} "
            "(saturated times excluded from excess/W1)"
        )
        plt.tight_layout()
        out = os.path.join(OUT, f"metrics_tau{tau:g}.png")
        plt.savefig(out, dpi=130)
        plt.close(fig)
        print("wrote", out)


def plot_psi_sweep(results):
    sel = [
        r
        for r in results
        if r["group"] == "13_psi_retune_postfix" and np.isfinite(r["psi"])
    ]
    byt = {}
    for r in sel:
        last = [w for w in r["rows"] if not w["saturated"]]
        if not last:
            continue
        w = last[-1]
        byt.setdefault(r["tau"], []).append(
            (r["psi"], w["excess"], w["W1"], w["rms_ratio"])
        )
    if not byt:
        return
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    for tau, rows in sorted(byt.items()):
        rows.sort()
        p, e, l1, rr = map(np.array, zip(*rows))
        for ax, v, ylab in zip(
            axes, (e, l1, rr), ("excess beyond front", "W1 to exact [kpc]", "RMS/exact")
        ):
            ax.semilogx(p, v, "o-", label=f"tau={tau:g}")
            ax.set_xlabel("psi (minmod bound)")
            ax.set_ylabel(ylab)
            ax.grid(alpha=0.25)
    axes[0].legend()
    fig.suptitle("psi sweep (13_psi_retune_postfix), last clean snapshot")
    plt.tight_layout()
    out = os.path.join(OUT, "psi_sweep.png")
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print("wrote", out)


def plot_ranking(results):
    rows = []
    for r in results:
        if r["ndim"] != 3:
            continue
        clean = [w for w in r["rows"] if not w["saturated"] and np.isfinite(w["W1"])]
        if not clean:
            continue
        w = clean[-1]
        rows.append(
            (
                r["tau"],
                f"{r['group'][:2]}/{r['name']}",
                w["W1"],
                w["excess"],
                family_of(r)[1],
                w["t"],
            )
        )
    if not rows:
        return
    rows.sort(key=lambda x: (x[0], x[2]))
    fig, ax = plt.subplots(figsize=(9, 0.32 * len(rows) + 1.5))
    y = np.arange(len(rows))
    ax.barh(y, [x[2] for x in rows], color=[x[4] for x in rows])
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"tau={x[0]:g}  {x[1]}  (t={x[5]:.1f})" for x in rows], fontsize=7
    )
    ax.invert_yaxis()
    ax.set_xlabel(
        "W1 distance to exact (dipole-free) [kpc] at last clean "
        "time -- NOTE: times differ per run, see labels"
    )
    ax.grid(alpha=0.25, axis="x")
    plt.tight_layout()
    out = os.path.join(OUT, "ranking_W1.png")
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compute", action="store_true")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--nproc", type=int, default=8)
    args = ap.parse_args()
    if args.compute:
        results = compute(args.nproc)
    else:
        with open(os.path.join(OUT, "cache.pkl"), "rb") as fh:
            results = pickle.load(fh)
    if args.plot or not args.compute:
        plot_profiles(results)
        plot_metrics(results)
        plot_psi_sweep(results)
        plot_ranking(results)


if __name__ == "__main__":
    main()
