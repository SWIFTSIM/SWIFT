#!/usr/bin/env python3
"""
Compare the results of the hydro test suite run by run_suite.sh for several
hydro schemes.

Usage: python3 plot_comparison.py <suite_output_dir> [test ...]

For every test, a figure <suite_output_dir>/plots/<test>.png is produced with
one column per scheme (raw particles + binned profile + analytic solution
where one exists) and a last column overlaying the binned profiles of all the
schemes. A summary table (steps, wall-clock time, energy conservation and
error norms w.r.t. the analytic solutions) is written to
<suite_output_dir>/plots/summary.md and summary.json.
"""

import glob
import json
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.special import gamma as Gamma  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))  # examples/HydroTests/riemannSolver.py

GAS_GAMMA = 5.0 / 3.0
COLORS = ["C0", "C3", "C2", "C1", "C4", "C5", "C6", "C7", "C8", "C9"]
ALL_TESTS = ["sod", "sedov", "noh", "gresho", "evrard", "kh", "square", "keplerian", "keplerian2d", "zeldovich",
             "zeldovich_glass", "zeldovich_pert", "blob"]  # "nfw" on request

scatter_props = dict(marker=".", s=1, alpha=0.15, rasterized=True, linewidths=0)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def snapshot_files(run_dir, basename):
    files = sorted(glob.glob(os.path.join(run_dir, f"{basename}_[0-9][0-9][0-9][0-9].hdf5")))
    return files


def load(filename):
    """Load the gas particle data of a snapshot into a dict."""
    d = {}
    with h5py.File(filename, "r") as f:
        d["boxsize"] = np.atleast_1d(f["/Header"].attrs["BoxSize"])
        if d["boxsize"].size == 1:
            d["boxsize"] = np.repeat(d["boxsize"], 3)
        d["time"] = float(np.atleast_1d(f["/Header"].attrs.get("Time", 0.0))[0])
        d["a"] = float(np.atleast_1d(f["/Header"].attrs.get("Scale-factor", 1.0))[0])
        d["z"] = float(np.atleast_1d(f["/Header"].attrs.get("Redshift", 0.0))[0])
        d["scheme"] = str(f["/HydroScheme"].attrs.get("Scheme", "")) if "/HydroScheme" in f else ""
        g = f["/PartType0"]
        d["pos"] = g["Coordinates"][:]
        d["vel"] = g["Velocities"][:]
        d["m"] = g["Masses"][:].reshape(-1)
        d["ids"] = g["ParticleIDs"][:].reshape(-1)
        # Fields absent from initial-condition files
        for key, name in [("rho", "Densities"), ("u", "InternalEnergies"), ("P", "Pressures"),
                          ("S", "Entropies"), ("h", "SmoothingLengths")]:
            d[key] = g[name][:].reshape(-1) if name in g else None
        if "/Units" in f:
            d["U_L"] = float(np.atleast_1d(f["/Units"].attrs["Unit length in cgs (U_L)"])[0])
            d["U_M"] = float(np.atleast_1d(f["/Units"].attrs["Unit mass in cgs (U_M)"])[0])
            d["U_t"] = float(np.atleast_1d(f["/Units"].attrs["Unit time in cgs (U_t)"])[0])
        if "/Cosmology" in f and "H0 [internal units]" in f["/Cosmology"].attrs:
            d["H0"] = float(np.atleast_1d(f["/Cosmology"].attrs["H0 [internal units]"])[0])
    return d


def binned(x, y, edges):
    """Mean and standard deviation of y in bins of x."""
    mean, _, _ = stats.binned_statistic(x, y, statistic="mean", bins=edges)
    sq, _, _ = stats.binned_statistic(x, y * y, statistic="mean", bins=edges)
    std = np.sqrt(np.maximum(sq - mean * mean, 0.0))
    return 0.5 * (edges[1:] + edges[:-1]), mean, std


def l1(x, y, x_ref, y_ref, mask, norm=None):
    """Relative L1 error of the particle values y(x) w.r.t. the reference curve."""
    ref = np.interp(x[mask], x_ref, y_ref)
    ok = np.isfinite(ref) & np.isfinite(y[mask])
    if norm is None:
        norm = np.mean(np.abs(ref[ok]))
    return float(np.mean(np.abs(y[mask][ok] - ref[ok])) / norm)


def run_stats(run_dir):
    """Number of steps, wall-clock time and energy conservation of a run."""
    out = {}
    ts = glob.glob(os.path.join(run_dir, "timesteps*.txt"))
    if ts:
        rows = [l.split() for l in open(ts[0]) if not l.startswith("#")]
        data = np.array([r for r in rows if len(r) == 15], dtype=float)
        out["steps"] = int(data.shape[0])
        out["wallclock_s"] = float(np.sum(data[:, 12]) / 1000.0)
        out["updates"] = float(np.sum(data[:, 7]))
    wt = os.path.join(run_dir, "walltime_s")
    if os.path.exists(wt):
        out["total_walltime_s"] = float(open(wt).read())
    st = os.path.join(run_dir, "statistics.txt")
    if os.path.exists(st):
        data = np.loadtxt(st, comments="#", ndmin=2)
        E = data[:, 13] + data[:, 14] + data[:, 15]
        cosmological = np.any(data[:, 2] != 1.0)
        out["energy_drift"] = float((E[-1] - E[0]) / abs(E[0])) if (E[0] != 0 and not cosmological) else float("nan")
        out["E_kin_end"] = float(data[-1, 13])
        out["E_int_end"] = float(data[-1, 14])
    # Fallback diagnostics of the MAGMA scheme (final snapshot)
    snaps = sorted(glob.glob(os.path.join(run_dir, "*_[0-9][0-9][0-9][0-9].hdf5")))
    if snaps:
        with h5py.File(snaps[-1], "r") as f:
            g = f["/PartType0"]
            if "FallbackFlags" in g:
                flags = g["FallbackFlags"][:]
                # Bits 1 and 2 switch the particle to base SPH; bit 4 only
                # records a regularised (planar/filamentary) C-matrix.
                out["frac_base_SPH_final"] = float(np.mean((flags & 3) != 0))
                out["frac_cond_fallback_final"] = float(np.mean((flags & 1) != 0))
                out["frac_regularised_final"] = float(np.mean((flags & 4) != 0))
            if "GradientFallbackPairs" in g:
                out["mean_pair_fallbacks_final"] = float(np.mean(g["GradientFallbackPairs"][:]))
    log = os.path.join(run_dir, "output.log")
    if os.path.exists(log):
        n_base = 0
        n_nongb = 0
        with open(log, errors="replace") as f:
            for line in f:
                if "will use base SPH" in line:
                    n_base += 1
                if "treated as having no neighbours" in line:
                    n_nongb += 1
        out["base_SPH_warnings"] = n_base
        out["no_neighbour_warnings"] = n_nongb
    return out


def setup_figure(schemes, quantities, extra_rows=0, width=3.6, height=2.8):
    ncol = len(schemes) + 1
    nrow = len(quantities) + extra_rows
    fig, axes = plt.subplots(nrow, ncol, figsize=(width * ncol, height * nrow), squeeze=False)
    return fig, axes


def profile_panels(axes, row, schemes, data, xkey_fn, ykey_fn, edges, exact=None, xlabel="",
                   ylabel="", xlim=None, ylim=None, logx=False, logy=False):
    """One row of panels: raw + binned per scheme, plus an overlay of the binned profiles."""
    over = axes[row, -1]
    for j, s in enumerate(schemes):
        ax = axes[row, j]
        d = data[s]
        if d is None:
            ax.text(0.5, 0.5, "missing", ha="center", transform=ax.transAxes)
            continue
        x = xkey_fn(d)
        y = ykey_fn(d)
        ax.scatter(x, y, color="0.4", **scatter_props)
        xb, yb, sb = binned(x, y, edges)
        ax.errorbar(xb, yb, yerr=sb, fmt=".", color=COLORS[j % len(COLORS)], ms=4, lw=1, zorder=3)
        over.plot(xb, yb, "-", color=COLORS[j % len(COLORS)], lw=1.3, label=s)
        if exact is not None:
            ax.plot(exact[0], exact[1], "-", color="k", lw=1, alpha=0.8, zorder=2)
        if row == 0:
            ax.set_title(s)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        if logx:
            ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
    if exact is not None:
        over.plot(exact[0], exact[1], "-", color="k", lw=1, alpha=0.8, label="exact / reference")
    if row == 0:
        over.set_title("binned overlay")
        over.legend(fontsize=7, loc="best")
    over.set_ylabel(ylabel)
    over.set_xlabel(xlabel)
    if xlim:
        over.set_xlim(xlim)
    if ylim:
        over.set_ylim(ylim)
    if logx:
        over.set_xscale("log")
    if logy:
        over.set_yscale("log")


def image_panels(axes, row, schemes, data, ckey, clim, title_fn=None, cmap="viridis", log=False):
    """One row of density maps (one per scheme); the last column is left for a curve."""
    for j, s in enumerate(schemes):
        ax = axes[row, j]
        d = data[s]
        if d is None:
            ax.text(0.5, 0.5, "missing", ha="center", transform=ax.transAxes)
            continue
        c = d[ckey]
        if log:
            c = np.log10(c)
        sc = ax.scatter(d["pos"][:, 0], d["pos"][:, 1], c=c, s=1.5, vmin=clim[0], vmax=clim[1],
                        cmap=cmap, linewidths=0, rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlim(0, d["boxsize"][0])
        ax.set_ylim(0, d["boxsize"][1])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{s}, " + (title_fn(d) if title_fn else f"t={d['time']:.2f}"), fontsize=9)
    plt.colorbar(sc, ax=axes[row, len(schemes) - 1], fraction=0.046)


def pick_snapshot(run_dir, basename, index=None, z=None, time=None):
    files = snapshot_files(run_dir, basename)
    if not files:
        return None
    if index is not None:
        return files[index] if index < len(files) else files[-1]
    if z is not None or time is not None:
        best, best_d = None, np.inf
        for f in files:
            with h5py.File(f, "r") as h:
                val = float(np.atleast_1d(h["/Header"].attrs["Redshift" if z is not None else "Time"])[0])
            dist = abs(val - (z if z is not None else time))
            if dist < best_d:
                best, best_d = f, dist
        return best
    return files[-1]


def load_all(root, schemes, test, basename, **kw):
    data = {}
    for s in schemes:
        run_dir = os.path.join(root, "runs", s, test)
        f = pick_snapshot(run_dir, basename, **kw)
        data[s] = load(f) if f else None
    return data


# ---------------------------------------------------------------------------
# Analytic solutions
# ---------------------------------------------------------------------------


def sedov_solution(t, E0, rho0, g, n=1000, nu=3):
    """Sedov-Taylor similarity solution (from examples/HydroTests/SedovBlast_3D)."""
    a = [0] * 8
    a[0] = 2.0 / (nu + 2)
    a[2] = (1 - g) / (2 * (g - 1) + nu)
    a[3] = nu / (2 * (g - 1) + nu)
    a[5] = 2 / (g - 2)
    a[6] = g / (2 * (g - 1) + nu)
    a[1] = (((nu + 2) * g) / (2.0 + nu * (g - 1.0))) * (
        (2.0 * nu * (2.0 - g)) / (g * (nu + 2.0) ** 2) - a[2])
    a[4] = a[1] * (nu + 2) / (2 - g)
    a[7] = (2 + nu * (g - 1)) * a[1] / (nu * (2 - g))

    v_min = 2.0 / ((nu + 2) * g)
    v_max = 4.0 / ((nu + 2) * (g + 1))
    v = v_min + np.arange(n) * (v_max - v_min) / (n - 1.0)

    beta = (nu + 2) * (g + 1) * np.array((
        0.25, (g / (g - 1)) * 0.5,
        -(2 + nu * (g - 1)) / 2.0 / ((nu + 2) * (g + 1) - 2 * (2 + nu * (g - 1))),
        -0.5 / (g - 1)), dtype=np.float64)
    beta = np.outer(beta, v)
    beta += (g + 1) * np.array((
        0.0, -1.0 / (g - 1), (nu + 2) / ((nu + 2) * (g + 1) - 2.0 * (2 + nu * (g - 1))),
        1.0 / (g - 1)), dtype=np.float64).reshape((4, 1))
    lbeta = np.log(beta)

    r = np.exp(-a[0] * lbeta[0] - a[2] * lbeta[1] - a[1] * lbeta[2])
    rho = ((g + 1.0) / (g - 1.0)) * np.exp(a[3] * lbeta[1] + a[5] * lbeta[3] + a[4] * lbeta[2])
    p = np.exp(nu * a[0] * lbeta[0] + (a[5] + 1) * lbeta[3] + (a[4] - 2 * a[1]) * lbeta[2])
    u = beta[0] * r * 4.0 / ((g + 1) * (nu + 2))
    p *= 8.0 / ((g + 1) * (nu + 2) * (nu + 2))
    u[0] = 0.0
    rho[0] = 0.0
    r[0] = 0.0
    p[0] = p[1]

    vol = (np.pi ** (nu / 2.0) / Gamma(nu / 2.0 + 1)) * np.power(r, nu)
    de = rho * u * u * 0.5 + p / (g - 1)
    q = np.inner(de[1:] + de[:-1], np.diff(vol)) * 0.5
    fac = (q * (t ** nu) * rho0 / E0) ** (-1.0 / (nu + 2))
    shock_speed = fac * (2.0 / (nu + 2))
    r_s = shock_speed * t * (nu + 2) / 2.0
    r *= fac * t
    u *= fac
    p *= fac * fac * rho0
    rho *= rho0
    return r, p, rho, u, r_s


def gresho_solution(r):
    P = np.where(r < 0.2, 5.0 + 12.5 * r ** 2,
                 np.where(r < 0.4, 9.0 + 12.5 * r ** 2 - 20.0 * r + 4.0 * np.log(np.maximum(r, 1e-10) / 0.2),
                          3.0 + 4.0 * np.log(2.0)))
    v = np.where(r < 0.2, 5.0 * r, np.where(r < 0.4, 2.0 - 5.0 * r, 0.0))
    return P, v


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_sod(root, schemes, out):
    from riemannSolver import RiemannSolver

    data = load_all(root, schemes, "sod", "sodShock", index=1)
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t = ref["time"]
    solver = RiemannSolver(GAS_GAMMA)
    x_s = np.linspace(-0.5, 0.5, 1000)
    rho_s, v_s, P_s, _ = solver.solve(1.0, 0.0, 1.0, 0.125, 0.0, 0.1, x_s / t)
    u_s = P_s / (rho_s * (GAS_GAMMA - 1.0))

    def xf(d):
        return d["pos"][:, 0] - 1.0

    edges = np.arange(-0.6, 0.6, 0.02)
    quantities = [("rho", lambda d: d["rho"], rho_s, r"$\rho$"),
                  ("v", lambda d: d["vel"][:, 0], v_s, r"$v_x$"),
                  ("P", lambda d: d["P"], P_s, r"$P$"),
                  ("u", lambda d: d["u"], u_s, r"$u$")]
    fig, axes = setup_figure(schemes, quantities)
    metrics = {}
    for i, (name, yf, ex, lab) in enumerate(quantities):
        profile_panels(axes, i, schemes, data, xf, yf, edges, exact=(x_s, ex), xlabel="x",
                       ylabel=lab, xlim=(-0.6, 0.6))
        for s in schemes:
            if data[s] is None:
                continue
            x = xf(data[s])
            mask = np.abs(x) < 0.5
            metrics.setdefault(s, {})[f"L1_{name}"] = l1(x, yf(data[s]), x_s, ex, mask)
    fig.suptitle(f"Sod shock (3D), t = {t:.3f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "sod.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_sedov(root, schemes, out):
    data = load_all(root, schemes, "sedov", "sedov", index=5)
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t = ref["time"]
    r_s, P_s, rho_s, v_s, r_shock = sedov_solution(t, 1.0, 1.0, GAS_GAMMA)
    r_s = np.append(r_s, [r_shock, r_shock * 1.5])
    rho_s = np.append(rho_s, [1.0, 1.0])
    P_s = np.append(P_s, [1e-6, 1e-6])
    v_s = np.append(v_s, [0.0, 0.0])
    with np.errstate(divide="ignore", invalid="ignore"):
        u_s = P_s / (rho_s * (GAS_GAMMA - 1.0))
    u_s[~np.isfinite(u_s)] = np.max(u_s[np.isfinite(u_s)])  # rho -> 0 at the centre

    def rf(d):
        c = 0.5 * d["boxsize"]
        return np.sqrt(np.sum((d["pos"] - c) ** 2, axis=1))

    def vr(d):
        c = 0.5 * d["boxsize"]
        dx = d["pos"] - c
        r = np.sqrt(np.sum(dx ** 2, axis=1))
        return np.sum(dx * d["vel"], axis=1) / np.maximum(r, 1e-10)

    edges = np.arange(0.0, 0.5, 0.01)
    quantities = [("rho", lambda d: d["rho"], rho_s, r"$\rho$"),
                  ("v", vr, v_s, r"$v_r$"),
                  ("P", lambda d: d["P"], P_s, r"$P$"),
                  ("u", lambda d: d["u"], u_s, r"$u$")]
    fig, axes = setup_figure(schemes, quantities)
    metrics = {}
    for i, (name, yf, ex, lab) in enumerate(quantities):
        profile_panels(axes, i, schemes, data, rf, yf, edges, exact=(r_s, ex), xlabel="r",
                       ylabel=lab, xlim=(0, 0.5), ylim=(0, None) if name != "u" else None,
                       logy=(name == "u"))
        for s in schemes:
            if data[s] is None:
                continue
            r = rf(data[s])
            mask = r < 1.3 * r_shock
            # normalised by the post-shock value of the exact solution
            metrics.setdefault(s, {})[f"L1_{name}"] = l1(r, yf(data[s]), r_s, ex, mask,
                                                       norm=np.max(np.abs(ex)))
    for s in schemes:
        if data[s] is not None:
            r = rf(data[s])
            xb, yb, _ = binned(r, data[s]["rho"], edges)
            metrics[s]["rho_peak_over_exact"] = float(np.nanmax(yb) / np.max(rho_s))
    fig.suptitle(f"Sedov blast (3D), t = {t:.3f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "sedov.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_noh(root, schemes, out):
    data = load_all(root, schemes, "noh", "noh", index=12)
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t = ref["time"]
    g = GAS_GAMMA
    r_s = np.linspace(1e-3, 1.0, 1000)
    rs = 0.5 * (g - 1) * t
    rho_s = np.where(r_s < rs, ((g + 1) / (g - 1)) ** 3, (1 + t / r_s) ** 2)
    v_s = np.where(r_s < rs, 0.0, -1.0)
    u_s = np.where(r_s < rs, 0.5, 1e-6 / (g - 1))

    def rf(d):
        c = 0.5 * d["boxsize"]
        return np.sqrt(np.sum((d["pos"] - c) ** 2, axis=1))

    def vr(d):
        c = 0.5 * d["boxsize"]
        dx = d["pos"] - c
        r = np.sqrt(np.sum(dx ** 2, axis=1))
        return np.sum(dx * d["vel"], axis=1) / np.maximum(r, 1e-10)

    edges = np.arange(0.0, 1.0, 0.02)
    quantities = [("rho", lambda d: d["rho"], rho_s, r"$\rho$"),
                  ("v", vr, v_s, r"$v_r$"),
                  ("u", lambda d: d["u"], u_s, r"$u$")]
    fig, axes = setup_figure(schemes, quantities)
    metrics = {}
    for i, (name, yf, ex, lab) in enumerate(quantities):
        profile_panels(axes, i, schemes, data, rf, yf, edges, exact=(r_s, ex), xlabel="r",
                       ylabel=lab, xlim=(0, 1.0))
        for s in schemes:
            if data[s] is None:
                continue
            r = rf(data[s])
            mask = r < 0.8
            metrics.setdefault(s, {})[f"L1_{name}"] = l1(r, yf(data[s]), r_s, ex, mask,
                                                       norm=np.max(np.abs(ex)))
    fig.suptitle(f"Noh implosion (3D), t = {t:.3f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "noh.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_gresho(root, schemes, out):
    data = load_all(root, schemes, "gresho", "gresho", index=10)
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t = ref["time"]
    r_s = np.linspace(0, 0.8, 500)
    P_s, v_s = gresho_solution(r_s)
    rho_s = np.ones_like(r_s)

    def rf(d):
        c = 0.5 * d["boxsize"]
        return np.sqrt((d["pos"][:, 0] - c[0]) ** 2 + (d["pos"][:, 1] - c[1]) ** 2)

    def vphi(d):
        c = 0.5 * d["boxsize"]
        x = d["pos"][:, 0] - c[0]
        y = d["pos"][:, 1] - c[1]
        r = np.sqrt(x * x + y * y)
        return (-y * d["vel"][:, 0] + x * d["vel"][:, 1]) / np.maximum(r, 1e-10)

    def vrad(d):
        c = 0.5 * d["boxsize"]
        x = d["pos"][:, 0] - c[0]
        y = d["pos"][:, 1] - c[1]
        r = np.sqrt(x * x + y * y)
        return (x * d["vel"][:, 0] + y * d["vel"][:, 1]) / np.maximum(r, 1e-10)

    edges = np.arange(0.0, 0.8, 0.02)
    quantities = [("vphi", vphi, v_s, r"$v_\phi$"),
                  ("vr", vrad, np.zeros_like(r_s), r"$v_r$"),
                  ("P", lambda d: d["P"], P_s, r"$P$"),
                  ("rho", lambda d: d["rho"], rho_s, r"$\rho$")]
    fig, axes = setup_figure(schemes, quantities)
    metrics = {}
    for i, (name, yf, ex, lab) in enumerate(quantities):
        profile_panels(axes, i, schemes, data, rf, yf, edges, exact=(r_s, ex), xlabel="r",
                       ylabel=lab, xlim=(0, 0.8))
        for s in schemes:
            if data[s] is None:
                continue
            r = rf(data[s])
            mask = r < 0.8
            metrics.setdefault(s, {})[f"L1_{name}"] = l1(r, yf(data[s]), r_s, ex, mask,
                                                       norm=1.0 if name in ("vphi", "vr") else None)
    fig.suptitle(f"Gresho vortex (2D), t = {t:.3f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "gresho.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_evrard(root, schemes, out):
    data = load_all(root, schemes, "evrard", "evrard", index=8)
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t = ref["time"]
    reffile = None
    for s in schemes:
        f = os.path.join(root, "runs", s, "evrard", "evrardCollapse3D_exact.txt")
        if os.path.exists(f):
            reffile = f
    exact = np.loadtxt(reffile) if reffile else None

    def rf(d):
        c = 0.5 * d["boxsize"]
        return np.sqrt(np.sum((d["pos"] - c) ** 2, axis=1))

    def vr(d):
        c = 0.5 * d["boxsize"]
        dx = d["pos"] - c
        r = np.sqrt(np.sum(dx ** 2, axis=1))
        return np.sum(dx * d["vel"], axis=1) / np.maximum(r, 1e-10)

    edges = np.logspace(-2.5, np.log10(2.0), 60)
    quantities = [("rho", lambda d: d["rho"], 1, r"$\rho$", True),
                  ("v", vr, 2, r"$v_r$", False),
                  ("P", lambda d: d["P"], 3, r"$P$", True),
                  ("S", lambda d: d["S"], None, r"$P/\rho^\gamma$", True)]
    fig, axes = setup_figure(schemes, quantities)
    metrics = {}
    for i, (name, yf, col, lab, logy) in enumerate(quantities):
        ex = None
        if exact is not None:
            if col is not None:
                ex = (exact[:, 0], exact[:, col])
            else:
                ex = (exact[:, 0], exact[:, 3] / exact[:, 1] ** GAS_GAMMA)
        profile_panels(axes, i, schemes, data, rf, yf, edges, exact=ex, xlabel="r",
                       ylabel=lab, xlim=(3e-3, 2.0), logx=True, logy=logy)
        for s in schemes:
            if data[s] is None or ex is None:
                continue
            r = rf(data[s])
            mask = (r > 0.01) & (r < 1.0)
            y = yf(data[s])
            if logy:
                metrics.setdefault(s, {})[f"L1_log_{name}"] = float(np.mean(np.abs(
                    np.log10(np.maximum(y[mask], 1e-30)) - np.log10(np.maximum(np.interp(r[mask], ex[0], ex[1]), 1e-30)))))
            else:
                metrics.setdefault(s, {})[f"L1_{name}"] = l1(r, y, ex[0], ex[1], mask)
    fig.suptitle(f"Evrard collapse (3D), t = {t:.3f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "evrard.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_kh(root, schemes, out):
    times = [1.5, 3.0, 4.5]
    fig, axes = plt.subplots(len(times), len(schemes) + 1,
                             figsize=(3.6 * (len(schemes) + 1), 3.4 * len(times)), squeeze=False)
    metrics = {}
    for i, t in enumerate(times):
        data = load_all(root, schemes, "kh", "kelvinHelmholtz", time=t)
        image_panels(axes, i, schemes, data, "rho", (0.8, 2.2), cmap="RdBu_r")
        axes[i, -1].axis("off")
    # Growth of the instability: rms of v_y versus time
    ax = axes[0, -1]
    ax.axis("on")
    for j, s in enumerate(schemes):
        files = snapshot_files(os.path.join(root, "runs", s, "kh"), "kelvinHelmholtz")
        tt, vy = [], []
        for f in files:
            d = load(f)
            tt.append(d["time"])
            vy.append(np.sqrt(np.mean(d["vel"][:, 1] ** 2)))
        if tt:
            ax.semilogy(tt, vy, "-o", ms=3, color=COLORS[j % len(COLORS)], label=s)
            metrics[s] = {"rms_vy_final": float(vy[-1]), "rms_vy_t1.5": float(np.interp(1.5, tt, vy))}
    ax.set_xlabel("t")
    ax.set_ylabel(r"rms $v_y$")
    ax.legend(fontsize=7)
    fig.suptitle("Kelvin-Helmholtz (2D), density")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "kh.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_square(root, schemes, out):
    times = [1.0, 4.0]
    fig, axes = plt.subplots(len(times) + 1, len(schemes) + 1,
                             figsize=(3.6 * (len(schemes) + 1), 3.4 * (len(times) + 1)), squeeze=False)
    metrics = {}
    for i, t in enumerate(times):
        data = load_all(root, schemes, "square", "square", time=t)
        image_panels(axes, i, schemes, data, "rho", (0.5, 4.5))
        for j, s in enumerate(schemes):
            axes[i, j].add_patch(plt.Rectangle((0.25, 0.25), 0.5, 0.5, fill=False, ec="k", lw=0.8))
        axes[i, -1].axis("off")
    # Density along a slice through the middle at the final time
    data = load_all(root, schemes, "square", "square", time=times[-1])
    ic = None
    for s in schemes:
        f = os.path.join(root, "runs", s, "square", "square.hdf5")
        if os.path.exists(f):
            ic = load(f)
            break

    def xf(d):
        return d["pos"][:, 0]

    edges = np.arange(0.0, 1.0, 0.01)
    row = len(times)
    for j, s in enumerate(schemes):
        d = data[s]
        ax = axes[row, j]
        if d is None:
            continue
        sl = np.abs(d["pos"][:, 1] - 0.5) < 0.03
        ax.scatter(xf(d)[sl], d["rho"][sl], color=COLORS[j % len(COLORS)], s=4)
        xb, yb, _ = binned(xf(d)[sl], d["rho"][sl], edges)
        axes[row, -1].plot(xb, yb, "-", color=COLORS[j % len(COLORS)], label=s)
        ax.set_xlabel("x (|y-0.5|<0.03)")
        ax.set_ylabel(r"$\rho$")
        ax.set_ylim(0, 5)
        if ic is not None:
            # target density from the initial positions
            order = np.argsort(ic["ids"])
            pos0 = ic["pos"][order][np.searchsorted(ic["ids"][order], d["ids"])]
            inside = (np.abs(pos0[:, 0] - 0.5) < 0.25) & (np.abs(pos0[:, 1] - 0.5) < 0.25)
            target = np.where(inside, 4.0, 1.0)
            metrics[s] = {"L1_rho_vs_initial": float(np.mean(np.abs(d["rho"] / target - 1.0)))}
    axes[row, -1].plot([0, 0.25, 0.25, 0.75, 0.75, 1], [1, 1, 4, 4, 1, 1], "k-", lw=0.8, label="initial")
    axes[row, -1].set_xlabel("x")
    axes[row, -1].set_ylabel(r"$\rho$")
    axes[row, -1].legend(fontsize=7)
    fig.suptitle("Square test (2D), density")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "square.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_keplerian(root, schemes, out, test="keplerian"):
    times = [0.0, 10.0, 25.0, 50.0]
    fig, axes = plt.subplots(len(schemes) + 1, len(times),
                             figsize=(3.4 * len(times), 3.4 * (len(schemes) + 1)), squeeze=False)
    metrics = {}
    clim = (-1.0, 0.3)
    edges = np.linspace(0, 5, 51)
    areas = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    for i, t in enumerate(times):
        data = load_all(root, schemes, test, "keplerian_ring", time=t)
        for j, s in enumerate(schemes):
            ax = axes[j, i]
            d = data[s]
            if d is None:
                ax.text(0.5, 0.5, "missing", ha="center", transform=ax.transAxes)
                continue
            c = 0.5 * d["boxsize"]
            x = d["pos"][:, 0] - c[0]
            y = d["pos"][:, 1] - c[1]
            if i == 0 and j == 0:
                # colour range from the initial ring (same for all panels)
                clim = np.log10(np.percentile(d["rho"], [1, 99.5]))
            sc = ax.scatter(x, y, c=np.log10(d["rho"]), s=1.5, cmap="viridis", linewidths=0,
                            rasterized=True, vmin=clim[0], vmax=clim[1])
            ax.set_aspect("equal")
            ax.set_xlim(-5, 5)
            ax.set_ylim(-5, 5)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{s}, t={d['time']:.1f}", fontsize=9)
            r = np.sqrt(x * x + y * y)
            mass, _ = np.histogram(r, bins=edges, weights=d["m"])
            sigma = mass / areas
            axes[-1, i].plot(0.5 * (edges[1:] + edges[:-1]), sigma, "-", color=COLORS[j % len(COLORS)], label=s)
            # Fraction of the mass that left the initial ring region (1 < r < 3.5)
            metrics.setdefault(s, {})[f"mass_frac_outside_ring_t{t:.0f}"] = float(
                np.sum(d["m"][(r < 1.0) | (r > 3.5)]) / np.sum(d["m"]))
        axes[-1, i].set_xlabel("r")
        axes[-1, i].set_ylabel(r"$\Sigma$")
        axes[-1, i].set_title(f"surface density, t={t:.0f}", fontsize=9)
        axes[-1, i].legend(fontsize=7)
    fig.suptitle(f"Keplerian ring ({'2D code' if test.endswith('2d') else '3D planar'}), log density")
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{test}.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_zeldovich(root, schemes, out, test="zeldovich"):
    T_i, z_c, z_i = 100.0, 1.0, 100.0
    mH, kB = 1.6737236e-27, 1.38064852e-23  # SI
    rows = [("z = 3 (before the caustic)", dict(z=3.0)), ("final", dict())]
    quantities = ["rho", "v", "T", "S"]
    fig, axes = plt.subplots(len(rows) * len(quantities), len(schemes) + 1,
                             figsize=(3.6 * (len(schemes) + 1), 2.6 * len(rows) * len(quantities)),
                             squeeze=False)
    metrics = {}
    # Reference entropy from the first snapshot (for the adiabaticity check)
    S0 = {}
    for s in schemes:
        f = pick_snapshot(os.path.join(root, "runs", s, test), "zeldovichPancake", index=0)
        if f:
            d0 = load(f)
            order = np.argsort(d0["ids"])
            S0[s] = (d0["ids"][order], d0["S"][order])
    for ri, (label, sel) in enumerate(rows):
        data = load_all(root, schemes, test, "zeldovichPancake", **sel)
        ref = next((d for d in data.values() if d), None)
        if ref is None:
            continue
        z = ref["z"]
        a = ref["a"]
        L = ref["boxsize"][0]
        rho_0 = np.sum(ref["m"]) / L ** 3
        H0 = ref["H0"]
        # Analytic (Zel'dovich) solution, valid before the caustic only
        q = np.linspace(-0.5 * L, 0.5 * L, 512)
        k = 2.0 * np.pi / L
        zfac = (1.0 + z_c) / (1.0 + z)
        x_s = q - zfac * np.sin(k * q) / k
        rho_s = rho_0 / (1.0 - zfac * np.cos(k * q))
        v_s = -H0 * (1.0 + z_c) / np.sqrt(1.0 + z) * np.sin(k * q) / k
        T_s = T_i * (((1.0 + z) / (1.0 + z_i)) ** 3 * rho_s / rho_0) ** (2.0 / 3.0)
        u_unit = (0.01 * ref["U_L"] / ref["U_t"]) ** 2  # SI

        def xf(d):
            return d["pos"][:, 0] - 0.5 * L

        def temp(d):
            u_phys = d["u"] * d["a"] ** (-3.0 * (GAS_GAMMA - 1.0)) * u_unit
            return u_phys * (GAS_GAMMA - 1.0) * mH / kB

        pre_caustic = z > z_c
        qs = [("rho", lambda d: d["rho"] / rho_0, (x_s, rho_s / rho_0) if pre_caustic else None, r"$\rho/\rho_0$ (comoving)", True),
              ("v", lambda d: d["vel"][:, 0], (x_s, v_s) if pre_caustic else None, r"$v_x$ (peculiar)", False),
              ("T", temp, (x_s, T_s) if pre_caustic else None, r"$T$ [K]", True),
              ("S", lambda d: d["S"], None, r"$A = P/\rho^\gamma$ (comoving)", True)]
        edges = np.linspace(-0.5 * L, 0.5 * L, 129)
        for qi, (name, yf, ex, lab, logy) in enumerate(qs):
            row = ri * len(quantities) + qi
            profile_panels(axes, row, schemes, data, xf, yf, edges, exact=ex, xlabel="x [Mpc, comoving]",
                           ylabel=lab, xlim=(-0.5 * L, 0.5 * L), logy=logy)
            if qi == 0:
                for ax in axes[row, :]:
                    ax.set_title(ax.get_title() + f"\n{label}, z = {z:.2f}", fontsize=9)
            for s in schemes:
                d = data[s]
                if d is None or ex is None:
                    continue
                x = xf(d)
                # order the analytic solution in x for interpolation (single-valued before caustic)
                o = np.argsort(ex[0])
                mask = np.ones_like(x, dtype=bool)
                metrics.setdefault(s, {})[f"L1_{name}_z{z:.1f}"] = l1(x, yf(d), ex[0][o], ex[1][o], mask)
        # Adiabaticity: entropy change since the ICs (should be 0 before the caustic)
        for s in schemes:
            d = data[s]
            if d is None or s not in S0:
                continue
            ids0, s0 = S0[s]
            idx = np.searchsorted(ids0, d["ids"])
            dS = d["S"] / s0[idx] - 1.0
            key = f"z{z:.1f}"
            metrics.setdefault(s, {})[f"max_dS_over_S_{key}"] = float(np.max(np.abs(dS)))
            metrics[s][f"median_dS_over_S_{key}"] = float(np.median(np.abs(dS)))
            metrics[s][f"frac_dS_gt_1pc_{key}"] = float(np.mean(np.abs(dS) > 0.01))
    variant = {"zeldovich_glass": " -- glass IC", "zeldovich_pert": " -- lattice symmetry broken"}.get(test, "")
    fig.suptitle("Zel'dovich pancake (3D, cosmological)" + variant)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{test}.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_blob(root, schemes, out):
    """Blob in a Mach 2.7 wind: projected density maps and the surviving blob
    mass fraction (Agertz et al. 2007 criterion: rho > 0.64 rho_blob and
    u < 0.9 u_ambient) as a function of time in units of tau_KH."""
    # The example's wind has speed 1 (its sound speed is set for Mach 2.7)
    tau_kh = (1.0 + 10.0) * 0.1 / (1.0 * np.sqrt(10.0))
    times = [1.0, 2.0, 4.0]  # in tau_KH
    fig, axes = plt.subplots(len(times), len(schemes) + 1,
                             figsize=(3.6 * (len(schemes) + 1), 2.0 * len(times) + 1.0), squeeze=False)
    metrics = {}
    for i, t in enumerate(times):
        data = load_all(root, schemes, "blob", "blob", time=t * tau_kh)
        for j, s in enumerate(schemes):
            ax = axes[i, j]
            d = data[s]
            if d is None:
                ax.text(0.5, 0.5, "missing", ha="center", transform=ax.transAxes)
                continue
            # Mass-weighted projection along z
            # Mass-weighted projection along z of the region around the blob and its wake
            xmax = 0.6 * d["boxsize"][0]
            H, xe, ye = np.histogram2d(d["pos"][:, 0], d["pos"][:, 1], bins=(300, 125),
                                       range=[[0, xmax], [0, d["boxsize"][1]]], weights=d["m"])
            ax.imshow(np.log10(H.T + 1e-10 * H.max()), origin="lower", extent=(0, xmax, 0, d["boxsize"][1]),
                      cmap="viridis", vmin=np.log10(H.max()) - 1.5, vmax=np.log10(H.max()))
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{s}, t={d['time'] / tau_kh:.1f} " + r"$\tau_{\rm KH}$", fontsize=9)
        axes[i, -1].axis("off")
    # Surviving blob mass fraction
    ax = axes[0, -1]
    ax.axis("on")
    for j, s in enumerate(schemes):
        run_dir = os.path.join(root, "runs", s, "blob")
        ic_file = os.path.join(run_dir, "blob.hdf5")
        files = snapshot_files(run_dir, "blob")
        if not files or not os.path.exists(ic_file):
            continue
        ic = load(ic_file)
        centre = np.array([0.5, 0.5, 0.5])
        blob_ids = np.sort(ic["ids"][np.sum((ic["pos"] - centre) ** 2, axis=1) < 0.1 ** 2])
        # Reference blob density and ambient internal energy from the first snapshot
        # (the example's units give rho_bg = 1024, not 1)
        d0 = load(files[0])
        in_blob0 = np.isin(d0["ids"], blob_ids)
        rho_blob = float(np.median(d0["rho"][in_blob0]))
        u_bg = float(np.median(d0["u"][~in_blob0]))
        tt, frac = [], []
        for f in files:
            d = load(f)
            in_blob = np.isin(d["ids"], blob_ids)
            cold_dense = (d["rho"] > 0.64 * rho_blob) & (d["u"] < 0.9 * u_bg)
            tt.append(d["time"] / tau_kh)
            frac.append(float(np.sum(in_blob & cold_dense) / max(np.sum(in_blob), 1)))
        ax.plot(tt, frac, "-o", ms=3, color=COLORS[j % len(COLORS)], label=s)
        metrics[s] = {"blob_mass_frac_2tau": float(np.interp(2.0, tt, frac)),
                      "blob_mass_frac_final": float(frac[-1]), "t_final_tauKH": float(tt[-1])}
    ax.set_xlabel(r"$t / \tau_{\rm KH}$")
    ax.set_ylabel("blob mass fraction")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7)
    fig.suptitle("Blob test (3D), projected density")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "blob.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_nfw(root, schemes, out):
    """Gas in hydrostatic equilibrium in an NFW potential: density profile at
    the end compared with the initial one."""
    fig, axes = plt.subplots(2, len(schemes) + 1, figsize=(3.6 * (len(schemes) + 1), 6.0), squeeze=False)
    metrics = {}
    edges = np.logspace(np.log10(0.5), np.log10(100.0), 36)
    rb = np.sqrt(edges[1:] * edges[:-1])

    def profile(d):
        c = 0.5 * d["boxsize"]
        r = np.sqrt(np.sum((d["pos"] - c) ** 2, axis=1))
        _, rho, _ = binned(r, d["rho"], edges)
        return rb, rho

    for j, s in enumerate(schemes):
        run_dir = os.path.join(root, "runs", s, "nfw")
        ic_file = os.path.join(run_dir, "nfw.hdf5")
        f = pick_snapshot(run_dir, "snapshot")
        if f is None or not os.path.exists(ic_file):
            axes[0, j].text(0.5, 0.5, "missing", ha="center", transform=axes[0, j].transAxes)
            continue
        ic, d = load(ic_file), load(f)
        r0, rho0 = profile(ic)
        r1, rho1 = profile(d)
        ax = axes[0, j]
        ax.loglog(r0, rho0, "k-", lw=1, label="initial")
        ax.loglog(r1, rho1, "-", color=COLORS[j % len(COLORS)], lw=1.3, label=f"t={d['time']:.2f}")
        ax.set_title(s)
        ax.set_ylabel(r"$\rho$")
        ax.set_xlabel("r")
        ax.legend(fontsize=7)
        ax = axes[1, j]
        ax.semilogx(r1, rho1 / rho0, "-", color=COLORS[j % len(COLORS)], lw=1.3)
        ax.axhline(1.0, color="k", lw=0.8)
        ax.set_ylim(0.5, 1.5)
        ax.set_xlabel("r")
        ax.set_ylabel(r"$\rho / \rho_{\rm initial}$")
        axes[0, -1].loglog(r1, rho1, "-", color=COLORS[j % len(COLORS)], lw=1.3, label=s)
        axes[1, -1].semilogx(r1, rho1 / rho0, "-", color=COLORS[j % len(COLORS)], lw=1.3, label=s)
        sel = (rb > 1.0) & (rb < 50.0) & np.isfinite(rho1) & np.isfinite(rho0)
        metrics[s] = {"L1_logrho_1-50": float(np.mean(np.abs(np.log10(rho1[sel] / rho0[sel])))),
                      "rho_ratio_r<2": float(np.nanmean((rho1 / rho0)[rb < 2.0]))}
        if j == 0:
            axes[0, -1].loglog(r0, rho0, "k-", lw=1, label="initial")
    axes[0, -1].set_title("overlay")
    axes[0, -1].legend(fontsize=7)
    axes[1, -1].axhline(1.0, color="k", lw=0.8)
    axes[1, -1].set_ylim(0.5, 1.5)
    axes[1, -1].legend(fontsize=7)
    fig.suptitle("NFW hydrostatic halo (3D), gas density profile")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "nfw.png"), dpi=150)
    plt.close(fig)
    return metrics


TESTS = {"sod": test_sod, "sedov": test_sedov, "noh": test_noh, "gresho": test_gresho,
         "evrard": test_evrard, "kh": test_kh, "square": test_square,
         "keplerian": test_keplerian,
         "keplerian2d": lambda root, schemes, out: test_keplerian(root, schemes, out, test="keplerian2d"),
         "zeldovich": test_zeldovich,
         "zeldovich_glass": lambda root, schemes, out: test_zeldovich(root, schemes, out, test="zeldovich_glass"),
         "zeldovich_pert": lambda root, schemes, out: test_zeldovich(root, schemes, out, test="zeldovich_pert"),
         "blob": test_blob, "nfw": test_nfw}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    root = sys.argv[1]
    tests = sys.argv[2:] or ALL_TESTS
    schemes = sorted(os.listdir(os.path.join(root, "runs")))
    out = os.path.join(root, "plots")
    os.makedirs(out, exist_ok=True)
    summary = {}
    for t in tests:
        if t not in TESTS:
            print(f"Unknown test {t}")
            continue
        present = [s for s in schemes if os.path.exists(os.path.join(root, "runs", s, t, "DONE"))]
        if not present:
            print(f"{t}: no finished run")
            continue
        print(f"{t}: {present}")
        try:
            metrics = TESTS[t](root, present, out)
        except Exception as e:  # keep going with the other tests
            import traceback
            traceback.print_exc()
            print(f"{t}: plotting failed: {e}")
            metrics = {}
        summary[t] = {}
        for s in present:
            summary[t][s] = run_stats(os.path.join(root, "runs", s, t))
            summary[t][s].update(metrics.get(s, {}))

    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)

    lines = ["# Hydro test suite summary", ""]
    for t, per_scheme in summary.items():
        keys = sorted({k for m in per_scheme.values() for k in m})
        lines.append(f"## {t}")
        lines.append("")
        lines.append("| quantity | " + " | ".join(per_scheme) + " |")
        lines.append("|---|" + "---|" * len(per_scheme))
        for k in keys:
            vals = []
            for s in per_scheme:
                v = per_scheme[s].get(k, "")
                if isinstance(v, float):
                    vals.append(f"{v:.4g}")
                else:
                    vals.append(str(v))
            lines.append(f"| {k} | " + " | ".join(vals) + " |")
        lines.append("")
    with open(os.path.join(out, "summary.md"), "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
