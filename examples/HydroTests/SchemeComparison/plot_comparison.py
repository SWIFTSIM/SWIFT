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
# Fixed colours for the schemes that appear on many figures, so that a scheme
# keeps its colour whatever the set of schemes on a given figure; other
# schemes (parameter-study variants) take the next free colour by position.
SCHEME_COLORS = {"sphenix": "C2", "mi3": "C3", "mi2": "C0", "mi1": "C1"}


def scheme_color(scheme, j):
    return SCHEME_COLORS.get(scheme, COLORS[j % len(COLORS)])

ALL_TESTS = ["sod", "sedov", "noh", "gresho", "evrard", "kh", "square", "square_mass", "keplerian", "keplerian2d",
             "keplerian_cd10", "zeldovich",
             "zeldovich_glass", "zeldovich_pert", "blob", "taylorgreen", "shearing", "blast", "rt",
             "khgr", "ibw", "soundwave", "wavedamping", "thermal", "gresho_boost", "kh8",
             "sod_weak", "sod_strong"]  # "nfw" on request

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


# Schemes drawn only in the panels that overlay all the curves (growth of
# the KH rms velocity, the square's density cut, the blob's mass fraction),
# without a column of their own in the image-based figures; set with
# --overlay-only=a,b. Their metrics are computed as for the others.
OVERLAY_ONLY = set()


def column_schemes(schemes):
    return [s for s in schemes if s not in OVERLAY_ONLY]


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
        ax.errorbar(xb, yb, yerr=sb, fmt=".", color=scheme_color(s, j), ms=4, lw=1, zorder=3)
        over.plot(xb, yb, "-", color=scheme_color(s, j), lw=1.3, label=s)
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


def test_sod(root, schemes, out, test="sod"):
    """Sod shock; test = sod_weak / sod_strong are the same tube with another right-state pressure
    (read from the IC header). A second figure zooms on the contact discontinuity."""
    from riemannSolver import RiemannSolver

    data = load_all(root, schemes, test, "sodShock", index=1)
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t = ref["time"]
    states = [1.0, 0.0, 1.0, 0.125, 0.0, 0.1]
    for s in schemes:
        att = ic_attrs(os.path.join(root, "runs", s, test), "sodShock.hdf5", "SodShockStates")
        if att is not None:
            states = [float(v) for v in att[0]]
            break
    solver = RiemannSolver(GAS_GAMMA)
    x_s = np.linspace(-1.0, 1.0, 8000)  # the whole tube (the hot right state of the weak tube sends its shock past 0.5)
    rho_s, v_s, P_s, _ = solver.solve(*states, x_s / t)
    u_s = P_s / (rho_s * (GAS_GAMMA - 1.0))
    # Contact: the density jump without a pressure jump; shock: the largest pressure jump to its right
    jumps = np.abs(np.diff(rho_s))
    dP = np.abs(np.diff(P_s))
    contact_candidates = np.where((jumps > 0.05 * jumps.max()) & (dP < 1e-3 * P_s.max()))[0]
    i_c = contact_candidates[np.argmax(jumps[contact_candidates])] if len(contact_candidates) else int(np.argmax(jumps))
    x_c = 0.5 * (x_s[i_c] + x_s[i_c + 1])
    right = np.arange(i_c + 1, len(dP))
    i_sh = right[np.argmax(dP[right])] if len(right) else i_c
    x_sh = 0.5 * (x_s[i_sh] + x_s[i_sh + 1])
    w = float(min(0.12, 0.6 * (x_sh - x_c))) if x_sh > x_c + 0.02 else 0.12
    P_star = float(P_s[i_c])
    print(f"{test}: contact at x = {x_c:.3f}, shock at x = {x_sh:.3f}, window +-{w:.3f}")
    rho_star_L, rho_star_R = float(rho_s[i_c - 5]), float(rho_s[i_c + 6])

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
    label = {"sod": "Sod shock", "sod_weak": "weak Sod shock", "sod_strong": "strong Sod shock"}.get(test, test)
    fig.suptitle(f"{label} (3D, P_R = {states[5]:g}), t = {t:.3f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{test}.png"), dpi=150)
    plt.close(fig)

    # Zoom on the contact discontinuity: profiles within +-w of x_c, L1 there, pressure blip and
    # the width of the density transition (10-90 per cent of the star-state jump)
    edges_c = np.arange(x_c - w, x_c + w + 1e-9, 0.01)
    fig, axes = setup_figure(schemes, quantities)
    in_w = np.abs(x_s - x_c) < w
    for i, (name, yf, ex, lab) in enumerate(quantities):
        lo, hi = float(np.min(ex[in_w])), float(np.max(ex[in_w]))
        pad = max(0.3 * (hi - lo), 0.6 * abs(hi))  # the contact's pressure blip can reach +-70 per cent
        profile_panels(axes, i, schemes, data, xf, yf, edges_c, exact=(x_s, ex), xlabel="x",
                       ylabel=lab, xlim=(x_c - w, x_c + w), ylim=(lo - pad, hi + pad))
        for s in schemes:
            if data[s] is None:
                continue
            x = xf(data[s])
            mask = np.abs(x - x_c) < w
            metrics[s][f"contact_L1_{name}"] = l1(x, yf(data[s]), x_s, ex, mask)
    for s in schemes:
        d = data[s]
        if d is None:
            continue
        x = xf(d)
        mask = np.abs(x - x_c) < w
        xb, Pb, _ = binned(x[mask], d["P"][mask], edges_c)
        metrics[s]["contact_P_blip"] = float(np.nanmax(np.abs(Pb / P_star - 1.0)))
        xb, rb, _ = binned(x[mask], d["rho"][mask], edges_c)
        jump = rho_star_L - rho_star_R
        hi = np.where(rb < rho_star_L - 0.1 * jump)[0]
        lo = np.where(rb < rho_star_R + 0.1 * jump)[0]
        metrics[s]["contact_width_rho"] = float(xb[lo[0]] - xb[hi[0]]) if len(hi) and len(lo) and lo[0] >= hi[0] else float("nan")
        axes[0, [j for j, ss in enumerate(schemes) if ss == s][0]].axvline(x_c, color="k", lw=0.5, ls=":")
    fig.suptitle(f"{label}: the contact discontinuity (x_c = {x_c:.3f}, window +-{w:.2f}), t = {t:.3f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{test}_contact.png"), dpi=150)
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


def test_gresho(root, schemes, out, test="gresho"):
    """Gresho vortex; test="gresho_boost" is the same vortex advected by a bulk velocity (Galilean
    invariance): the bulk velocity is the mass-weighted mean velocity of the box, the vortex centre
    has moved by it, and the metrics are also compared with the static run of the same scheme."""
    data = load_all(root, schemes, test, "gresho", index=10)
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t = ref["time"]
    r_s = np.linspace(0, 0.8, 500)
    P_s, v_s = gresho_solution(r_s)
    rho_s = np.ones_like(r_s)

    def frame(d):
        """Position relative to the (advected) vortex centre and the velocity in its rest frame."""
        if "frame" not in d:
            L = d["boxsize"]
            vb = np.sum(d["m"][:, None] * d["vel"], axis=0) / np.sum(d["m"])
            c = np.mod(0.5 * L + vb * d["time"], np.where(L > 0, L, 1.0))
            dx = d["pos"] - c
            dx -= L * np.round(dx / np.where(L > 0, L, 1.0))
            d["frame"] = (dx, d["vel"] - vb)
        return d["frame"]

    def rf(d):
        dx, _ = frame(d)
        return np.sqrt(dx[:, 0] ** 2 + dx[:, 1] ** 2)

    def vphi(d):
        dx, v = frame(d)
        r = np.sqrt(dx[:, 0] ** 2 + dx[:, 1] ** 2)
        return (-dx[:, 1] * v[:, 0] + dx[:, 0] * v[:, 1]) / np.maximum(r, 1e-10)

    def vrad(d):
        dx, v = frame(d)
        r = np.sqrt(dx[:, 0] ** 2 + dx[:, 1] ** 2)
        return (dx[:, 0] * v[:, 0] + dx[:, 1] * v[:, 1]) / np.maximum(r, 1e-10)

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
    if test == "gresho_boost":
        # Galilean invariance: the same error as the static vortex of the same scheme?
        static = load_all(root, schemes, "gresho", "gresho", index=10)
        for s in schemes:
            if data[s] is None or static[s] is None:
                continue
            r = rf(static[s])
            L1_static = l1(r, vphi(static[s]), r_s, v_s, r < 0.8, norm=1.0)
            metrics[s]["L1_vphi_static"] = L1_static
            metrics[s]["dL1_vphi"] = metrics[s]["L1_vphi"] - L1_static
            vb = np.sum(data[s]["m"][:, None] * data[s]["vel"], axis=0) / np.sum(data[s]["m"])
            metrics[s]["boost"] = float(np.sqrt(np.sum(vb**2)))
        title = f"Gresho vortex (2D) advected at |v| = {metrics[schemes[0]].get('boost', 0):.1f}, t = {t:.3f} (vortex frame)"
    else:
        title = f"Gresho vortex (2D), t = {t:.3f}"
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{test}.png"), dpi=150)
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
        # The exact profiles run to densities and pressures 1e-15 below the
        # gas's in the outer vacuum: keep the y range to where the particles are.
        ylim = None
        if logy:
            lo, hi = [], []
            for s in schemes:
                if data[s] is None:
                    continue
                r = rf(data[s])
                y = yf(data[s])[(r > 3e-3) & (r < 2.0)]
                y = y[y > 0]
                if y.size:
                    lo.append(np.percentile(y, 0.5))
                    hi.append(np.percentile(y, 99.9))
            if lo:
                ylim = (0.3 * min(lo), 3.0 * max(hi))
        profile_panels(axes, i, schemes, data, rf, yf, edges, exact=ex, xlabel="r",
                       ylabel=lab, xlim=(3e-3, 2.0), ylim=ylim, logx=True, logy=logy)
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


def test_kh(root, schemes, out, test="kh"):
    """Kelvin-Helmholtz; test="kh8" is the same set-up at a density contrast of 8 (longer run)."""
    last = max((load(f)["time"] for s in schemes
                for f in snapshot_files(os.path.join(root, "runs", s, test), "kelvinHelmholtz")[-1:]), default=4.5)
    times = [last / 3.0, 2.0 * last / 3.0, last]
    cols = column_schemes(schemes)
    fig, axes = plt.subplots(len(times), len(cols) + 1,
                             figsize=(3.6 * (len(cols) + 1), 3.4 * len(times)), squeeze=False)
    metrics = {}
    rho_max, contrast = 2.0, 2
    for i, t in enumerate(times):
        data = load_all(root, cols, test, "kelvinHelmholtz", time=t)
        ref = next((d for d in data.values() if d), None)
        if ref is not None and i == 0:
            rho_max = float(np.percentile(ref["rho"], 99.5))
            # Density contrast of the initial conditions (first snapshot of the reference run)
            first = next((snapshot_files(os.path.join(root, "runs", s, test), "kelvinHelmholtz")[:1]
                          for s in cols if data[s] is not None), [])
            rho_ic = load(first[0])["rho"] if first else ref["rho"]
            contrast = round(float(np.percentile(rho_ic, 90)) / float(np.percentile(rho_ic, 10)))
        log = rho_max > 3.0
        clim = (np.log10(0.8), np.log10(rho_max + 0.2)) if log else (0.8, rho_max + 0.2)
        image_panels(axes, i, cols, data, "rho", clim, cmap="RdBu_r", log=log)
        axes[i, -1].axis("off")
    # Growth of the instability: rms of v_y versus time (all schemes)
    ax = axes[0, -1]
    ax.axis("on")
    for j, s in enumerate(schemes):
        files = snapshot_files(os.path.join(root, "runs", s, test), "kelvinHelmholtz")
        tt, vy = [], []
        for f in files:
            d = load(f)
            tt.append(d["time"])
            vy.append(np.sqrt(np.mean(d["vel"][:, 1] ** 2)))
        if tt:
            ax.semilogy(tt, vy, "-o", ms=3, color=scheme_color(s, j), label=s)
            metrics[s] = {"rms_vy_final": float(vy[-1]), f"rms_vy_t{times[0]:.3g}": float(np.interp(times[0], tt, vy))}
    ax.set_xlabel("t")
    ax.set_ylabel(r"rms $v_y$")
    ax.legend(fontsize=7)
    fig.suptitle(f"Kelvin-Helmholtz (2D, density contrast {contrast}), " + ("log10 density" if rho_max > 3.0 else "density"))
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{test}.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_square(root, schemes, out, test="square"):
    times = [1.0, 4.0]
    cols = column_schemes(schemes)
    fig, axes = plt.subplots(len(times) + 1, len(cols) + 1,
                             figsize=(3.6 * (len(cols) + 1), 3.4 * (len(times) + 1)), squeeze=False)
    metrics = {}
    for i, t in enumerate(times):
        data = load_all(root, cols, test, "square", time=t)
        image_panels(axes, i, cols, data, "rho", (0.5, 4.5))
        for j, s in enumerate(cols):
            axes[i, j].add_patch(plt.Rectangle((0.25, 0.25), 0.5, 0.5, fill=False, ec="k", lw=0.8))
        axes[i, -1].axis("off")
    # Density along a slice through the middle at the final time
    data = load_all(root, schemes, test, "square", time=times[-1])
    ic = None
    for s in schemes:
        f = os.path.join(root, "runs", s, test, "square.hdf5")
        if os.path.exists(f):
            ic = load(f)
            break

    def xf(d):
        return d["pos"][:, 0]

    edges = np.arange(0.0, 1.0, 0.01)
    row = len(times)
    for j, s in enumerate(schemes):
        d = data[s]
        if d is None:
            continue
        sl = np.abs(d["pos"][:, 1] - 0.5) < 0.03
        xb, yb, _ = binned(xf(d)[sl], d["rho"][sl], edges)
        axes[row, -1].plot(xb, yb, "-", color=scheme_color(s, j), label=s)
        if s in cols:
            ax = axes[row, cols.index(s)]
            ax.scatter(xf(d)[sl], d["rho"][sl], color=scheme_color(s, j), s=4)
            ax.set_xlabel("x (|y-0.5|<0.03)")
            ax.set_ylabel(r"$\rho$")
            ax.set_ylim(0, 5)
        if ic is not None:
            # target density from the initial positions
            order = np.argsort(ic["ids"])
            pos0 = ic["pos"][order][np.searchsorted(ic["ids"][order], d["ids"])]
            inside = (np.abs(pos0[:, 0] - 0.5) < 0.25) & (np.abs(pos0[:, 1] - 0.5) < 0.25)
            target = np.where(inside, 4.0, 1.0)
            # Shape: dense gas that left the initial square (rotation into a
            # diamond, fingers); the per-particle L1 is blind to it.
            now_inside = (np.abs(d["pos"][:, 0] - 0.5) < 0.25) & (np.abs(d["pos"][:, 1] - 0.5) < 0.25)
            dense = d["rho"] > 2.5
            metrics[s] = {"L1_rho_vs_initial": float(np.mean(np.abs(d["rho"] / target - 1.0))),
                          "frac_dense_outside_square": float(np.sum(dense & ~now_inside) / max(np.sum(dense), 1))}
    axes[row, -1].plot([0, 0.25, 0.25, 0.75, 0.75, 1], [1, 1, 4, 4, 1, 1], "k-", lw=0.8, label="initial")
    axes[row, -1].set_xlabel("x")
    axes[row, -1].set_ylabel(r"$\rho$")
    axes[row, -1].legend(fontsize=7)
    fig.suptitle("Square test (2D), density" + (" -- equal spacing, 4:1 masses" if test == "square_mass" else ""))
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{test}.png"), dpi=150)
    plt.close(fig)
    return metrics


def ring_positions_figure(runs, times, out_path, half_width=2.5, title=None):
    """Particle positions of the Keplerian ring at several times, in the style
    of Fig. 8 of Cullen & Dehnen (2010): one column per run, one row per time,
    black dots on white, centred on the point mass.

    runs: list of (label, run_dir) pairs; times: list of output times.
    """
    fig, axes = plt.subplots(len(times), len(runs), figsize=(2.6 * len(runs), 2.6 * len(times)),
                             squeeze=False, sharex=True, sharey=True)
    for j, (label, run_dir) in enumerate(runs):
        axes[0, j].set_title(label, fontsize=9)
        for i, t in enumerate(times):
            ax = axes[i, j]
            f = pick_snapshot(run_dir, "keplerian_ring", time=t)
            if f is None:
                ax.text(0.5, 0.5, "no output", transform=ax.transAxes, ha="center", va="center")
                continue
            with h5py.File(f, "r") as h:
                pos = h["/PartType0/Coordinates"][:]
                box = np.atleast_1d(h["/Header"].attrs["BoxSize"])
                t_snap = float(np.atleast_1d(h["/Header"].attrs["Time"])[0])
            x = pos[:, 0] - 0.5 * box[0]
            y = pos[:, 1] - 0.5 * box[min(1, len(box) - 1)]
            ax.plot(x, y, ",", color="k", rasterized=True)
            ax.text(0.04, 0.95, f"t={t_snap:.0f}", transform=ax.transAxes, ha="left", va="top", fontsize=9)
            ax.set_xlim(-half_width, half_width)
            ax.set_ylim(-half_width, half_width)
            ax.set_aspect("equal")
            if i == len(times) - 1:
                ax.set_xlabel("x")
            if j == 0:
                ax.set_ylabel("y")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


# Geometry of the ring tests: the Hopkins (2015) disc occupies 0.5 < r < 2 and
# the Cullen & Dehnen (2010) ring is a Gaussian of centre 10 and width 2.5
# (annulus = +-2 widths). "plunge" radii define the plunger metrics, "half"
# the window of the position plots.
RING_GEOMETRY = {
    "keplerian": dict(r_in=0.5, r_out=2.0, plunge=(0.25, 0.1), half=2.5,
                      times=[0.0, 10.0, 25.0, 50.0], pos_times=[5.0, 10.0, 20.0, 30.0]),
    "keplerian2d": dict(r_in=0.5, r_out=2.0, plunge=(0.25, 0.1), half=2.5,
                        times=[0.0, 10.0, 25.0, 50.0], pos_times=[5.0, 10.0, 20.0, 30.0]),
    # the Gaussian's inner tail reaches r ~ 1.7 in the IC: plunge radii below it
    "keplerian_cd10": dict(r_in=5.0, r_out=15.0, plunge=(1.5, 0.5), half=16.0,
                           times=[0.0, 10.0, 20.0, 40.0], pos_times=[10.0, 15.0, 25.0, 40.0]),
}


def test_keplerian(root, schemes, out, test="keplerian"):
    geo = RING_GEOMETRY[test]
    times = geo["times"]
    # Initial mass in the ring annulus, for the ring-mass metric
    ring_mass_0 = {}
    for s in schemes:
        f0 = pick_snapshot(os.path.join(root, "runs", s, test), "keplerian_ring", index=0)
        if f0:
            d0 = load(f0)
            r0 = np.hypot(d0["pos"][:, 0] - 0.5 * d0["boxsize"][0], d0["pos"][:, 1] - 0.5 * d0["boxsize"][1])
            ring_mass_0[s] = float(np.sum(d0["m"][(r0 > geo["r_in"]) & (r0 < geo["r_out"])]))
    fig, axes = plt.subplots(len(schemes) + 1, len(times),
                             figsize=(3.4 * len(times), 3.4 * (len(schemes) + 1)), squeeze=False)
    metrics = {}
    clim = (-1.0, 0.3)
    edges = np.linspace(0, 2.0 * geo["half"], 51)
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
            ax.set_xlim(-2.0 * geo["half"], 2.0 * geo["half"])
            ax.set_ylim(-2.0 * geo["half"], 2.0 * geo["half"])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{s}, t={d['time']:.1f}", fontsize=9)
            r = np.sqrt(x * x + y * y)
            mass, _ = np.histogram(r, bins=edges, weights=d["m"])
            sigma = mass / areas
            axes[-1, i].plot(0.5 * (edges[1:] + edges[:-1]), sigma, "-", color=scheme_color(s, j), label=s)
            # Fraction of the ring's mass that left the initial ring annulus
            # (0.5 < r < 2.0 in the example's IC; the gas beyond r = 2 is a
            # 100:1 lower-density background filling the periodic box)
            in_ring = (r > geo["r_in"]) & (r < geo["r_out"])
            metrics.setdefault(s, {})[f"mass_frac_outside_ring_t{t:.0f}"] = float(
                1.0 - np.sum(d["m"][in_ring]) / ring_mass_0[s]) if ring_mass_0.get(s) else float("nan")
        axes[-1, i].set_xlabel("r")
        axes[-1, i].set_ylabel(r"$\Sigma$")
        axes[-1, i].set_title(f"surface density, t={t:.0f}", fontsize=9)
        axes[-1, i].legend(fontsize=7)
    # Plungers: particles that fall well inside the inner edge of the ring
    # (r < 0.25, the inner edge is at r = 0.5) towards the point mass. They
    # end the useful life of the test (they set the time step and seed the
    # chaotic divergence), so the time before the first one appears is a
    # figure of merit: the longer, the better. Resolution: the snapshot cadence.
    # inside the inner edge; deep in the potential (see RING_GEOMETRY)
    r_plunge = {"": geo["plunge"][0], f"_r{geo['plunge'][1]:g}": geo["plunge"][1]}
    for s in schemes:
        files = sorted(glob.glob(os.path.join(root, "runs", s, test, "keplerian_ring_[0-9]*.hdf5")))
        t_first = {k: None for k in r_plunge}
        n_final = {k: 0 for k in r_plunge}
        mass_final = {k: 0.0 for k in r_plunge}
        t_last, r_min_final, m_tot = 0.0, np.inf, 1.0
        for f in files:
            with h5py.File(f, "r") as h:
                t_last = float(np.atleast_1d(h["/Header"].attrs["Time"])[0])
                pos = h["/PartType0/Coordinates"][:]
                mass = h["/PartType0/Masses"][:].reshape(-1)
                box = np.atleast_1d(h["/Header"].attrs["BoxSize"])
            r = np.hypot(pos[:, 0] - 0.5 * box[0], pos[:, 1] - 0.5 * box[min(1, len(box) - 1)])
            r_min_final = float(r.min())
            m_tot = float(mass.sum())
            for k, rp in r_plunge.items():
                sel = r < rp
                n_final[k] = int(np.sum(sel))
                mass_final[k] = float(mass[sel].sum())
                if n_final[k] > 0 and t_first[k] is None:
                    t_first[k] = t_last
        if files:
            m = metrics.setdefault(s, {})
            for k in r_plunge:
                # No plunger by the end: report the end time (a lower bound)
                m[f"t_first_plunger{k}"] = float(t_first[k] if t_first[k] is not None else t_last)
                m[f"plunger_free_to_end{k}"] = float(t_first[k] is None)
                m[f"n_plungers_final{k}"] = n_final[k]
                # Mass fraction, comparable between resolutions (the particle
                # masses follow the ring's density profile)
                m[f"plunger_mass_frac_final{k}"] = mass_final[k] / m_tot
            m["r_min_final"] = r_min_final
    fig.suptitle(f"Keplerian ring ({'2D code' if test.endswith('2d') else '3D planar'}), log density")
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{test}.png"), dpi=150)
    plt.close(fig)
    # Particle positions, Cullen & Dehnen (2010) Fig. 8 style
    runs = [(s, os.path.join(root, "runs", s, test)) for s in schemes
            if glob.glob(os.path.join(root, "runs", s, test, "keplerian_ring_[0-9]*.hdf5"))]
    if runs:
        ring_positions_figure(runs, geo["pos_times"], os.path.join(out, f"{test}_positions.png"), half_width=geo["half"],
                              title=f"Keplerian ring ({test}), particle positions")
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
    cols = column_schemes(schemes)
    fig, axes = plt.subplots(len(times), len(cols) + 1,
                             figsize=(3.6 * (len(cols) + 1), 2.0 * len(times) + 1.0), squeeze=False)
    metrics = {}
    for i, t in enumerate(times):
        data = load_all(root, cols, "blob", "blob", time=t * tau_kh)
        for j, s in enumerate(cols):
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
        ax.plot(tt, frac, "-o", ms=3, color=scheme_color(s, j), label=s)
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
        ax.loglog(r1, rho1, "-", color=scheme_color(s, j), lw=1.3, label=f"t={d['time']:.2f}")
        ax.set_title(s)
        ax.set_ylabel(r"$\rho$")
        ax.set_xlabel("r")
        ax.legend(fontsize=7)
        ax = axes[1, j]
        ax.semilogx(r1, rho1 / rho0, "-", color=scheme_color(s, j), lw=1.3)
        ax.axhline(1.0, color="k", lw=0.8)
        ax.set_ylim(0.5, 1.5)
        ax.set_xlabel("r")
        ax.set_ylabel(r"$\rho / \rho_{\rm initial}$")
        axes[0, -1].loglog(r1, rho1, "-", color=scheme_color(s, j), lw=1.3, label=s)
        axes[1, -1].semilogx(r1, rho1 / rho0, "-", color=scheme_color(s, j), lw=1.3, label=s)
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


# ---------------------------------------------------------------------------
# Tests added 2026-10-12 (viscosity-oriented: Taylor-Green, shearing shock,
# strong blast wave, Rayleigh-Taylor, KH growth rate, interacting blast waves)
# ---------------------------------------------------------------------------


def ic_arrays(run_dir, ic_name, *datasets):
    """Arrays of the initial-condition file linked into a run directory (None if absent)."""
    f = os.path.join(run_dir, ic_name)
    if not os.path.exists(f):
        return None
    with h5py.File(f, "r") as h:
        arrays = [h["/PartType0/" + d][:] for d in datasets]
    return [a.reshape(-1) if a.ndim == 2 and a.shape[1] == 1 else a for a in arrays]


def energy_history(run_dir):
    """Time and kinetic energy from statistics.txt (None if absent)."""
    st = os.path.join(run_dir, "statistics.txt")
    if not os.path.exists(st):
        return None
    data = np.loadtxt(st, comments="#", ndmin=2)
    return data[:, 1], data[:, 13]


def snapshot_series(run_dir, basename, fn):
    """Apply fn(d) to every snapshot of a run; returns (times, values)."""
    tt, vv = [], []
    for f in snapshot_files(run_dir, basename):
        d = load(f)
        tt.append(d["time"])
        vv.append(fn(d))
    return np.array(tt), np.array(vv)


def growth_rate(t, a, t0, t1):
    """Exponential growth rate of a(t) fitted on ln a between t0 and t1 (nan if too few points)."""
    ok = (t >= t0) & (t <= t1) & (a > 0)
    if np.sum(ok) < 3:
        return float("nan")
    return float(np.polyfit(t[ok], np.log(a[ok]), 1)[0])


def test_taylorgreen(root, schemes, out):
    """Stationary Taylor-Green vortex: kinetic-energy decay (effective viscosity) and the error of the final field."""
    v0, k = 0.1, 2.0 * np.pi
    data = load_all(root, schemes, "taylorgreen", "taylorGreen")
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t_end = ref["time"]
    fig, axes = plt.subplots(2, len(schemes) + 1, figsize=(3.6 * (len(schemes) + 1), 3.2 * 2), squeeze=False)
    metrics = {}
    axE = axes[0, -1]
    for j, s in enumerate(schemes):
        run_dir = os.path.join(root, "runs", s, "taylorgreen")
        hist = energy_history(run_dir)
        if hist is not None:
            t, E = hist
            ok = E > 0
            slope = np.polyfit(t[ok], np.log(E[ok] / E[0]), 1)[0]
            nu = -slope / (16.0 * np.pi**2)
            metrics.setdefault(s, {})["nu_eff"] = float(nu)
            metrics[s]["E_kin_ratio_end"] = float(E[-1] / E[0])
            axE.semilogy(t, E / E[0], "-", color=scheme_color(s, j), label=f"{s}: nu_eff = {nu:.1e}")
        d = data[s]
        if d is None:
            continue
        x, y = d["pos"][:, 0], d["pos"][:, 1]
        vx_ex = v0 * np.sin(k * x) * np.cos(k * y)
        vy_ex = -v0 * np.cos(k * x) * np.sin(k * y)
        dv = np.sqrt((d["vel"][:, 0] - vx_ex) ** 2 + (d["vel"][:, 1] - vy_ex) ** 2)
        metrics.setdefault(s, {})["L1_v"] = float(np.mean(dv) / v0)
        metrics[s]["L1_rho"] = float(np.mean(np.abs(d["rho"] - 1.0)))
        # Velocity magnitude map and the v_x profile along y = 0 (exact: v0 sin kx)
        d["vmag"] = np.sqrt(d["vel"][:, 0] ** 2 + d["vel"][:, 1] ** 2) / v0
        strip = np.minimum(y, 1.0 - y) < 0.02
        ax = axes[1, j]
        ax.scatter(x[strip], d["vel"][strip, 0] / v0, color="0.4", **scatter_props)
        xs = np.linspace(0, 1, 200)
        ax.plot(xs, np.sin(k * xs), "k-", lw=1)
        xb, vb, sb = binned(x[strip], d["vel"][strip, 0] / v0, np.linspace(0, 1, 26))
        ax.errorbar(xb, vb, yerr=sb, fmt=".", color=scheme_color(s, j), ms=4, lw=1, zorder=3)
        axes[1, -1].plot(xb, vb, "-", color=scheme_color(s, j), lw=1.3, label=s)
        ax.set_xlabel("x (|y| < 0.02)")
        ax.set_ylabel(r"$v_x / v_0$")
        ax.set_ylim(-1.2, 1.2)
    image_panels(axes, 0, schemes, data, "vmag", (0, 1), cmap="viridis")
    axE.axis("on")
    axE.set_xlabel("t")
    axE.set_ylabel(r"$E_{\rm kin}(t) / E_{\rm kin}(0)$")
    axE.legend(fontsize=6)
    axes[1, -1].plot(xs, np.sin(k * xs), "k-", lw=1, label="exact")
    axes[1, -1].set_xlabel("x")
    axes[1, -1].set_ylabel(r"$v_x / v_0$")
    axes[1, -1].legend(fontsize=7)
    fig.suptitle(f"Taylor-Green vortex (2D, Mach 0.1), t = {t_end:.2f}: |v| / v0 and kinetic-energy decay")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "taylorgreen.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_shearing(root, schemes, out):
    """Shearing shock: profiles against the exact 1D solution, post-shock noise and the diffusion of v_y."""
    from riemannSolver import RiemannSolver

    data = load_all(root, schemes, "shearing", "shearingShock")
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t = ref["time"]
    L = float(ref["boxsize"][0])
    xc = 0.5 * L
    # Stream speed, shear amplitude and pressure from the initial conditions
    first = next(s for s in schemes if data[s] is not None)
    ic = ic_arrays(os.path.join(root, "runs", first, "shearing"), "shearingShock.hdf5",
                   "Velocities", "ParticleIDs", "InternalEnergy")
    if ic is None:
        return {}
    v_ic, ids_ic, u_ic = ic
    dv = float(np.max(v_ic[:, 0]))
    s_amp = float(np.max(np.abs(v_ic[:, 1])))
    P0 = (GAS_GAMMA - 1.0) * float(u_ic[0])
    order_ic = np.argsort(ids_ic)
    vy0_sorted = v_ic[order_ic, 1]
    ids_sorted = ids_ic[order_ic]

    solver = RiemannSolver(GAS_GAMMA)
    x_s = np.linspace(0, L, 2000)
    rho_c, v_c, P_c, _ = solver.solve(1.0, dv, P0, 1.0, -dv, P0, (x_s - xc) / t)
    x_edge = np.where(x_s < xc, x_s, x_s - L)
    rho_e, v_e, P_e, _ = solver.solve(1.0, -dv, P0, 1.0, dv, P0, x_edge / t)
    central = np.abs(x_s - xc) < 0.25 * L
    rho_s = np.where(central, rho_c, rho_e)
    v_s = np.where(central, v_c, v_e)
    P_s = np.where(central, P_c, P_e)
    # Post-shock region: where the exact density is above half the jump
    shocked = central & (rho_s > 0.5 * (rho_s.max() + 1.0))
    x_shock = float(np.max(np.abs(x_s[shocked] - xc))) if np.any(shocked) else 0.0

    def xf(d):
        return d["pos"][:, 0]

    edges = np.linspace(0, L, 161)
    quantities = [("rho", lambda d: d["rho"], rho_s, r"$\rho$"),
                  ("vx", lambda d: d["vel"][:, 0], v_s, r"$v_x$"),
                  ("P", lambda d: d["P"], P_s, r"$P$")]
    fig, axes = setup_figure(schemes, quantities, extra_rows=1)
    metrics = {}
    for i, (name, yf, ex, lab) in enumerate(quantities):
        profile_panels(axes, i, schemes, data, xf, yf, edges, exact=(x_s, ex), xlabel="x", ylabel=lab, xlim=(0, L))
        for s in schemes:
            if data[s] is None:
                continue
            x = xf(data[s])
            metrics.setdefault(s, {})[f"L1_{name}"] = l1(x, yf(data[s]), x_s, ex, np.ones_like(x, dtype=bool))
    row = len(quantities)
    for j, s in enumerate(schemes):
        d = data[s]
        ax = axes[row, j]
        if d is None:
            continue
        x = xf(d)
        # Lagrangian change of v_y (exact: none) and post-shock noise of v_x
        vy0 = vy0_sorted[np.searchsorted(ids_sorted, d["ids"])]
        dvy = d["vel"][:, 1] - vy0
        metrics[s]["L1_vy_lagrangian"] = float(np.mean(np.abs(dvy)) / s_amp)
        post = np.abs(x - xc) < 0.8 * x_shock
        metrics[s]["vx_noise_postshock"] = float(np.std(d["vel"][post, 0]) / dv) if np.any(post) else float("nan")
        metrics[s]["x_shock_exact"] = x_shock
        ax.scatter(x, dvy / s_amp, color="0.4", **scatter_props)
        xb, yb, sb = binned(x, dvy / s_amp, edges)
        ax.errorbar(xb, yb, yerr=sb, fmt=".", color=scheme_color(s, j), ms=4, lw=1, zorder=3)
        axes[row, -1].plot(xb, sb, "-", color=scheme_color(s, j), lw=1.3, label=s)
        ax.axhline(0, color="k", lw=1)
        ax.set_xlim(0, L)
        ax.set_ylim(-0.05, 0.05)
        ax.set_xlabel("x")
        ax.set_ylabel(r"$(v_y - v_{y,0}) / s$")
    axes[row, -1].set_xlabel("x")
    axes[row, -1].set_ylabel(r"rms$(v_y - v_{y,0}) / s$ per bin")
    axes[row, -1].set_yscale("log")
    axes[row, -1].legend(fontsize=7)
    fig.suptitle(f"Shearing shock (3D), t = {t:.3f}, streams at Mach {dv:g}, shear amplitude {s_amp:g}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "shearing.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_blast(root, schemes, out):
    """Strong blast wave: profiles against the exact Riemann solution (P jump 1e4, uniform density)."""
    from riemannSolver import RiemannSolver

    data = load_all(root, schemes, "blast", "blastWave")
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t = ref["time"]
    L = float(ref["boxsize"][0])
    first = next(s for s in schemes if data[s] is not None)
    ic = ic_arrays(os.path.join(root, "runs", first, "blast"), "blastWave.hdf5", "InternalEnergy")
    if ic is None:
        return {}
    P_L, P_R = (GAS_GAMMA - 1.0) * float(ic[0].max()), (GAS_GAMMA - 1.0) * float(ic[0].min())
    # The periodic boundary carries the same jump: only the middle half of the
    # box, |x - L/2| < L/4, is free of the waves from it until the end time.
    solver = RiemannSolver(GAS_GAMMA)
    x_s = np.linspace(0.25 * L, 0.75 * L, 4000)
    rho_s, v_s, P_s, _ = solver.solve(1.0, 0.0, P_L, 1.0, 0.0, P_R, (x_s - 0.5 * L) / t)
    u_s = P_s / (rho_s * (GAS_GAMMA - 1.0))

    def xf(d):
        return d["pos"][:, 0]

    edges = np.linspace(0.25 * L, 0.75 * L, 101)
    quantities = [("rho", lambda d: d["rho"], rho_s, r"$\rho$", False),
                  ("v", lambda d: d["vel"][:, 0], v_s, r"$v_x$", False),
                  ("P", lambda d: d["P"], P_s, r"$P$", True),
                  ("u", lambda d: d["u"], u_s, r"$u$", True)]
    fig, axes = setup_figure(schemes, quantities)
    metrics = {}
    for i, (name, yf, ex, lab, logy) in enumerate(quantities):
        profile_panels(axes, i, schemes, data, xf, yf, edges, exact=(x_s, ex), xlabel="x", ylabel=lab,
                       xlim=(0.25 * L, 0.75 * L), logy=logy)
        for s in schemes:
            if data[s] is None:
                continue
            x = xf(data[s])
            metrics.setdefault(s, {})[f"L1_{name}"] = l1(x, yf(data[s]), x_s, ex, np.abs(x - 0.5 * L) < 0.25 * L)
    # Particle interpenetration: velocity noise behind the shock
    shocked = (rho_s > 0.5 * (rho_s.max() + 1.0)) & (x_s > 0.5 * L)
    if np.any(shocked):
        x0, x1 = float(x_s[shocked].min()), float(x_s[shocked].max())
        for s in schemes:
            if data[s] is None:
                continue
            x = xf(data[s])
            post = (x > x0 + 0.2 * (x1 - x0)) & (x < x1 - 0.2 * (x1 - x0))
            vref = float(np.interp(0.5 * (x0 + x1), x_s, v_s))
            metrics[s]["v_noise_postshock"] = float(np.std(data[s]["vel"][post, 0]) / vref) if np.any(post) else float("nan")
    fig.suptitle(f"Strong blast wave (3D), t = {t:.4f}, P_L / P_R = {P_L / P_R:g}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "blast.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_rt(root, schemes, out):
    """Rayleigh-Taylor: growth of the seeded mode against linear theory, bubble and spike heights."""
    gamma, g, k = 1.4, 0.5, 4.0 * np.pi  # the seed is cos(4 pi x): wavelength 1/2
    y_if = 0.75
    times = [1.0, 2.0, 3.0, 4.0]
    fig, axes = plt.subplots(len(times), len(schemes) + 1,
                             figsize=(3.0 * (len(schemes) + 1), 4.2 * len(times)), squeeze=False)
    metrics = {}
    sigma_th = None
    for i, tt in enumerate(times):
        data = load_all(root, schemes, "rt", "rayleigh_taylor", time=tt)
        image_panels(axes, i, schemes, data, "rho", (0.8, 2.2), cmap="RdBu_r")
        axes[i, -1].axis("off")
    ax = axes[0, -1]
    ax.axis("on")
    for j, s in enumerate(schemes):
        run_dir = os.path.join(root, "runs", s, "rt")
        ic = ic_arrays(run_dir, "rayleigh_taylor.hdf5", "Coordinates", "ParticleIDs")
        if ic is None:
            continue
        pos0, ids0 = ic
        order = np.argsort(ids0)
        y0_sorted = pos0[order, 1]
        ids_sorted = ids0[order]
        if sigma_th is None:
            rho_top, rho_bot = 2.0, 1.0
            A = (rho_top - rho_bot) / (rho_top + rho_bot)
            sigma_th = float(np.sqrt(A * g * k))

        def analyse(d):
            y0 = y0_sorted[np.searchsorted(ids_sorted, d["ids"])]
            free = (y0 > 0.1) & (y0 < 1.4)
            heavy = free & (y0 > y_if)
            light = free & (y0 < y_if)
            x, y = d["pos"][:, 0], d["pos"][:, 1]
            band = free & (np.abs(y0 - y_if) < 0.2)
            # amplitude of the seeded mode in v_y (projection on cos 4 pi x)
            amp = 2.0 * np.abs(np.mean(d["vel"][band, 1] * np.cos(k * x[band])))
            bubble = float(np.percentile(y[light], 99.5) - y_if)
            spike = float(y_if - np.percentile(y[heavy], 0.5))
            return amp, bubble, spike

        t, vals = snapshot_series(run_dir, "rayleigh_taylor", analyse)
        if len(t) == 0:
            continue
        amp, bubble, spike = vals[:, 0], vals[:, 1], vals[:, 2]
        # The seed relaxes until t ~ 0.5; the mode then grows cleanly until the
        # mushrooms form (t ~ 2.5). sigma_theory is the incompressible
        # sharp-interface value sqrt(A g k); the stratified, smoothed
        # interface of the example grows more slowly, so the ratio is a
        # comparison between schemes rather than an absolute error.
        sigma = growth_rate(t, amp, 0.5, 2.0)
        m = metrics.setdefault(s, {})
        m["sigma_fit"] = sigma
        m["sigma_theory"] = sigma_th
        m["sigma_ratio"] = float(sigma / sigma_th)
        for tq in (2.0, 4.0):
            if t[-1] >= tq - 0.05:
                m[f"bubble_height_t{tq:g}"] = float(np.interp(tq, t, bubble))
                m[f"spike_depth_t{tq:g}"] = float(np.interp(tq, t, spike))
        m["mixing_width_end"] = float(bubble[-1] + spike[-1])
        ax.semilogy(t, amp, "-", color=scheme_color(s, j), label=f"{s}: {sigma / sigma_th:.2f} sigma_th")
        axes[1, -1].axis("on")
        axes[1, -1].plot(t, bubble, "-", color=scheme_color(s, j), label=s)
        axes[1, -1].plot(t, spike, "--", color=scheme_color(s, j))
    if sigma_th is not None:
        tt = np.linspace(0, 1.2, 20)
        ax.semilogy(tt, 0.02 * np.exp(sigma_th * tt), "k:",
                    label=rf"$e^{{\sigma t}}$, $\sigma$ = {sigma_th:.2f}")
    ax.set_xlabel("t")
    ax.set_ylabel("mode amplitude of $v_y$")
    ax.legend(fontsize=6)
    axes[1, -1].set_xlabel("t")
    axes[1, -1].set_ylabel("bubble height (solid), spike depth (dashed)")
    axes[1, -1].legend(fontsize=6)
    fig.suptitle("Rayleigh-Taylor (2D), density")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "rt.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_khgr(root, schemes, out):
    """KH growth rate: E_y(t) against the analytic growth rate of the piecewise-linear shear layer."""
    U0, d_ramp, k = 0.5, 0.0317, 4.0 * np.pi
    # Rayleigh (1880): growth rate of a shear layer of half-thickness d with a linear velocity ramp
    arg = np.exp(-4.0 * k * d_ramp) - (2.0 * k * d_ramp - 1.0) ** 2
    sigma_th = float(U0 / (2.0 * d_ramp) * np.sqrt(max(arg, 0.0)))
    times = [1.0, 2.0]
    fig, axes = plt.subplots(len(times), len(schemes) + 1,
                             figsize=(3.4 * (len(schemes) + 1), 3.4 * len(times)), squeeze=False)
    metrics = {}
    for i, tt in enumerate(times):
        data = load_all(root, schemes, "khgr", "kelvinHelmholtzGrowthRate", time=tt)
        for d in data.values():
            if d is not None:
                d["vy"] = d["vel"][:, 1]
        image_panels(axes, i, schemes, data, "vy", (-0.3, 0.3), cmap="RdBu_r")
        axes[i, -1].axis("off")
    ax = axes[0, -1]
    ax.axis("on")

    def mode_amplitude(d):
        """Amplitude of the seeded v_y mode (McNally et al. 2012, eq. 10), averaged over the two interfaces."""
        x, y, vy = d["pos"][:, 0], d["pos"][:, 1], d["vel"][:, 1]
        amps = []
        for y_if in (0.25, 0.75):
            w = np.exp(-k * np.abs(y - y_if))
            si = np.sum(vy * np.sin(k * x) * w)
            ci = np.sum(vy * np.cos(k * x) * w)
            amps.append(2.0 * np.sqrt(si * si + ci * ci) / np.sum(w))
        return float(np.mean(amps))

    series = {}
    for j, s in enumerate(schemes):
        run_dir = os.path.join(root, "runs", s, "khgr")
        t, vals = snapshot_series(run_dir, "kelvinHelmholtzGrowthRate",
                                  lambda d: (mode_amplitude(d), float(np.sum(0.5 * d["m"] * d["vel"][:, 1] ** 2))))
        if len(t) == 0:
            continue
        M, Ey = vals[:, 0], vals[:, 1]
        series[s] = (t, M)
        # linear phase: the mode grows from the seed until it reaches ~ 0.1 of the shear
        t_lin = t[M < 0.1 * U0]
        t1 = float(min(t_lin[-1], 1.0)) if len(t_lin) else 0.6
        sigma = growth_rate(t, M, 0.1, t1)
        m = metrics.setdefault(s, {})
        m["sigma_fit"] = sigma
        m["sigma_theory"] = sigma_th
        m["sigma_ratio"] = float(sigma / sigma_th)
        m["mode_amp_t0.5"] = float(np.interp(0.5, t, M))
        m["mode_amp_max"] = float(np.max(M))
        m["t_mode_amp_max"] = float(t[np.argmax(M)])
        # Incoherent v_y generated by the shear flow itself (gradient noise), early on
        m["E_y_noise_t0.2"] = float(np.interp(0.2, t, Ey))
        ax.semilogy(t, M, "-", color=scheme_color(s, j), label=f"{s}: {sigma / sigma_th:.2f} sigma_th")
    tt = np.linspace(0, 0.8, 20)
    if series:
        M0 = float(np.mean([np.interp(0.1, *ser) for ser in series.values()]))
        ax.semilogy(tt + 0.1, M0 * np.exp(sigma_th * tt), "k:", label=rf"$e^{{\sigma t}}$, $\sigma$ = {sigma_th:.2f}")
    ax.set_xlabel("t")
    ax.set_ylabel("mode amplitude of $v_y$")
    ax.legend(fontsize=6)
    fig.suptitle("Kelvin-Helmholtz growth rate (2D, uniform density): $v_y$")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "khgr.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_ibw(root, schemes, out):
    """Interacting blast waves (1D, gamma = 7/5) against the high-resolution reference solution."""
    gamma = 1.4
    data = load_all(root, schemes, "ibw", "interactingBlastWaves")
    ref = next((d for d in data.values() if d), None)
    if ref is None:
        return {}
    t = ref["time"]
    first = next(s for s in schemes if data[s] is not None)
    ref_file = os.path.join(root, "runs", first, "ibw", "interactingBlastWaves1D_exact.txt")
    if not os.path.exists(ref_file):
        return {}
    exact = np.loadtxt(ref_file)
    x_s, rho_s, v_s, P_s = exact[:, 0], exact[:, 1], exact[:, 2], exact[:, 3]
    u_s = P_s / (rho_s * (gamma - 1.0))

    # The periodic box of length 2 holds the problem on [0, 1] and its mirror
    # image: fold both halves onto the reference domain.
    def xf(d):
        x = d["pos"][:, 0]
        return np.minimum(x, 2.0 - x)

    edges = np.linspace(0, 1, 201)
    quantities = [("rho", lambda d: d["rho"], rho_s, r"$\rho$", False),
                  # the mirror image has the opposite velocity
                  ("v", lambda d: d["vel"][:, 0] * np.where(d["pos"][:, 0] < 1.0, 1.0, -1.0), v_s, r"$v_x$", False),
                  ("P", lambda d: d["P"], P_s, r"$P$", True),
                  ("u", lambda d: d["u"], u_s, r"$u$", True)]
    fig, axes = setup_figure(schemes, quantities)
    metrics = {}
    for i, (name, yf, ex, lab, logy) in enumerate(quantities):
        profile_panels(axes, i, schemes, data, xf, yf, edges, exact=(x_s, ex), xlabel="min(x, 2 - x)", ylabel=lab,
                       xlim=(0.4, 1.0) if name != "u" else (0, 1), logy=logy)
        for s in schemes:
            if data[s] is None:
                continue
            x = xf(data[s])
            metrics.setdefault(s, {})[f"L1_{name}"] = l1(x, yf(data[s]), x_s, ex, np.ones_like(x, dtype=bool))
    fig.suptitle(f"Interacting blast waves (1D), t = {t:.3f}; the two mirrored halves of the box are folded together")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "ibw.png"), dpi=150)
    plt.close(fig)
    return metrics


def ic_attrs(run_dir, ic_name, *attrs):
    """Header attributes of the initial-condition file linked into a run directory (None if absent)."""
    f = os.path.join(run_dir, ic_name)
    if not os.path.exists(f):
        return None
    with h5py.File(f, "r") as h:
        if any(a not in h["/Header"].attrs for a in attrs):
            return None
        return [h["/Header"].attrs[a] for a in attrs]


def test_soundwave(root, schemes, out, test="soundwave"):
    """Linear sound wave in 3D: travelling (phase speed, damping, L1 of the field, transverse noise)
    or standing (frequency and damping of the mode energy over many periods)."""
    metrics = {}
    series = {}
    final = {}
    params = {}
    for s in schemes:
        run_dir = os.path.join(root, "runs", s, test)
        att = ic_attrs(run_dir, "soundWave.hdf5", "SoundWaveAmplitude", "SoundWaveNWave",
                       "SoundWaveStanding", "SoundWaveBoost", "SoundWaveSoundSpeed", "BoxSize")
        if att is None:
            continue
        A, n, standing, boost, cs, L = (float(att[0]), int(att[1]), bool(att[2]), float(att[3]),
                                        float(att[4]), float(np.atleast_1d(att[5])[0]))
        k = 2.0 * np.pi * n / L
        omega = cs * k
        params[s] = (A, k, omega, standing, boost, cs, L)

        def project(d):
            x = np.mod(d["pos"][:, 0] - boost * d["time"], L)
            w = d["m"] / np.sum(d["m"])
            e = np.exp(-1j * k * x)
            Zv = 2.0 * np.sum(w * (d["vel"][:, 0] - boost) * e)
            Zr = 2.0 * np.sum(w * (d["rho"] / np.sum(w * d["rho"]) - 1.0) * e)
            return [Zv.real, Zv.imag, Zr.real, Zr.imag]

        tt, zz = snapshot_series(run_dir, "soundWave", project)
        if len(tt) < 4:
            continue
        Zv = zz[:, 0] + 1j * zz[:, 1]
        Zr = zz[:, 2] + 1j * zz[:, 3]
        series[s] = (tt, Zv, Zr)
        d = load(snapshot_files(run_dir, "soundWave")[-1])
        final[s] = d
        x = np.mod(d["pos"][:, 0] - boost * d["time"], L)
        if standing:
            v_ex = A * cs * np.sin(k * x) * np.sin(omega * d["time"])
        else:
            v_ex = A * cs * np.cos(k * x - omega * d["time"])
        if standing:
            r_ex = A * np.cos(k * x) * np.cos(omega * d["time"])
        else:
            r_ex = A * np.cos(k * x - omega * d["time"])
        m = metrics.setdefault(s, {})
        m["L1_v_end"] = float(np.mean(np.abs(d["vel"][:, 0] - boost - v_ex)) / (A * cs))
        # Particle scatter about the exact wave, in units of the amplitude:
        # the noise a scheme adds to a linear wave (1 = as large as the wave)
        m["scatter_v_end"] = float(np.sqrt(np.mean((d["vel"][:, 0] - boost - v_ex) ** 2)) / (A * cs))
        m["scatter_rho_end"] = float(np.sqrt(np.mean((d["rho"] / np.mean(d["rho"]) - 1.0 - r_ex) ** 2)) / A)
        m["transverse_noise_end"] = float(np.sqrt(np.mean(d["vel"][:, 1] ** 2 + d["vel"][:, 2] ** 2)) / (A * cs))
        if standing:
            E = np.abs(Zv) ** 2 + cs**2 * np.abs(Zr) ** 2
            gam = -0.5 * np.polyfit(tt, np.log(E / (A * cs) ** 2), 1)[0]
            a = Zr.real / A
            crossings = np.where(np.diff(np.sign(a)) != 0)[0]
            if len(crossings) > 2:
                # linear interpolation of the zero crossings
                tc = tt[crossings] - a[crossings] * (tt[crossings + 1] - tt[crossings]) / (a[crossings + 1] - a[crossings])
                omega_meas = np.pi / np.mean(np.diff(tc))
            else:
                omega_meas = float("nan")
            m["omega_ratio"] = float(omega_meas / omega)
            m["damping_per_period"] = float(gam * 2.0 * np.pi / omega)
            m["E_mode_end"] = float(E[-1] / (A * cs) ** 2)
        else:
            phase = np.unwrap(np.angle(Zv))
            omega_meas = -np.polyfit(tt, phase, 1)[0]
            gam = -np.polyfit(tt, np.log(np.abs(Zv) / (A * cs)), 1)[0]
            m["phase_speed_ratio"] = float(omega_meas / omega)
            m["damping_per_period"] = float(gam * 2.0 * np.pi / omega)
            m["amplitude_end"] = float(np.abs(Zv[-1]) / (A * cs))
    if not series:
        return {}
    present = [s for s in schemes if s in series]
    A, k, omega, standing, boost, cs, L = params[present[0]]
    fig, axes = plt.subplots(2, len(present) + 1, figsize=(3.6 * (len(present) + 1), 3.0 * 2), squeeze=False)
    xs = np.linspace(0, L, 400)
    for j, s in enumerate(present):
        d = final[s]
        t_end = d["time"]
        x = np.mod(d["pos"][:, 0] - boost * d["time"], L)
        if standing:
            v_ex = np.sin(k * xs) * np.sin(omega * t_end)
            r_ex = np.cos(k * xs) * np.cos(omega * t_end)
        else:
            v_ex = np.cos(k * xs - omega * t_end)
            r_ex = v_ex
        for row, (y, ex, lab) in enumerate([((d["vel"][:, 0] - boost) / (A * cs), v_ex, r"$v_x / (A c_s)$"),
                                            ((d["rho"] / np.mean(d["rho"]) - 1.0) / A, r_ex, r"$\delta\rho / (A \rho_0)$")]):
            ax = axes[row, j]
            ax.scatter(x, y, color="0.4", **scatter_props)
            xb, yb, sb = binned(x, y, np.linspace(0, L, 33))
            ax.errorbar(xb, yb, yerr=sb, fmt=".", color=scheme_color(s, j), ms=4, lw=1, zorder=3)
            ax.plot(xs, ex, "k-", lw=1)
            ax.set_ylim(-2.5, 2.5)
            ax.set_xlabel("x")
            ax.set_ylabel(lab)
            sc = metrics[s]["scatter_v_end" if row == 0 else "scatter_rho_end"]
            ax.text(0.03, 0.95, f"scatter {sc:.2f}", transform=ax.transAxes, fontsize=8, va="top",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))
            if row == 0:
                ax.set_title(s, fontsize=9)
        tt, Zv, Zr = series[s]
        c = scheme_color(s, j)
        if standing:
            E = np.abs(Zv) ** 2 + cs**2 * np.abs(Zr) ** 2
            axes[0, -1].semilogy(tt, E / (A * cs) ** 2, "-", color=c, lw=1.2,
                                 label=f"{s}: {metrics[s]['damping_per_period']:.1e} / period")
            axes[1, -1].plot(tt, Zr.real / A, "-", color=c, lw=1, label=f"{s}: omega / (c k) = {metrics[s]['omega_ratio']:.4f}")
        else:
            axes[0, -1].plot(tt, np.abs(Zv) / (A * cs), "-", color=c, lw=1.2,
                             label=f"{s}: {metrics[s]['damping_per_period']:.1e} / period")
            axes[1, -1].plot(tt, np.unwrap(np.angle(Zv)) + omega * tt, "-", color=c, lw=1.2,
                             label=f"{s}: c / c_s = {metrics[s]['phase_speed_ratio']:.4f}")
    if standing:
        axes[1, -1].plot(xs * 0 + tt[0], xs * 0, "k:", lw=0.5)
        axes[0, -1].set_ylabel(r"$E_{\rm mode} / (A c_s)^2$")
        axes[1, -1].set_ylabel(r"Re $Z_\rho / A$")
        kind = "standing"
    else:
        axes[0, -1].set_ylabel(r"$|Z_v| / (A c_s)$")
        axes[1, -1].set_ylabel("phase $-$ exact [rad]")
        kind = "travelling"
    for ax in axes[:, -1]:
        ax.set_xlabel("t")
        ax.legend(fontsize=6)
    fig.suptitle(f"{kind} sound wave (3D, A = {A:g}, {2 * np.pi * n / L / k * round(len(final[present[0]]['m']) ** (1 / 3)):.0f} particles per wavelength), t = {t_end:.2f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"{test}.png"), dpi=150)
    plt.close(fig)
    return metrics


def test_thermal(root, schemes, out):
    """Thermal discontinuity in pressure equilibrium: width of the u profiles (effective conductivity),
    pressure blip and spurious velocities."""
    from scipy.optimize import curve_fit
    from scipy.special import erf

    metrics = {}
    series = {}
    profiles = {}
    for s in schemes:
        run_dir = os.path.join(root, "runs", s, "thermal")
        att = ic_attrs(run_dir, "thermalDiffusion.hdf5", "ThermalDiffusionContrast", "ThermalDiffusionPressure",
                       "ThermalDiffusionInterfaces", "BoxSize")
        if att is None:
            continue
        C, P0 = float(att[0]), float(att[1])
        x1, x2 = [float(v) for v in att[2]]
        Lx = float(np.atleast_1d(att[3])[0])
        u_hot = P0 / (GAS_GAMMA - 1.0)
        u_cold = u_hot / C
        cs_cold = np.sqrt(GAS_GAMMA * P0 / C)
        edges = np.linspace(0, Lx, 201)

        def width(x, u, xi, sign):
            sel = np.abs(x - xi) < 0.3
            theta = (u[sel] - u_cold) / (u_hot - u_cold)
            model = lambda xx, w, x0: 0.5 * (1.0 + sign * erf((xx - x0) / (np.sqrt(2.0) * w)))
            try:
                p, _ = curve_fit(model, x[sel], theta, p0=[0.02, xi], bounds=([1e-4, xi - 0.1], [0.3, xi + 0.1]))
                return p[0]
            except RuntimeError:
                return float("nan")

        def measure(d):
            x = d["pos"][:, 0]
            xb, ub, _ = binned(x, d["u"], edges)
            _, Pb, _ = binned(x, d["P"], edges)
            vmag = np.sqrt(np.sum(d["vel"] ** 2, axis=1))
            return [0.5 * (width(x, d["u"], x1, -1) + width(x, d["u"], x2, +1)),
                    np.nanmax(np.abs(Pb / P0 - 1.0)), np.max(vmag) / cs_cold,
                    np.sqrt(np.mean(vmag**2)) / cs_cold] + list(ub / u_hot) + list(Pb / P0)

        tt, vals = snapshot_series(run_dir, "thermalDiffusion", measure)
        if len(tt) < 2:
            continue
        sig, blip, vmax, vrms = vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3]
        nb = len(edges) - 1
        profiles[s] = (0.5 * (edges[1:] + edges[:-1]), vals[:, 4:4 + nb], vals[:, 4 + nb:4 + 2 * nb])
        series[s] = (tt, sig, blip, vmax, vrms)
        ok = np.isfinite(sig)
        kappa = np.polyfit(tt[ok], sig[ok] ** 2, 1)[0] / 2.0 if np.sum(ok) > 2 else float("nan")
        m = metrics.setdefault(s, {})
        m["sigma_0"] = float(sig[0])
        m["sigma_end"] = float(sig[-1])
        m["kappa_eff"] = float(kappa)
        m["P_blip_max"] = float(np.max(blip))
        m["P_blip_end"] = float(blip[-1])
        m["v_max_max"] = float(np.max(vmax))
        m["v_rms_end"] = float(vrms[-1])
    if not series:
        return {}
    present = [s for s in schemes if s in series]
    fig, axes = plt.subplots(2, len(present) + 1, figsize=(3.6 * (len(present) + 1), 3.0 * 2), squeeze=False)
    for j, s in enumerate(present):
        xb, up, Pp = profiles[s]
        tt = series[s][0]
        for i in range(len(tt)):
            c = plt.cm.viridis(i / max(len(tt) - 1, 1))
            axes[0, j].plot(xb, up[i], color=c, lw=0.8, label=f"t = {tt[i]:.2f}" if i in (0, len(tt) - 1) else None)
            axes[1, j].plot(xb, Pp[i], color=c, lw=0.8)
        axes[0, j].set_title(s, fontsize=9)
        axes[0, j].set_ylabel(r"$u / u_{\rm hot}$")
        axes[0, j].legend(fontsize=6)
        axes[1, j].set_ylabel(r"$P / P_0$")
        axes[1, j].set_ylim(0.8, 1.2)
        for ax in axes[:, j]:
            ax.set_xlabel("x")
        tt, sig, blip, vmax, vrms = series[s]
        c = scheme_color(s, j)
        axes[0, -1].plot(tt, sig, "o-", color=c, ms=3, lw=1, label=f"{s}: kappa = {metrics[s]['kappa_eff']:.1e}")
        axes[1, -1].semilogy(tt, vrms, "-", color=c, lw=1.2, label=f"{s} rms")
        axes[1, -1].semilogy(tt, blip, "--", color=c, lw=1, label=f"{s} P blip")
    axes[0, -1].set_ylabel(r"interface width $\sigma$")
    axes[1, -1].set_ylabel(r"$v_{\rm rms} / c_{s,\rm cold}$, max $|P/P_0 - 1|$")
    for ax in axes[:, -1]:
        ax.set_xlabel("t")
        ax.legend(fontsize=6)
    fig.suptitle(f"Thermal discontinuity in pressure equilibrium (2D, contrast {C:g}), t = {tt[-1]:.2f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "thermal.png"), dpi=150)
    plt.close(fig)
    return metrics


TESTS = {"sod": test_sod, "sedov": test_sedov, "noh": test_noh, "gresho": test_gresho,
         "evrard": test_evrard, "kh": test_kh, "square": test_square,
         "square_mass": lambda root, schemes, out: test_square(root, schemes, out, test="square_mass"),
         "keplerian": test_keplerian,
         "keplerian2d": lambda root, schemes, out: test_keplerian(root, schemes, out, test="keplerian2d"),
         "keplerian_cd10": lambda root, schemes, out: test_keplerian(root, schemes, out, test="keplerian_cd10"),
         "zeldovich": test_zeldovich,
         "zeldovich_glass": lambda root, schemes, out: test_zeldovich(root, schemes, out, test="zeldovich_glass"),
         "zeldovich_pert": lambda root, schemes, out: test_zeldovich(root, schemes, out, test="zeldovich_pert"),
         "blob": test_blob, "nfw": test_nfw,
         "taylorgreen": test_taylorgreen, "shearing": test_shearing, "blast": test_blast,
         "rt": test_rt, "khgr": test_khgr, "ibw": test_ibw,
         "soundwave": test_soundwave,
         "wavedamping": lambda root, schemes, out: test_soundwave(root, schemes, out, test="wavedamping"),
         "thermal": test_thermal,
         "gresho_boost": lambda root, schemes, out: test_gresho(root, schemes, out, test="gresho_boost"),
         "kh8": lambda root, schemes, out: test_kh(root, schemes, out, test="kh8"),
         "sod_weak": lambda root, schemes, out: test_sod(root, schemes, out, test="sod_weak"),
         "sod_strong": lambda root, schemes, out: test_sod(root, schemes, out, test="sod_strong")}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    # plot_comparison.py ROOT [--schemes=a,b,c] [--out=DIR] [--overlay-only=d] [test ...]
    # --schemes restricts the comparison to a subset of runs/<scheme> (for
    # large parameter studies), --out puts the figures and summary elsewhere
    # than ROOT/plots, --overlay-only lists schemes that get no column of
    # their own in the image-based figures (KH, squares, blob), only their
    # curve in the panel that overlays all schemes.
    root = sys.argv[1]
    args = sys.argv[2:]
    only = None
    out = os.path.join(root, "plots")
    tests = []
    for a in args:
        if a.startswith("--schemes="):
            only = a[len("--schemes="):].split(",")
        elif a.startswith("--out="):
            out = a[len("--out="):]
        elif a.startswith("--overlay-only="):
            OVERLAY_ONLY.update(a[len("--overlay-only="):].split(","))
        else:
            tests.append(a)
    tests = tests or ALL_TESTS
    schemes = sorted(os.listdir(os.path.join(root, "runs")))
    if only is not None:
        schemes = [s for s in schemes if s in only]
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
