"""Field around one star: the optically thin profile and the shadow of a clump.

Usage: python3 isrf_fig_thin.py <runs directory with shadow0/ and shadow1/>
"""

import json
import sys

import h5py
import numpy as np
import isrf_style as st
import matplotlib.pyplot as plt

PC = 3.0857e18
C = 2.99792458e10
M_H = 1.6726e-24
MU_H = 1.4
Z_SUN = 0.01295
SIGMA_PE = 9e-22


def load(path: str) -> dict:
    """Read the PE field, density and positions of one snapshot, in cgs."""
    with h5py.File(path) as f:
        u = f["Units"].attrs
        L, M, T = (
            float(u[k][0])
            for k in (
                "Unit length in cgs (U_L)",
                "Unit mass in cgs (U_M)",
                "Unit time in cgs (U_t)",
            )
        )
        g, s = f["PartType0"], f["PartType4"]
        return dict(
            t=float(np.ravel(f["Header"].attrs["Time"])[0]) * T / 3.156e13,
            x=g["Coordinates"][:].astype(np.float64) * L,
            rho=g["Densities"][:].astype(np.float64) * M / L**3,
            u=g["PacketKernelMeanSpecificEnergies"][:, 0].astype(np.float64)
            * (L / T) ** 2,
            Z=g["SmoothedMetalMassFractions"][:, -1],
            star=s["Coordinates"][0].astype(np.float64) * L,
            L_pe=float(s["PELuminosities"][0]) * M * L**2 / T**3,
        )


def kappa(rho: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """Return the PE dust absorption rate of the chapter, 1/cm, at Grackle's default dust ratio."""
    return SIGMA_PE * (Z / Z_SUN) * rho / (MU_H * M_H)


def binned(x: np.ndarray, y: np.ndarray, edges: np.ndarray):
    """Return bin centres, medians and the 16 and 84 percentiles of y in bins of x."""
    idx = np.digitize(x, edges)
    out = []
    for k in range(1, len(edges)):
        sel = y[idx == k]
        if sel.size > 5:
            out.append(
                (0.5 * (edges[k - 1] + edges[k]), *np.percentile(sel, [50, 16, 84]))
            )
    return np.array(out).T


def panel_thin(ax, runs: str) -> None:
    """Field of one star in uniform gas against L e^{-kappa r} / (4 pi r^2)."""
    for snap, col, mk in [("0001", st.SKY, "s"), ("0002", st.BLUE, "o")]:
        d = load(f"{runs}/shadow0/snap/snapshot_{snap}.hdf5")
        r = np.linalg.norm(d["x"] - d["star"], axis=1)
        flux = C * d["rho"] * d["u"]
        c, med, lo, hi = binned(r / PC, flux, np.geomspace(0.15, 8.0, 26))
        ax.fill_between(c, lo, hi, color=col, alpha=0.25, lw=0)
        ax.plot(
            c, med, mk, color=col, ms=3, label=f"run, $t={d['t']:.2f}$ Myr (median)"
        )
    k = np.median(kappa(d["rho"], d["Z"]))
    rr = np.geomspace(0.15, 8.0, 200)
    ax.plot(
        rr,
        d["L_pe"] * np.exp(-k * rr * PC) / (4 * np.pi * (rr * PC) ** 2),
        color=st.BLACK,
        lw=1.0,
        label=r"$L_{\rm PE}\,e^{-\kappa r}/(4\pi r^2)$",
    )
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="distance to the star $r$ [pc]",
        ylabel=r"$c\,\rho\,\tilde u_{\rm PE}$ [erg s$^{-1}$ cm$^{-2}$]",
    )
    ax.set_ylim(1e-3, 3e1)
    ax.legend(frameon=False, loc="lower left")
    ax.set_title("(a) one star in uniform gas, PE band")


def panel_shadow(ax, runs: str) -> None:
    """Field with a clump over the field without it, along the axis through the clump."""
    geo = json.load(open(f"{runs}/shadow1/geometry.json"))
    star = np.array(geo["star_kpc"]) * 1e3 * PC
    axis = np.array(geo["axis"])
    centre = (np.array(geo["clump_centre_kpc"]) * 1e3 * PC - star) @ axis
    radius = geo["clump_radius_kpc"] * 1e3 * PC
    edges = np.arange(0.0, 4.6, 0.25) * PC
    med = {}
    for run in ("shadow0", "shadow1"):
        d = load(f"{runs}/{run}/snap/snapshot_0002.hdf5")
        x = d["x"] - star
        s = x @ axis
        perp = np.linalg.norm(x - np.outer(s, axis), axis=1)
        sel = perp < 0.4 * radius
        c, m, _, _ = binned(s[sel], (d["rho"] * d["u"])[sel], edges)
        med[run] = m
        if run == "shadow1":
            inside = np.linalg.norm(x - centre * axis, axis=1) < radius
            k_c = np.median(kappa(d["rho"][inside], d["Z"][inside]))
    ratio = med["shadow1"] / med["shadow0"]
    ax.axvspan((centre - radius) / PC, (centre + radius) / PC, color=st.LIGHT, lw=0)
    ax.text(centre / PC, 0.9, "clump", ha="center", fontsize=8)
    ax.plot(c / PC, ratio, "o", color=st.BLUE, ms=3, label="run, with / without clump")
    ss = np.linspace(centre - radius, centre + radius, 100)
    ax.plot(
        ss / PC,
        np.exp(-k_c * (ss - (centre - radius))),
        color=st.BLACK,
        lw=1.0,
        label=r"$e^{-\kappa_{\rm c}\,\ell}$, $\ell$: depth",
    )
    ax.set(yscale="log", xlabel="distance along the axis [pc]", ylabel="field ratio")
    ax.set_ylim(1e-4, 4)
    ax.set_xlim(0, 4.6)
    ax.legend(frameon=False, loc="lower left")
    ax.set_title("(b) shadow behind a dusty clump")


def main() -> None:
    """Draw the figure."""
    runs = (
        sys.argv[1]
        if len(sys.argv) > 1
        else "/run/media/darwinr/SimulationData/claude/ISRF/report_figs/runs"
    )
    st.setup()
    fig, (a, b) = plt.subplots(1, 2, figsize=(st.TEXT_WIDTH, 2.7))
    panel_thin(a, runs)
    panel_shadow(b, runs)
    fig.tight_layout(w_pad=1.2)
    st.save(fig, "isrf_thin.pdf")


if __name__ == "__main__":
    main()
