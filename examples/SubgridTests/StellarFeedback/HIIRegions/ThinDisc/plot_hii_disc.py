"""
Plot the gas disc face-on (x-y) and edge-on (x-z) with the HII-tagged gas
coloured and the star particles marked.

Usage:
    python3 plot_hii_disc.py                       # last snapshot in snap/
    python3 plot_hii_disc.py -s snap/snapshot -i 4 --half-width-pc 500
"""

import argparse
import glob

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


def parse_options():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-s", dest="base", default="snap/snapshot", help="Snapshot basename")
    p.add_argument(
        "-i",
        dest="index",
        type=int,
        default=None,
        help="Snapshot index (default: last)",
    )
    p.add_argument(
        "--half-width-pc", type=float, default=500.0, help="Face-on half width [pc]"
    )
    p.add_argument(
        "--half-height-pc", type=float, default=150.0, help="Edge-on half height [pc]"
    )
    p.add_argument("-o", dest="out", default=None, help="Output PNG")
    return p.parse_args()


def column_map(x, y, w, xr, yr, nx, ny):
    img, _, _ = np.histogram2d(x, y, bins=[nx, ny], range=[xr, yr], weights=w)
    return img.T


def main():
    opt = parse_options()
    if opt.index is None:
        files = sorted(glob.glob(f"{opt.base}_*.hdf5"))
        if not files:
            raise SystemExit(f"No snapshot found for {opt.base}")
        fname = files[-1]
    else:
        fname = f"{opt.base}_{opt.index:04d}.hdf5"

    with h5py.File(fname, "r") as f:
        units = f["Units"].attrs
        unit_l_pc = (
            float(np.ravel(units["Unit length in cgs (U_L)"])[0])
            / 3.0856775814913673e18
        )
        unit_m_msun = float(np.ravel(units["Unit mass in cgs (U_M)"])[0]) / 1.98841e33
        unit_t_myr = float(np.ravel(units["Unit time in cgs (U_t)"])[0]) / 3.15576e13
        t_myr = float(np.ravel(f["Header"].attrs["Time"])[0]) * unit_t_myr
        box = f["Header"].attrs["BoxSize"]
        centre = 0.5 * np.asarray(box)
        pos = (f["PartType0/Coordinates"][:] - centre) * unit_l_pc
        mass = f["PartType0/Masses"][:] * unit_m_msun
        tag = f["PartType0/IsIonizedFlags"][:] > 0
        star_pos = (f["PartType4/Coordinates"][:] - centre) * unit_l_pc

    hw, hh = opt.half_width_pc, opt.half_height_pc
    pix_pc = 2.0
    nx = int(2 * hw / pix_pc)
    nz = int(2 * hh / pix_pc)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(7.5, 7.5 * (1 + hh / hw) / 1.0),
        gridspec_kw={"height_ratios": [2 * hw, 2 * hh]},
        constrained_layout=True,
    )
    views = [
        (0, 1, hw, "x [pc]", "y [pc]", "Face-on"),
        (0, 2, hh, "x [pc]", "z [pc]", "Edge-on"),
    ]
    for ax, (ia, ib, hb, xl, yl, title) in zip(axes, views):
        xr, yr = [-hw, hw], [-hb, hb]
        ny = nx if ib == 1 else nz
        gas = column_map(pos[:, ia], pos[:, ib], mass, xr, yr, nx, ny) / pix_pc**2
        ion = (
            column_map(pos[tag, ia], pos[tag, ib], mass[tag], xr, yr, nx, ny)
            / pix_pc**2
        )
        ax.imshow(
            np.where(gas > 0, gas, np.nan),
            origin="lower",
            extent=xr + yr,
            cmap="Greys",
            norm=LogNorm(vmin=1e-1, vmax=1e2),
            aspect="equal",
        )
        im = ax.imshow(
            np.where(ion > 0, ion, np.nan),
            origin="lower",
            extent=xr + yr,
            cmap="plasma",
            norm=LogNorm(vmin=1e-1, vmax=1e2),
            aspect="equal",
            alpha=0.9,
        )
        ax.scatter(
            star_pos[:, ia],
            star_pos[:, ib],
            marker="*",
            s=40,
            c="cyan",
            edgecolors="k",
            linewidths=0.4,
        )
        ax.set_xlim(xr)
        ax.set_ylim(yr)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(f"{title}, t = {t_myr:.2f} Myr, {tag.sum()} tagged particles")
    cb = fig.colorbar(im, ax=axes, shrink=0.6)
    cb.set_label(r"Tagged gas surface density [M$_\odot$ pc$^{-2}$] (grey: all gas)")
    out = opt.out or fname.replace(".hdf5", "_hii_disc.png")
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
