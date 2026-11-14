###############################################################################
# This file is part of SWIFT.
# Copyright (c) ...
#               2026 Nikyta Shchutskyi (shchutskyi@lorentz.leidenuniv.nl)
#
# This file plots the magnetic-field patterns of the Roberts-flow dynamo eigensolutions
# and compares them with SWIFT / Pencil slices.
# (Roberts 1972, Phil. Trans. R. Soc. A 271, 411, sections 8 and 10)
# (Shchutskyi, Brandenburg & Schaller, Convergence of Roberts flow dynamos in Eulerian
#  and Lagrangian codes, 2026)
#
# Initial version developed with the assistance of Claude (Opus 5.5, Anthropic),
# November 2026; reviewed and validated by the authors above.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################
"""Magnetic-field patterns (r_i, theta_i; eqs. 2.7-2.8 of the paper) of the Roberts-flow eigensolution,
in the same convention and layout as FT_slice_pixel_show.py, so they can be compared pixel by pixel
with SWIFT / Pencil slices.

Convention (identical to FT_slice_pixel_show.py):
  * B is normalised by its box rms, B0 = B / Brms                                    (2.6)
  * b(x, y) = (1/Lz) int B0 e^{-i kz z} dz, i.e. the coefficient of e^{+i kz z}  (numpy fft / N)
    For the eigenmode B = Re[bt(x, y) e^{i kz z}] this is b = bt / (2 Brms).
  * one global phase factor makes Im<b_x>_xy = 0 and Re<b_x>_xy > 0                  (2.8)
  * r_i = |b_i|, theta_i = arg b_i on pixels x_i = i Lx/N, i = 0..N-1 (pixel centres of the
    SWIFT slices, which compensate the half-cell shift), then periodically continued.

The growing solution can consist of two modes with the same growth rate: a +-omega pair (flows II,
III) or a degenerate real eigenvalue (flow IV).  The paper excites only the mode with non-zero
plane-averaged B_x (initial condition 2.12); `polarization='x'` selects that mode (use 'y' for
the other one).

Usage:
    python eigen_patterns.py -flow I -eta 0.05 -kz 1.0 -res_xy 64 -outim eig_I.png -out eig_I.h5
    python eigen_patterns.py -flow I -eta 0.05 -kz 1.0 -res_xy 64 -compare slice.h5   # vs code slice

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import scipy.linalg as la, scipy.sparse.linalg as sla
from solver import converge, build, DENSE_MAX

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm"})


# ----------------------------------------------------------------------------- eigenfunction
def _full_coeffs(vec, s, kz):
    """Eigenvector of solver.build (Hx, Hy coefficients, |m|,|n| <= s) -> array (3, N, N) of the
    Fourier coefficients of all three components of bt(x, y), with Hz = -(m Hx + n Hy)/kz from
    div B = 0.  Index [component, m + s, n + s]."""
    N = 2 * s + 1
    m = np.arange(-s, s + 1)
    M, Nn = np.meshgrid(m, m, indexing="ij")
    hx, hy = vec.reshape(2, N, N)
    return np.array([hx, hy, -(M * hx + Nn * hy) / kz])


def eigenmode(name, eta, kz, w0=1.0, tol=1e-8, polarization="x", degen_tol=1e-8):
    """Converged fastest eigenmode at a fixed kz.

    Parameters
        name, eta, kz, w0 : flow, diffusivity, vertical wavenumber, vertical flow amplitude
        tol               : absolute tolerance on p passed to solver.converge
        polarization      : 'x' or 'y', the horizontal direction of the plane-averaged field
                            that the simulation's initial condition seeds (eq. 2.12 -> 'x')
        degen_tol         : relative distance below which two eigenvalues count as equal
    The leading growth rate can belong to more than one mode:
      * a +-omega pair (flows II, III): the member with the larger plane-averaged
        |B_polarization| is returned;
      * a degenerate real eigenvalue (flow IV): the combination with zero plane-averaged
        B in the other horizontal direction is returned, i.e. the mode an initial condition
        with only B_polarization (eq. 2.12) excites.
    Returns
        dict(p = growth rate + i*frequency, s, coeffs (3, N, N) from _full_coeffs, kz,
             converged, note = which selection was made, '' if the mode is unique)
    """
    r = converge(name, eta, kz, w0, tol=tol)
    s, pstar = r["s"], r["p"] + eta * kz**2
    N = 2 * s + 1
    A = build(name, eta, kz, s, w0)
    targets = [pstar] + ([np.conj(pstar)] if abs(pstar.imag) > 1e-8 else [])
    if A.shape[0] <= DENSE_MAX:
        ev, V = la.eig(A.toarray())
    else:
        ev, V = zip(*[sla.eigs(A, k=6, sigma=t) for t in targets])
        ev, V = np.concatenate(ev), np.concatenate(V, axis=1)
    i_pol, i_oth = {"x": (0, 1), "y": (1, 0)}[polarization]
    mean = lambda v: v.reshape(2, N, N)[:, s, s]  # plane-averaged (Hx, Hy)
    modes, note = [], ""
    for t in targets:
        idx = [
            i
            for i in np.argsort(abs(ev - t))
            if abs(ev[i] - t) < degen_tol * max(1, abs(t))
        ]
        if not idx:
            continue
        v = V[:, idx[0]]
        if len(idx) >= 2:  # degenerate: remove other component
            v1, v2 = V[:, idx[0]], V[:, idx[1]]
            v = mean(v2)[i_oth] * v1 - mean(v1)[i_oth] * v2
            v /= np.linalg.norm(v)
            note = f'degenerate eigenvalue (multiplicity {len(idx)}); <B_{"xy"[i_oth]}> = 0 member chosen'
        modes.append((ev[idx[0]] - eta * kz**2, _full_coeffs(v, s, kz)))
    if len(targets) == 2 and len(modes) == 2:
        note = note or "+-omega pair; member with larger <B_%s> chosen" % polarization
    p, c = max(modes, key=lambda pc: abs(pc[1][i_pol, s, s]))
    return dict(p=p, s=s, coeffs=c, kz=kz, converged=r["converged"], note=note)


def pattern(mode, res_xy=64, Lxy=2 * np.pi, x0=0.0):
    """Complex map of the mode on a pixel grid, in the FT_slice_pixel_show.py convention.

    Parameters
        mode   : dict from eigenmode()
        res_xy : pixels per direction (use the resolution of the code slice)
        Lxy    : horizontal box size
        x0     : coordinate of the first pixel centre (pixels at x_i = x0 + i Lxy/res_xy)
    Returns
        b : complex array (res_xy, res_xy, 3) indexed [x, y, component], with
            b = bt/(2 Brms) (the e^{+i kz z} coefficient of B/Brms) and the global phase fixed
            so that <b_x>_xy is real and positive.  r_i = |b_i|, theta_i = arg b_i.
    Evaluated exactly as a direct Fourier sum, so there is no aliasing at low res_xy.
    """
    c = mode["coeffs"]
    s = (c.shape[1] - 1) // 2
    m = np.arange(-s, s + 1)
    x = x0 + Lxy * np.arange(res_xy) / res_xy
    E = np.exp(1j * np.outer(x, m))  # (res_xy, N)
    bt = np.einsum("xm,cmn,yn->xyc", E, c, E)  # sum_mn c_mn e^{i(m x + n y)}
    Brms = np.sqrt(
        0.5 * np.sum(np.abs(c) ** 2)
    )  # <|Re(bt e^{ikz z})|^2>^(1/2), Parseval
    b = bt / (2 * Brms)
    bx_mean = c[0, s, s] / (2 * Brms)  # exact plane average of b_x
    if abs(bx_mean) > 1e-12 * np.abs(b).max():
        b = b / (bx_mean / abs(bx_mean))  # Im<b_x> = 0, Re<b_x> > 0
    return b


# ----------------------------------------------------------------------------- plotting (as FT_slice_pixel_show.py)
def periodic_continue(Q):
    """Append the first row and column at the end, (M, N) -> (M+1, N+1), so that pcolormesh
    with shading='nearest' shows the full periodic cell (same as FT_slice_pixel_show.py).
    """
    M, N = Q.shape
    Qp = np.zeros((M + 1, N + 1), dtype=Q.dtype)
    Qp[:-1, :-1] = Q
    Qp[-1, :-1] = Q[0, :]
    Qp[:-1, -1] = Q[:, 0]
    Qp[-1, -1] = Q[0, 0]
    return Qp


def plot_pattern(
    b,
    out_image=None,
    title=None,
    cmapr="plasma",
    cmaptheta="twilight_shifted",
    rmax=1.5,
):
    """6-panel figure: r_1..3 (top row) and theta_1..3 (bottom row), same style, colour maps and
    ticks as FT_slice_pixel_show.py.

    Parameters
        b          : complex map from pattern() or load_slice(), shape (res, res, 3)
        out_image  : file name to save (None -> not saved)
        title      : optional figure title
        cmapr, cmaptheta, rmax : colour maps and upper limit of the r colour scale
    Returns
        the matplotlib Figure
    """
    res_xy = b.shape[0]
    r, th = np.abs(b), np.angle(b)
    nx, ny = 3, 2
    fig, ax = plt.subplots(ny, nx, sharey=True, sharex=True, figsize=(5 * nx, 5 * ny))
    for a in ax.flat:
        a.tick_params(
            axis="both",
            which="major",
            direction="in",
            top=True,
            right=True,
            length=7,
            width=1.2,
            labelsize=16,
        )
        a.tick_params(
            axis="both",
            which="minor",
            direction="in",
            top=True,
            right=True,
            length=4,
            width=1.0,
        )
        a.xaxis.set_minor_locator(AutoMinorLocator())
        a.yaxis.set_minor_locator(AutoMinorLocator())
    ticks_r = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
    ticks_phase = [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi]
    ticklabels_phase = [
        r"$-\pi$",
        r"$-\frac{\pi}{2}$",
        r"$0$",
        r"$\frac{\pi}{2}$",
        r"$\pi$",
    ]
    for i in range(3):
        im = ax[0, i].pcolormesh(
            periodic_continue(r[:, :, i].T),
            cmap=cmapr,
            vmin=0,
            vmax=rmax,
            shading="nearest",
        )
        cb = plt.colorbar(im, ax=ax[0, i], fraction=0.046, pad=0.04, ticks=ticks_r)
        cb.set_ticklabels(ticks_r, fontsize=15)
        ax[0, i].set_title(rf"$r_{i+1} \;$", fontsize=30, pad=10)
        ax[0, i].tick_params(color="w", which="both")
        im = ax[1, i].pcolormesh(
            periodic_continue(th[:, :, i].T),
            cmap=cmaptheta,
            vmin=-np.pi,
            vmax=np.pi,
            shading="nearest",
        )
        cb = plt.colorbar(im, ax=ax[1, i], fraction=0.046, pad=0.04, ticks=ticks_phase)
        cb.set_ticklabels(ticklabels_phase, fontsize=20)
        ax[1, i].set_title(rf"$\theta_{i+1} \;$", fontsize=30, pad=10)
    locs = [res_xy / 4, res_xy / 2, 3 * res_xy / 4, res_xy]
    labels = [r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]
    for a in ax.flat:
        a.set_xticks(locs, labels, fontsize=25)
        a.set_yticks(locs, labels, fontsize=25)
        a.set_aspect("equal", adjustable="box")
        a.set_xlim(0, res_xy)
        a.set_ylim(0, res_xy)
    if title:
        fig.suptitle(title, fontsize=22)
    fig.tight_layout()
    if out_image:
        fig.savefig(out_image, dpi=220)
    return fig


def save_pattern(b, filename):
    """Write the map to an HDF5 file with the same datasets as FT_slice_pixel_show.py:
    r1..r3, theta1..theta3 (periodically continued and transposed, as plotted) and
    ReB1..ReB3, ImB1..ImB3 (raw (res_xy, res_xy) arrays indexed [x, y]).  Requires h5py.
    """
    import h5py

    r, th = np.abs(b), np.angle(b)
    with h5py.File(filename, "w") as f:
        for i in range(3):
            f.create_dataset(f"r{i+1}", data=periodic_continue(r[:, :, i].T))
            f.create_dataset(f"theta{i+1}", data=periodic_continue(th[:, :, i].T))
            f.create_dataset(f"ReB{i+1}", data=np.real(b[:, :, i]))
            f.create_dataset(f"ImB{i+1}", data=np.imag(b[:, :, i]))


# ----------------------------------------------------------------------------- comparison with a code slice
def load_slice(filename):
    """Read a slice .h5 written by FT_slice_pixel_show.py or save_pattern() and return the
    complex map b = ReB + i ImB, shape (res, res, 3) indexed [x, y, component]."""
    import h5py

    with h5py.File(filename, "r") as f:
        return np.stack(
            [f[f"ReB{i}"][()] + 1j * f[f"ImB{i}"][()] for i in (1, 2, 3)], axis=-1
        )


def compare(b_ref, b_code, out_image=None, labels=("eigensolution", "code")):
    """Compare a code slice with the eigensolution.

    Parameters
        b_ref, b_code : complex maps (res, res, 3) on the same pixels, same normalisation and
                        phase convention (eigensolution from pattern(), code from load_slice())
        out_image     : if given, a 3-panel figure |B~| (eq. 2.30) for both and |dB~| (2.31)
        labels        : panel labels for the two maps
    Returns
        <|dB~|>_rms over the plane, eq. (2.32)
    """
    d = np.sqrt(np.sum(np.abs(b_code - b_ref) ** 2, axis=-1))
    rms = np.sqrt(np.mean(d**2))
    if out_image:
        fig, ax = plt.subplots(1, 3, figsize=(16, 5))
        for a, q, t in zip(
            ax,
            [
                np.sqrt(np.sum(np.abs(b_ref) ** 2, -1)),
                np.sqrt(np.sum(np.abs(b_code) ** 2, -1)),
                d,
            ],
            [
                rf"$|\tilde B|$ {labels[0]}",
                rf"$|\tilde B|$ {labels[1]}",
                rf"$|\delta\tilde B|$, rms $= {rms:.2e}$",
            ],
        ):
            im = a.pcolormesh(periodic_continue(q.T), cmap="plasma", shading="nearest")
            plt.colorbar(im, ax=a, fraction=0.046, pad=0.04)
            a.set_title(t, fontsize=16)
            a.set_aspect("equal")
        fig.tight_layout()
        fig.savefig(out_image, dpi=200)
    return rms


if __name__ == "__main__":
    import argparse as ap

    pa = ap.ArgumentParser(
        description="Eigensolution magnetic-field patterns for Roberts flows"
    )
    pa.add_argument(
        "-flow", default="I", help="key of solver.FLOWS (I, II, III, IV, R72_I, ...)"
    )
    pa.add_argument("-eta", type=float, default=0.05)
    pa.add_argument(
        "-kz", type=float, default=1.0, help="vertical wavenumber of the mode"
    )
    pa.add_argument("-w0", type=float, default=1.0)
    pa.add_argument(
        "-res_xy", type=int, default=64, help="pixels in x and y (match the code slice)"
    )
    pa.add_argument(
        "-x0", type=float, default=0.0, help="coordinate of the first pixel centre"
    )
    pa.add_argument(
        "-pol", default="x", choices=["x", "y"], help="member of a +-omega pair"
    )
    pa.add_argument("-tol", type=float, default=1e-8)
    pa.add_argument("-outim", default="eig_pattern.png")
    pa.add_argument("-out", default="eig_slice.h5")
    pa.add_argument(
        "-compare", default="", help="slice .h5 from FT_slice_pixel_show.py"
    )
    a = pa.parse_args()

    # 1. eigenmode at the requested kz, 2. pixel map, 3. figure/file, 4. optional comparison
    mode = eigenmode(a.flow, a.eta, a.kz, a.w0, a.tol, a.pol)
    print(
        f"flow {a.flow}, eta={a.eta}, kz={a.kz}: lambda={mode['p'].real:.10f}, "
        f"omega={mode['p'].imag:+.10f}, s={mode['s']}, converged={mode['converged']}"
        + (f"\n  note: {mode['note']}" if mode["note"] else "")
    )
    b = pattern(mode, a.res_xy, x0=a.x0)
    plot_pattern(b, a.outim)
    if a.out:
        try:
            save_pattern(b, a.out)
        except ImportError:
            print("h5py not installed: slice file not written")
    if a.compare:
        rms = compare(b, load_slice(a.compare), a.outim.replace(".png", "_diff.png"))
        print(f"<|dB~|>_rms (eq. 2.32) = {rms:.4e}")
