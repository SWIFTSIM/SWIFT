###############################################################################
# This file is part of SWIFT.
# Copyright (c) ...
#               2026 Nikyta Shchutskyi (shchutskyi@lorentz.leidenuniv.nl)
#
# This file solves Roberts-flow kinematic dynamos I-IV (Roberts 1972) in a periodic box.
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
"""Eigenvalue solver for the kinematic Roberts-flow dynamos I-IV (Roberts 1972) in a periodic box.

Flow (k0 = 1, period 2*pi in x and y), defined in FLOWS (add new flows there):
    u = (v0 sin x cos y, -v0 cos x sin y, w0 * vz(x, y))
    vz = sin x sin y (I), cos x cos y (II), (cos 2x + cos 2y)/2 (III), sin x (IV)
    'R72_I', 'R72_II', 'R72_III', 'R72_IV': Roberts' (1972) original motions I-IV
Field  B = H(x, y) exp(i kz z + p t).  Only Hx, Hy are needed (Roberts 1972 reduction;
Hz follows from div B = 0):
    p* H = -i kz uz H - (u.grad)H + (H.grad)u + eta lap2D H,     p = p* - eta kz^2

Method: Fourier-Galerkin truncation |m|,|n| <= s (Roberts 1972, Phil. Trans. R. Soc. A 271, 411,
sections 8 and 10).  For each eta:
  1. kz = 2*pi*j/Lz is scanned upward with a coarse eigensolve at every kz, until the growth
     rate has clearly passed its maximum -> candidate (kz, branch) pairs;
  2. each promising candidate is refined by sparse shift-invert with s increased until the
     eigenvalue changes by less than `tol` between successive s and the edge spectrum z(s) is small;
  3. the converged mode is verified by a fresh search for the rightmost eigenvalue;
  4. the Table-6 diagnostics (Brms ratios, helicities, Lorentz work, ...) are computed from
     the eigenvector of the fastest mode.
Many eta values run in parallel processes (one eta per process).  See README.md.

Validation (w0 = v0, Lz = 2 Lx): eta_crit = 0.501 (I), 0.146 (II), 0.254 (III), 0.220 (IV),
cf. Pencil 64^2x128: 0.495, 0.145, 0.255, 0.220. Flow IV, eta = 0.011: lambda = 0.0337,
omega_mod = 2|Im p| = 0.1398 (Table 3: 0.0336, 0.1397).  Flow I, eta = 0.05, w0 = sqrt(2) v0:
lambda = 0.153398401832 and all Table-6 quantities agree with Pencil 256^2x512 to 1e-11.

Usage:
    from solver import growth, scan, run_jobs
    r = growth('I', eta=0.01, tol=1e-5)            # box Lz = 4 pi by default
    rows = scan('I', np.logspace(-3, np.log10(0.6), 40), tol=1e-4, workers=8)
"""

import os, time, multiprocessing as _mp

# Child processes: one BLAS/LAPACK thread each (avoids oversubscription; must precede numpy import).
if _mp.parent_process() is not None:
    for _v in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_v] = "1"

import numpy as np
import scipy.sparse as sp, scipy.linalg as la, scipy.sparse.linalg as sla
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

# ----------------------------------------------------------------------------- operator
# Velocity fields u(x, y) = (ux, uy, uz), 2*pi-periodic in x and y, independent of z.
# To study another flow, add an entry here; nothing else in the code needs to change.
# Requirement: each component must be a trigonometric polynomial of degree <= 4 in x and y
# (see _coeffs).  v0 scales the horizontal and w0 the vertical component.
FLOWS = {
    # Roberts flows I-IV, eq. (2.1) of the paper
    "I": lambda x, y, v0, w0: (
        v0 * np.sin(x) * np.cos(y),
        -v0 * np.cos(x) * np.sin(y),
        w0 * np.sin(x) * np.sin(y),
    ),
    "II": lambda x, y, v0, w0: (
        v0 * np.sin(x) * np.cos(y),
        -v0 * np.cos(x) * np.sin(y),
        w0 * np.cos(x) * np.cos(y),
    ),
    "III": lambda x, y, v0, w0: (
        v0 * np.sin(x) * np.cos(y),
        -v0 * np.cos(x) * np.sin(y),
        w0 * 0.5 * (np.cos(2 * x) + np.cos(2 * y)),
    ),
    "IV": lambda x, y, v0, w0: (
        v0 * np.sin(x) * np.cos(y),
        -v0 * np.cos(x) * np.sin(y),
        w0 * np.sin(x),
    ),
    # Roberts (1972) motion I, u = (cos y - cos z, sin z, sin y) with his invariant direction x
    # renamed z (his y, z -> our x, y).  See README, "Original Roberts (1972) setup".
    "R72_I": lambda x, y, v0, w0: (
        v0 * np.sin(y),
        v0 * np.sin(x),
        w0 * (np.cos(x) - np.cos(y)),
    ),
    # Roberts (1972) motions II-IV, eqs. (6.1)-(6.3), same relabelling of axes:
    #   (cos y + cos z, sin z, sin y), (2 cos y cos z, sin z, sin y), (sin(y+z), sin 2z, sin 2y)
    "R72_II": lambda x, y, v0, w0: (
        v0 * np.sin(y),
        v0 * np.sin(x),
        w0 * (np.cos(x) + np.cos(y)),
    ),
    "R72_III": lambda x, y, v0, w0: (
        v0 * np.sin(y),
        v0 * np.sin(x),
        w0 * 2 * np.cos(x) * np.cos(y),
    ),
    "R72_IV": lambda x, y, v0, w0: (
        v0 * np.sin(2 * y),
        v0 * np.sin(2 * x),
        w0 * np.sin(x + y),
    ),
}


def velocity(name, x, y, v0=1.0, w0=1.0):
    """Velocity components (ux, uy, uz) of flow `name` at points (x, y)."""
    return FLOWS[name](x, y, v0, w0)


def vrms(name, v0=1.0, w0=1.0):
    """Root-mean-square velocity over the box (exact, from the Fourier coefficients)."""
    return np.sqrt(
        sum(
            sum(
                abs(c) ** 2
                for c in _coeffs(lambda x, y: velocity(name, x, y, v0, w0)[i]).values()
            )
            for i in range(3)
        )
    )


def _coeffs(f, n=16, amax=4):
    """Fourier coefficients of a 2*pi-periodic function f(x, y).

    Parameters
        f    : callable f(X, Y) on a grid (a flow component or a derivative of one)
        n    : grid points per direction for the FFT (16 is exact for the flows used here)
        amax : largest |wavenumber| kept
    Returns
        dict {(a, b): c_ab} with f = sum c_ab e^{i(a x + b y)}, only non-zero coefficients.
    """
    x = 2 * np.pi * np.arange(n) / n
    X, Y = np.meshgrid(x, x, indexing="ij")
    F = np.fft.fft2(f(X, Y)) / n**2
    return {
        (a, b): F[a % n, b % n]
        for a in range(-amax, amax + 1)
        for b in range(-amax, amax + 1)
        if abs(F[a % n, b % n]) > 1e-12
    }


def _mult(c, s):
    """Sparse matrix of multiplication by sum_ab c_ab e^{i(ax+by)} acting on the Fourier
    coefficients |m|,|n| <= s (index order: m slow, n fast).  Products that fall outside the
    truncation are dropped (Galerkin truncation)."""
    N = 2 * s + 1
    # sp.eye(N, k=-a) shifts the index m -> m + a, i.e. multiplication by e^{i a x}
    M = sp.csr_matrix((N * N, N * N), dtype=complex)  # zero if c is empty
    return sum(
        (v * sp.kron(sp.eye(N, k=-a), sp.eye(N, k=-b)) for (a, b), v in c.items()), M
    ).tocsr()


def build(name, eta, kz, s, w0=1.0, v0=1.0):
    """Sparse matrix of the operator p* acting on the (Hx, Hy) Fourier coefficients.

    Parameters
        name   : key of FLOWS ('I', 'II', 'III', 'IV', 'R72_I', ... or a flow you added)
        eta    : magnetic diffusivity (units v0/k0)
        kz     : vertical wavenumber of the mode
        s      : truncation, modes |m|,|n| <= s; matrix size 2(2s+1)^2
        w0, v0 : amplitudes of the vertical and horizontal flow
    Returns
        scipy.sparse CSC matrix A, eigenvalues of A are p* = p + eta kz^2.
    """
    ux, uy, uz = [
        _coeffs(lambda x, y, i=i: velocity(name, x, y, v0, w0)[i]) for i in range(3)
    ]
    d = lambda c, i: {
        k: v * 1j * k[i] for k, v in c.items() if k[i] != 0
    }  # d/dx (i=0), d/dy (i=1)
    N = 2 * s + 1
    m = np.arange(-s, s + 1)
    Dx = sp.kron(sp.diags(1j * m), sp.eye(N))
    Dy = sp.kron(sp.eye(N), sp.diags(1j * m))
    # terms acting on each component alike: -i kz uz H - (u.grad)H + eta lap2D H
    adv = (
        -1j * kz * _mult(uz, s)
        - _mult(ux, s) @ Dx
        - _mult(uy, s) @ Dy
        + eta * (Dx @ Dx + Dy @ Dy)
    )
    # stretching (H.grad)u couples Hx and Hy:  x-eq: Hx du_x/dx + Hy du_x/dy, y-eq likewise
    return sp.bmat(
        [
            [adv + _mult(d(ux, 0), s), _mult(d(ux, 1), s)],
            [_mult(d(uy, 0), s), adv + _mult(d(uy, 1), s)],
        ]
    ).tocsc()


def edge_z(vec, s):
    """Roberts' resolution diagnostic z(s): largest coefficient magnitude on the truncation edge
    max(|m|,|n|) = s, relative to the largest coefficient.  The eigenvalue error is ~ z^2.
    """
    N = 2 * s + 1
    a = np.abs(vec).reshape(2, N, N).max(axis=0)
    edge = np.r_[a[0, :], a[-1, :], a[:, 0], a[:, -1]].max()
    return edge / a.max()


# ----------------------------------------------------------------------------- eigen-solves
DENSE_MAX = 2000  # use dense LAPACK below this matrix size


def _top(ev, n):
    """The n eigenvalues with the largest real part, with near-duplicates (|diff| < 1e-6) removed."""
    out = []
    for e in np.asarray(ev)[np.argsort(-np.real(ev))]:
        if all(abs(e - o) > 1e-6 for o in out):
            out.append(e)
        if len(out) == n:
            break
    return out


def _shifts(name, eta, kz, w0, re0=0.3, m=10, s_lo=12, band=1.0, dw=0.5):
    """Shift points re0 + i*w for ARPACK shift-invert.

    The frequencies w are taken from the top-m modes of a cheap dense solve at s_lo (frequencies
    are set by advection and are already good at low s), plus a coarse uniform band |w| <= band
    (step dw) as a safety net.  re0 must exceed the largest growth rate (0.3 > ~0.15 for flows I-IV).
    """
    top = _top(la.eigvals(build(name, eta, kz, s_lo, w0).toarray()), m)
    ws = {round(float(e.imag) * 4) / 4 for e in top} | set(
        np.round(np.arange(-band, band + 1e-9, dw), 2)
    )
    return [re0 + 1j * w for w in sorted(ws)]


def _candidates(name, eta, kz, s, w0, n, k=6, shifts=None, **shift_kw):
    """Top-n eigenvalues p* (largest real part) at truncation s.

    Dense LAPACK for small matrices; otherwise ARPACK shift-invert at several shifts (ARPACK
    cannot target 'largest real part' of this strongly non-normal operator reliably).
    k: eigenvalues per shift; shifts: explicit list, else from _shifts(**shift_kw).
    """
    A = build(name, eta, kz, s, w0)
    if A.shape[0] <= DENSE_MAX:
        return _top(la.eigvals(A.toarray()), n)
    shifts = shifts if shifts is not None else _shifts(name, eta, kz, w0, **shift_kw)
    ev = []
    for sig in shifts:
        try:
            ev.extend(sla.eigs(A, k=k, sigma=sig, return_eigenvectors=False))
        except sla.ArpackNoConvergence as e:  # keep whatever converged
            ev.extend(e.eigenvalues)
    return _top(ev, n)


def _solve_near(name, eta, kz, s, w0, sigma):
    """Eigenvalue p* closest to sigma at truncation s, with its eigenvector (Hx, Hy coefficients)."""
    A = build(name, eta, kz, s, w0)
    if A.shape[0] <= DENSE_MAX:
        ev, V = la.eig(A.toarray())
    else:
        ev, V = sla.eigs(A, k=4, sigma=sigma, which="LM")
    i = np.argmin(abs(ev - sigma))
    return ev[i], V[:, i]


def s_guess(eta, c=1.6, smin=6):
    """Initial truncation s = c/sqrt(eta) (at least smin), from the eta^(1/2) thickness of the
    magnetic boundary layers (cf. Roberts 1972, table 5)."""
    return max(smin, int(np.ceil(c / np.sqrt(eta))))


def converge(
    name,
    eta,
    kz,
    w0=1.0,
    sigma=None,
    tol=1e-4,
    s0=None,
    smax=200,
    growth_factor=1.3,
    ztol=None,
    verify=True,
    cand_kw=None,
):
    """Follow one eigen-branch at fixed kz, increasing s until it is converged.

    Parameters
        sigma         : starting guess for p* (default: rightmost eigenvalue at s0)
        tol           : absolute tolerance, |p(s) - p(s_prev)| < tol (complex p, units v0 k0)
        s0, smax      : first and largest truncation; s grows by growth_factor (at least +2)
        ztol          : required z(s) (default sqrt(tol), since the error is ~ z^2)
        verify        : at the final s, search again for the rightmost eigenvalue at this kz;
                        if another branch grows faster (by more than tol), follow that one
                        instead (up to 3 times)
        cand_kw       : keyword arguments passed to _candidates
    Returns
        dict(p = growth rate + i*frequency, s, z, kz, converged, vec = eigenvector at s)
    """
    cand_kw = cand_kw or {}
    s = s0 or s_guess(eta)
    if sigma is None:
        sigma = _candidates(name, eta, kz, s, w0, 1, **cand_kw)[0]
    ztol = ztol if ztol is not None else np.sqrt(tol)  # |dp| ~ z^2 (Roberts eq. 10.6)
    for attempt in range(3):
        prev, res = None, None
        while s <= smax:
            p, v = _solve_near(name, eta, kz, s, w0, sigma)
            z = edge_z(v, s)
            if prev is not None and abs(p - prev) < tol and z < ztol:
                res = dict(p=p - eta * kz**2, s=s, z=z, kz=kz, converged=True, vec=v)
                break
            prev, sigma, vprev, sprev, zprev = (
                p,
                p,
                v,
                s,
                z,
            )  # continuation: next shift = p
            s = max(s + 2, int(np.ceil(s * growth_factor)))
        if res is None:  # smax reached: return last estimate
            return dict(
                p=prev - eta * kz**2,
                s=sprev,
                z=zprev,
                kz=kz,
                converged=False,
                vec=vprev,
            )
        if not verify:
            return res
        top = _candidates(name, eta, kz, res["s"], w0, 1, **cand_kw)[0]
        if (top - eta * kz**2).real <= res["p"].real + tol:
            return res
        sigma = top  # a faster branch exists: follow it
    return res


# ----------------------------------------------------------------------------- growth at one eta
def growth(
    name,
    eta,
    w0=1.0,
    Lz=4 * np.pi,
    tol=1e-4,
    kzmax=10.0,
    n_branch=2,
    n_kz=4,
    s_coarse=None,
    c_coarse=2.5,
    kz_margin=0.02,
    smax=200,
    verify=True,
    cand_kw=None,
    k1=None,
):
    """Fastest-growing mode at resistivity eta in a periodic box (allowed kz = 2 pi j / Lz).

    Parameters
        name, eta, w0 : flow, diffusivity, vertical flow amplitude (v0 = 1)
        Lz            : box height (Lx = Ly = 2 pi)
        tol           : absolute accuracy target for the complex growth rate p (units v0 k0)
        s_coarse      : truncation for the kz/branch search (default c_coarse/sqrt(eta), at
                        least 12); it must be fine enough to rank kz correctly (flow II at
                        eta = 0.01 needs ~20)
        kzmax         : hard upper limit of the kz scan.  The scan stops earlier once the
                        coarse growth rate has decreased at 3 successive kz and is below the best
                        by kz_margin.  The optimal kz grows as eta decreases (flow I: kz ~ 3.5
                        at eta = 2.6e-3)
        n_kz, n_branch: refined to tolerance are all kz within kz_margin of the best (at most
                        n_kz of them), with n_branch branches each
        smax, verify, cand_kw : passed to converge()
        k1            : wavenumber used to normalise Jrms (default: kz of the fastest mode)
    Returns
        dict(p, s, z, kz, converged, eta, flow, Rm = vrms/(eta k0), kz_scanned, and the
        Table-6 quantities listed in QUANTITIES)
    """
    cand_kw = cand_kw or {}
    sc = s_coarse or s_guess(eta, c=c_coarse, smin=12)
    dk = 2 * np.pi / Lz
    coarse, gk = [], []  # (growth estimate, kz, p*), best per kz
    j = 1
    while True:  # upward kz scan at coarse resolution
        kz = j * dk
        cands = _candidates(name, eta, kz, sc, w0, n_branch, **cand_kw)
        coarse += [((c - eta * kz**2).real, kz, c) for c in cands]
        gk.append(max(g for g, k_, _ in coarse if k_ == kz))
        if kz + dk > kzmax + 1e-9:
            break
        if j >= 3 and gk[-1] < gk[-2] < gk[-3] and gk[-1] < max(gk) - kz_margin:
            break  # clearly past the maximum in kz
        j += 1
    kz_last = kz
    coarse.sort(key=lambda r: -r[0])
    best_kz = []
    for g, kz, c in coarse:  # promising kz values (within kz_margin)
        if g < coarse[0][0] - kz_margin:
            break
        if kz not in best_kz:
            best_kz.append(kz)
        if len(best_kz) == n_kz:
            break
    results = [
        converge(
            name,
            eta,
            kz,
            w0,
            sigma=c,
            tol=tol,
            s0=sc,
            smax=smax,
            verify=verify,
            cand_kw=cand_kw,
        )
        for g, kz, c in coarse
        if kz in best_kz
    ]
    best = max(results, key=lambda r: r["p"].real)
    best.update(
        mode_quantities(name, best["kz"], best.pop("vec"), best["s"], w0, k1=k1)
    )
    for r in results:
        r.pop("vec", None)  # eigenvectors are not returned
    best.update(eta=eta, flow=name, Rm=vrms(name, w0=w0) / eta, kz_scanned=kz_last)
    return best


# ----------------------------------------------------------------------------- mode diagnostics
# Column names of the Table-6 quantities (volume averages, mu0 = 1).
QUANTITIES = [
    "Bmean_rms/Brms",
    "Jrms/(k1 Brms)",
    "<J.B>/(Jrms Brms)",
    "<(J.B)^2>^1/2/(Jrms Brms)",
    "<(JxB)^2>^1/2/(Jrms Brms)",
    "-W_L/(urms Jrms Brms)",
    "<J^2 B^2>/(Jrms^2 Brms^2)",
    "-<A.B>/Brms^2",
]


def mode_fields(name, kz, vec, s, w0=1.0, v0=1.0, nz=8):
    """Real 3D fields of one eigenmode, B = Re[b(x,y) exp(i kz z)].

    The grid is fine enough that volume averages of products of up to four fields are exact
    (x, y: 4s+4 points; z: nz points over one wavelength 2 pi/kz, which gives the box average
    since the box holds whole wavelengths).  Growth factor and time phase drop out of all
    normalised averages (for a single mode a time shift is a z shift).
    Returns
        B, J, A (Coulomb gauge), u : arrays (3, ng, ng, nz);  Bmean (xy-average) : (3, nz)
    """
    N = 2 * s + 1
    m = np.arange(-s, s + 1)
    M, Nn = np.meshgrid(m, m, indexing="ij")
    hx, hy = vec.reshape(2, N, N)
    hz = -(M * hx + Nn * hy) / kz  # div B = 0
    Bh = np.array([hx, hy, hz])
    K = np.array([M, Nn, np.full_like(M, kz, dtype=float)])
    Jh = 1j * np.cross(K, Bh, axis=0)  # J = curl B   (mu0 = 1)
    Ah = 1j * np.cross(K, Bh, axis=0) / np.sum(K**2, axis=0)  # curl A = B, div A = 0
    ng = 4 * s + 4

    def to_grid(c):  # sum_mn c_mn e^{i(mx+ny)} on grid
        F = np.zeros((ng, ng), complex)
        F[np.ix_(m % ng, m % ng)] = c
        return np.fft.ifft2(F) * ng * ng

    z = 2 * np.pi / kz * np.arange(nz) / nz
    ph = np.exp(1j * kz * z)
    real3 = lambda ch: np.array([np.real(to_grid(c)[:, :, None] * ph) for c in ch])
    x = 2 * np.pi * np.arange(ng) / ng
    X, Y = np.meshgrid(x, x, indexing="ij")
    u = np.array(velocity(name, X, Y, v0, w0))
    u = np.repeat(u[..., None], nz, axis=-1)
    Bmean = np.real(Bh[:, s, s][:, None] * ph)  # plane average = (0,0) mode
    return real3(Bh), real3(Jh), real3(Ah), u, Bmean


def mode_quantities(name, kz, vec, s, w0=1.0, k1=None):
    """Table-6 quantities (volume averages) for one eigenmode, as a dict keyed by QUANTITIES.
    k1 normalises Jrms (default kz).  W_L = <u.(J x B)> is the work done by the Lorentz force.
    """
    B, J, A, u, Bm = mode_fields(name, kz, vec, s, w0)
    avg = lambda f: f.mean()
    dot = lambda a, b: np.sum(a * b, axis=0)
    Brms, Jrms, urms = (
        np.sqrt(avg(dot(B, B))),
        np.sqrt(avg(dot(J, J))),
        np.sqrt(avg(dot(u, u))),
    )
    JB, JxB = dot(J, B), np.cross(J, B, axis=0)
    k1 = kz if k1 is None else k1
    vals = [
        np.sqrt(np.mean(dot(Bm, Bm))) / Brms,
        Jrms / (k1 * Brms),
        avg(JB) / (Jrms * Brms),
        np.sqrt(avg(JB**2)) / (Jrms * Brms),
        np.sqrt(avg(dot(JxB, JxB))) / (Jrms * Brms),
        -avg(dot(u, JxB)) / (urms * Jrms * Brms),
        avg(dot(J, J) * dot(B, B)) / (Jrms**2 * Brms**2),
        -avg(dot(A, B)) / Brms**2,
    ]
    return dict(zip(QUANTITIES, vals))


# ----------------------------------------------------------------------------- parallel scan
def _job(args):
    """Worker: run growth() for one (run_id, flow, eta, kwargs) tuple, record run id and time."""
    run_id, name, eta, kw = args
    t0 = time.time()
    r = growth(name, eta, **kw)
    r.update(run_id=run_id, seconds=time.time() - t0)
    return r


def _flat(r):
    """Result dict -> CSV-friendly row (complex p split into growth rate and frequency)."""
    return {
        "run_id": r.get("run_id"),
        "flow": r["flow"],
        "eta": r["eta"],
        "Rm": r["Rm"],
        "growth_rate": r["p"].real,
        "omega": r["p"].imag,
        "omega_mod": 2 * abs(r["p"].imag),
        "kz": r["kz"],
        "kz_scanned": r.get("kz_scanned"),
        "s": r["s"],
        "z": r["z"],
        "converged": r["converged"],
        "seconds": r.get("seconds"),
        **{q: r.get(q) for q in QUANTITIES},
    }


def _fmt_time(sec):
    """Seconds -> 'h:mm:ss'."""
    sec = int(sec)
    return f"{sec//3600}:{sec%3600//60:02d}:{sec%60:02d}"


def run_jobs(jobs, workers=None, out=None, progress=True, **kw):
    """Run growth() for a list of jobs [(run_id, flow, eta), ...] in parallel processes.

    Parameters
        workers  : number of processes (None -> all cores)
        out      : tab-separated results file; each finished job is appended immediately, so
                   partial results survive an interrupted run (an existing file is replaced)
        progress : print one line per finished job (elapsed time, lambda, omega, kz, s)
        **kw     : passed to growth() (e.g. tol, w0, Lz)
    Returns
        list of flat result rows, in the order of `jobs`; failed jobs carry an 'error' entry.
    """
    # Smallest eta first: these are the slowest jobs, starting them early balances the load.
    order = sorted(range(len(jobs)), key=lambda i: jobs[i][2])
    rows = [None] * len(jobs)
    t0 = time.time()
    if out is not None and os.path.exists(out):
        os.remove(out)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_job, (*jobs[i], kw)): i for i in order}
        if progress:
            print(f"{len(jobs)} jobs on {ex._max_workers} workers", flush=True)
        for done, f in enumerate(as_completed(futs), 1):
            i = futs[f]
            run_id, name, eta = jobs[i]
            try:
                rows[i] = _flat(f.result())
            except Exception as e:  # keep going; record the failure
                rows[i] = {"run_id": run_id, "flow": name, "eta": eta, "error": repr(e)}
            if out is not None:
                pd.DataFrame([rows[i]]).to_csv(
                    out, sep="\t", index=False, mode="a", header=not os.path.exists(out)
                )
            if progress:
                el = time.time() - t0
                r = rows[i]
                msg = (
                    (
                        f"lambda={r['growth_rate']:+.5f} omega={r['omega']:+.5f} kz={r['kz']:.3g} "
                        f"s={r['s']} {'ok' if r['converged'] else 'NOT CONVERGED'} "
                        f"({r['seconds']:.0f} s)"
                    )
                    if "error" not in r
                    else f"FAILED: {r['error']}"
                )
                print(
                    f"[{done:>{len(str(len(jobs)))}}/{len(jobs)}] elapsed {_fmt_time(el)}  "
                    f"run {run_id} flow {name} eta={eta:.4g}: {msg}",
                    flush=True,
                )
    return rows


def scan(name, etas, workers=None, **kw):
    """growth() for every eta of one flow, in parallel; returns the raw result dicts
    (no progress output, no file).  **kw is passed to growth()."""
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return [r for r in ex.map(_job, [(None, name, e, kw) for e in etas])]


# ----------------------------------------------------------------------------- parameter file
FLOW_CODES = {1: "I", 2: "II", 3: "III", 4: "IV"}  # Flow_kind column -> flow name


def read_parameter_csv(path, status="done"):
    """Read the SWIFT run table -> list of (run_id, flow_name, eta).

    Required columns: 'Run #', 'eta', 'Flow_kind' (1-4 or I-IV) and, unless status is None,
    'Status'.  The separator is auto-detected (tab or comma).  Rows whose status equals `status`
    (case- and space-insensitive; None -> all rows) and with a valid eta and flow are kept;
    skipped rows are reported by run number.
    """
    data = pd.read_csv(path, sep=None, engine="python")
    data.columns = data.columns.str.strip()
    if status is not None:
        data = data[
            data["Status"].astype(str).str.strip().str.lower() == status.lower()
        ]
    eta = pd.to_numeric(data["eta"], errors="coerce")
    flow = data["Flow_kind"].map(_flow_name)
    ok = eta.notna() & (eta > 0) & flow.notna()
    if (~ok).any():
        print(
            f"read_parameter_csv: skipping {(~ok).sum()} row(s) with missing/invalid eta or "
            f"Flow_kind: runs {data.loc[~ok, 'Run #'].tolist()}"
        )
    return list(zip(data.loc[ok, "Run #"], flow[ok], eta[ok].astype(float)))


def _flow_name(v):
    """Flow_kind value -> flow name (accepts 1-4, '1'-'4', 1.0, or a FLOWS key); None if invalid."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    sv = str(v).strip().upper()
    if sv in FLOWS:
        return sv
    try:
        return FLOW_CODES.get(int(float(sv)))
    except ValueError:
        return None


if __name__ == "__main__":
    t = time.time()
    jobs = read_parameter_csv("./test_run_parameters.csv")
    # tol = 1e-3 is loose near onset and in flow II (lambda ~ 0.005); 1e-4 is recommended there.
    rows = run_jobs(jobs, out="results.csv", tol=1e-4)
    print(
        f"done: {len(rows)} runs in {_fmt_time(time.time()-t)}, results in results.csv"
    )
