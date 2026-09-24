"""Re-derive the ISRF artificial-dissipation parameters for the flux-first
update order with e-weighted dissipation.

Shipped update, per band (e = exp(-a), phi = (1 - e)/a):

    F^{n+1} = e F^n - c_hyp^2 dt phi G[u^n]                    (extra ghost)
    u^{n+1} = e (u^n + dt phi diss[u^n]) + dt phi (r s - div F^{n+1})

Parameters under test:

    alpha_max  ISRF_dissipation_alpha_max             (0.5)
    eps_1      ISRF_dissipation_negativity_threshold  (0.01)
    guard      6.2 alpha C_hyp + 0.70 C_hyp^2 <= 2,   alpha = max(alpha_max,
               alpha_floor), checked at start-up.

Linearised model, conventions of verify_isrf_dissipation.py Part I: one
Fourier mode on a uniform lattice, state (u, F), transport number nu =
c_hyp K dt phi, dissipation number a_d = dt phi Gamma(k). Under the closure
c_hyp = C_hyp h/dt with x = h/lambda, a = C_hyp x, so a_d,max = G alpha
C_hyp phi(a) and nu_max = 1.18 C_hyp phi(a), with G = 2 I_W = 6.2 (every
neighbour out of phase, the enforced constant) or the lattice zone maximum.

Part A reproduces the constants and the guard of the previous
derivation (dissipation applied to the transported u, matrix G_r).
Part B proves that the flux-first e-weighted matrix G_A_ew has the same
characteristic polynomial as G_r at every (e, a_d, nu), derives its Schur
conditions, and shows e = 1 is the binding case under the closure.
Part C bisects the largest stable alpha on the closure scan for both
matrices. Part D computes the exact joint-k lattice bound. Part E measures
the transient growth max_n ||G^n||_2. Part F quantifies the e-weighting of
the dissipation in thick gas. Part G states what the linear model says
about eps_1.
"""

import numpy as np
import sympy as sp
from scipy.integrate import quad

GAMMA_3D = 1.936492  # kernel_gamma, Wendland C2, 3D (src/kernel_hydro.h)
ETA = 1.2348  # resolution_eta in every shipped SubgridRadiation example
ALPHA_SHIPPED = 0.5
EPS_1_SHIPPED = 0.01
C_HYP_SHIPPED = 0.5
GUARD_G = 6.2
GUARD_NU2 = 0.70


def guard_alpha(C_hyp):
    """Largest alpha the start-up guard admits at C_hyp.

    Parameters
    ----------
    C_hyp : float
        Closure coefficient of the hyperbolic speed.

    Returns
    -------
    float
        (2 - 0.70 C_hyp^2) / (6.2 C_hyp).
    """
    return (2.0 - GUARD_NU2 * C_hyp**2) / (GUARD_G * C_hyp)


def guard_C(alpha):
    """Largest C_hyp the start-up guard admits at alpha.

    Parameters
    ----------
    alpha : float
        Dissipation coefficient ceiling.

    Returns
    -------
    float
        Positive root of 0.70 C^2 + 6.2 alpha C - 2 = 0.
    """
    b = GUARD_G * alpha
    return (-b + np.sqrt(b**2 + 8.0 * GUARD_NU2)) / (2.0 * GUARD_NU2)


# ---------------------------------------------------------------------------
# Kernel and lattice symbols
# ---------------------------------------------------------------------------
def wc2_dwdr(r, H):
    """dW/dr of the 3D Wendland C2 kernel with support radius H.

    Parameters
    ----------
    r : numpy.ndarray
        Distances.
    H : float
        Support radius.

    Returns
    -------
    numpy.ndarray
        Radial derivative, zero outside the support.
    """
    r = np.atleast_1d(np.asarray(r, dtype=np.float64))
    q = r / H
    out = np.zeros_like(q)
    inside = q < 1.0
    qi = q[inside]
    out[inside] = (21.0 / (2.0 * np.pi)) * (-20.0 * qi * (1.0 - qi) ** 3) / H**4
    return out


def lattice_neighbours(h):
    """Cubic-lattice offsets (dx = 1) inside the support of smoothing length h.

    Parameters
    ----------
    h : float
        Smoothing length in units of the lattice spacing.

    Returns
    -------
    tuple of numpy.ndarray
        Offsets (N, 3) and their norms (N,).
    """
    H = GAMMA_3D * h
    n = int(np.ceil(H)) + 1
    g = np.arange(-n, n + 1)
    xx, yy, zz = np.meshgrid(g, g, g, indexing="ij")
    pos = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1).astype(np.float64)
    r = np.linalg.norm(pos, axis=1)
    keep = (r > 0) & (r < H)
    return pos[keep], r[keep]


def lattice_symbols(k_vecs, h=ETA):
    """Dimensionless dissipation and transport symbols on the lattice.

    Parameters
    ----------
    k_vecs : numpy.ndarray
        Wave vectors (M, 3) in units of 1/dx.
    h : float
        Smoothing length in units of dx.

    Returns
    -------
    tuple of numpy.ndarray
        Gamma(k) h / v_sig = h sum |W'| (1 - cos k.r), and |K(k)| h with
        K(k) = sum W' (r/|r|) sin(k.r), both of shape (M,).
    """
    pos, r = lattice_neighbours(h)
    wdr = wc2_dwdr(r, GAMMA_3D * h)
    phase = k_vecs @ pos.T
    gam = h * (np.abs(wdr)[None, :] * (1.0 - np.cos(phase))).sum(axis=1)
    kvec = (np.sin(phase) * wdr[None, :]) @ (pos / r[:, None])
    return gam, h * np.linalg.norm(kvec, axis=1)


# ---------------------------------------------------------------------------
# Amplification matrices, batched over the last two axes
# ---------------------------------------------------------------------------
def G_relocated(e, a_d, nu):
    """Previous order: u* = e u - i nu F, u' = (1 - a_d) u*, F' = e F - i nu u*.

    Parameters
    ----------
    e, a_d, nu : numpy.ndarray
        Absorption factor, dissipation number, transport number (broadcast).

    Returns
    -------
    numpy.ndarray
        Complex amplification matrices, shape (..., 2, 2).
    """
    e, a_d, nu = np.broadcast_arrays(e, a_d, nu)
    M = np.empty(e.shape + (2, 2), dtype=complex)
    M[..., 0, 0] = (1.0 - a_d) * e
    M[..., 0, 1] = -1j * (1.0 - a_d) * nu
    M[..., 1, 0] = -1j * nu * e
    M[..., 1, 1] = e - nu**2
    return M


def G_A_ew(e, a_d, nu):
    """Shipped order: F' = e F - i nu u, u' = e (1 - a_d) u - i nu F'.

    Parameters
    ----------
    e, a_d, nu : numpy.ndarray
        Absorption factor, dissipation number, transport number (broadcast).

    Returns
    -------
    numpy.ndarray
        Complex amplification matrices, shape (..., 2, 2).
    """
    e, a_d, nu = np.broadcast_arrays(e, a_d, nu)
    M = np.empty(e.shape + (2, 2), dtype=complex)
    M[..., 0, 0] = e * (1.0 - a_d) - nu**2
    M[..., 0, 1] = -1j * e * nu
    M[..., 1, 0] = -1j * nu
    M[..., 1, 1] = e
    return M


def G_A_plain(e, a_d, nu):
    """Flux-first order without e-weighting: u' = (e - a_d) u - i nu F'.

    Parameters
    ----------
    e, a_d, nu : numpy.ndarray
        Absorption factor, dissipation number, transport number (broadcast).

    Returns
    -------
    numpy.ndarray
        Complex amplification matrices, shape (..., 2, 2).
    """
    e, a_d, nu = np.broadcast_arrays(e, a_d, nu)
    M = np.empty(e.shape + (2, 2), dtype=complex)
    M[..., 0, 0] = e - a_d - nu**2
    M[..., 0, 1] = -1j * e * nu
    M[..., 1, 0] = -1j * nu
    M[..., 1, 1] = e
    return M


def spectral_radius(M):
    """Largest eigenvalue modulus of a batch of matrices.

    Parameters
    ----------
    M : numpy.ndarray
        Matrices, shape (..., 2, 2).

    Returns
    -------
    numpy.ndarray
        Spectral radii, shape (...).
    """
    return np.abs(np.linalg.eigvals(M)).max(axis=-1)


def closure_grid(G_const, alpha, C_hyp, x_values, n_frac):
    """(e, a_d, nu) samples on the closure, full (a_d, nu) rectangle.

    Parameters
    ----------
    G_const : float
        Zone constant of the dissipation symbol.
    alpha : float
        Dissipation coefficient ceiling.
    C_hyp : float
        Closure coefficient.
    x_values : numpy.ndarray
        h/lambda samples; 0 is the thin limit e = 1.
    n_frac : int
        Samples of each of a_d and nu inside their maxima.

    Returns
    -------
    tuple of numpy.ndarray
        e, a_d, nu of shape (len(x_values), n_frac, n_frac).
    """
    a = C_hyp * np.asarray(x_values, dtype=np.float64)
    e = np.exp(-a)
    phi = np.where(a > 0, -np.expm1(-a) / np.where(a > 0, a, 1.0), 1.0)
    fr = np.linspace(0.0, 1.0, n_frac)
    e3 = e[:, None, None] * np.ones((1, n_frac, n_frac))
    a_d = (G_const * alpha * C_hyp * phi)[:, None, None] * fr[None, :, None]
    nu = (1.18 * C_hyp * phi)[:, None, None] * fr[None, None, :]
    return np.broadcast_arrays(e3, a_d, nu)


X_SCAN = np.concatenate(([0.0], np.logspace(-6.0, 3.0, 181)))


def closure_rho(G_fn, G_const, alpha, C_hyp, n_frac=33):
    """Maximum spectral radius of G_fn along the closure.

    Parameters
    ----------
    G_fn : callable
        Matrix constructor.
    G_const : float
        Zone constant of the dissipation symbol.
    alpha : float
        Dissipation coefficient ceiling.
    C_hyp : float
        Closure coefficient.
    n_frac : int
        Samples of each of a_d and nu inside their maxima.

    Returns
    -------
    float
        Maximum spectral radius.
    """
    e, a_d, nu = closure_grid(G_const, alpha, C_hyp, X_SCAN, n_frac)
    return float(spectral_radius(G_fn(e, a_d, nu)).max())


def bisect_alpha(G_fn, G_const, C_hyp, tol_rho=1e-9):
    """Largest alpha with closure spectral radius <= 1 + tol_rho.

    Parameters
    ----------
    G_fn : callable
        Matrix constructor.
    G_const : float
        Zone constant of the dissipation symbol.
    C_hyp : float
        Closure coefficient.
    tol_rho : float
        Tolerance on the spectral radius.

    Returns
    -------
    float
        Bisected alpha, to 1e-7.
    """
    lo, hi = 0.0, 4.0
    assert closure_rho(G_fn, G_const, lo, C_hyp) <= 1.0 + tol_rho
    assert closure_rho(G_fn, G_const, hi, C_hyp) > 1.0 + tol_rho
    while hi - lo > 1e-7:
        mid = 0.5 * (lo + hi)
        if closure_rho(G_fn, G_const, mid, C_hyp) <= 1.0 + tol_rho:
            lo = mid
        else:
            hi = mid
    return lo


# ---------------------------------------------------------------------------
# Part A: previous derivation reproduced
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part A: previous derivation (dissipation on the transported u)")
print("=" * 78)


def _integrand_IW(r_):
    return 4.0 * np.pi * r_**2 * abs(wc2_dwdr(np.array([r_]), GAMMA_3D)[0])


I_W = quad(_integrand_IW, 0.0, GAMMA_3D, limit=200)[0]
print(f"  I_W = h int |W'| dV = {I_W:.4f}  (chapter 3.10); 2 I_W = {2 * I_W:.4f}")
assert abs(I_W - 3.10) < 5e-3
assert abs(2 * I_W - GUARD_G) < 1e-2

axes = np.array([[1.0, 0, 0], [1.0, 1.0, 0], [1.0, 1.0, 1.0]])
k_line = np.linspace(1e-3, np.pi, 400)
Kh_max = 0.0
for ax in axes:
    kh = ax / np.linalg.norm(ax)
    _, Kh = lattice_symbols(k_line[:, None] * kh[None, :])
    Kh_max = max(Kh_max, Kh.max())
print(f"  (K h)_max along axis/face/body diagonals = {Kh_max:.4f} (1.18)")
print(f"  (K h)_max^2 / 2 = {Kh_max**2 / 2:.4f} (enforced 0.70)")
assert abs(Kh_max - 1.18) < 0.02
assert abs(Kh_max**2 / 2 - GUARD_NU2) < 0.01

edge_gam, _ = lattice_symbols(np.array([[np.pi, np.pi, 0.0]]))
corner_gam, _ = lattice_symbols(np.array([[np.pi, np.pi, np.pi]]))
print(
    f"  lattice Gamma h / v_sig: corner {corner_gam[0]:.4f} (3.15), edge"
    f" {edge_gam[0]:.4f} (3.31)"
)
assert abs(corner_gam[0] - 3.15) < 0.05 and abs(edge_gam[0] - 3.31) < 0.05

e_s, ad_s, nu_s = sp.symbols("e a_d nu", real=True)
Gc = sp.Matrix([[1 - ad_s, -sp.I * nu_s], [-sp.I * nu_s * (1 - ad_s), 1 - nu_s**2]])
tr_c, det_c = sp.expand(Gc.trace()), sp.expand(Gc.det())
assert sp.simplify(tr_c - (2 - ad_s - nu_s**2)) == 0
assert sp.simplify(det_c - (1 - ad_s)) == 0
print("  density-loop placement at e = 1: tr = 2 - a_d - nu^2, det = 1 - a_d")
print("  Schur: |det| <= 1 and 1 + tr + det >= 0  <=>  a_d + nu^2/2 <= 2")
print("  with a_d,max = 2 I_W alpha C_hyp, nu_max = 1.18 C_hyp:")
print("        6.2 alpha C_hyp + 0.70 C_hyp^2 <= 2")
old = {
    "alpha(C=0.5)": guard_alpha(0.5),
    "C*(alpha=0.5)": guard_C(0.5),
    "C*(alpha=0)": np.sqrt(2.0 / GUARD_NU2),
}
for k_, v_ in old.items():
    print(f"  guard {k_} = {v_:.4f}")
assert abs(old["alpha(C=0.5)"] - 0.5887) < 1e-4
assert abs(old["C*(alpha=0.5)"] - 0.5714) < 1e-4
assert abs(old["C*(alpha=0)"] - 1.6903) < 1e-4
print(
    f"  shipped alpha_max = {ALPHA_SHIPPED} sits at {ALPHA_SHIPPED / old['alpha(C=0.5)']:.3f}"
    " of the guard at C_hyp = 0.5 (modal headroom"
    f" {old['alpha(C=0.5)'] / ALPHA_SHIPPED:.3f}x)"
)
print()

# ---------------------------------------------------------------------------
# Part B: the flux-first e-weighted map, exact invariants and Schur conditions
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part B: G_A_ew invariants and Schur conditions, symbolic")
print("=" * 78)
Gr_s = sp.Matrix(
    [
        [(1 - ad_s) * e_s, -sp.I * (1 - ad_s) * nu_s],
        [-sp.I * nu_s * e_s, e_s - nu_s**2],
    ]
)
Gew_s = sp.Matrix(
    [[e_s * (1 - ad_s) - nu_s**2, -sp.I * e_s * nu_s], [-sp.I * nu_s, e_s]]
)
lam = sp.symbols("lambda")
cp_r = sp.expand(Gr_s.charpoly(lam).as_expr())
cp_ew = sp.expand(Gew_s.charpoly(lam).as_expr())
assert sp.simplify(cp_r - cp_ew) == 0
T_ew = sp.expand(Gew_s.trace())
D_ew = sp.expand(Gew_s.det())
assert sp.simplify(T_ew - (e_s * (2 - ad_s) - nu_s**2)) == 0
assert sp.simplify(D_ew - e_s**2 * (1 - ad_s)) == 0
print("  char. polynomial of G_A_ew == that of G_r at every (e, a_d, nu):")
print("    tr = e (2 - a_d) - nu^2,   det = e^2 (1 - a_d)")
lower = sp.factor(sp.expand(1 + T_ew + D_ew))
upper = sp.factor(sp.expand(1 - T_ew + D_ew))
print(f"  1 + tr + det = {lower}")
print(f"  1 - tr + det = {upper}")
assert (
    sp.simplify(1 + T_ew + D_ew - ((1 + e_s) ** 2 - e_s * (1 + e_s) * ad_s - nu_s**2))
    == 0
)
print("  Schur conditions (spectral radius <= 1):")
print("    (S1) e a_d (1 + e) + nu^2 <= (1 + e)^2")
print("    (S2) (1 - e)^2 + e (1 - e) a_d + nu^2 >= 0   (always, a_d >= 0)")
print("    (S3) e^2 (1 - a_d) >= -1  <=>  a_d <= 1 + 1/e^2")
print("  at e = 1, (S1) is a_d + nu^2/2 <= 2: the density-loop condition.")

# Under the closure a_d = G alpha C phi(a), nu = 1.18 C phi(a), e = exp(-a).
# (S1)/(1+e)^2 = G alpha C [e phi/(1+e)] + (1.18 C)^2 [phi/(1+e)]^2, and both
# brackets decrease monotonically in a, so a = 0 binds.
a_grid = np.concatenate(([0.0], np.logspace(-8, 3, 4000)))
e_g = np.exp(-a_grid)
phi_g = np.where(a_grid > 0, -np.expm1(-a_grid) / np.where(a_grid > 0, a_grid, 1), 1.0)
b1 = e_g * phi_g / (1 + e_g)
b2 = (phi_g / (1 + e_g)) ** 2
assert np.all(np.diff(b1) <= 1e-15) and np.all(np.diff(b2) <= 1e-15)
assert abs(b1[0] - 0.5) < 1e-12 and abs(b2[0] - 0.25) < 1e-12
print("  under the closure, e phi/(1+e) and (phi/(1+e))^2 decrease monotonically")
print("  in a = C_hyp h/lambda (checked on 4000 points, a in [0, 1e3]):")
print("  the thin limit e = 1 binds, so the A_ew guard is")
print("        2 I_W alpha C_hyp + (K h)_max^2 C_hyp^2 / 2 <= 2")
print("  identical to the shipped 6.2 alpha C_hyp + 0.70 C_hyp^2 <= 2.")
print()

# ---------------------------------------------------------------------------
# Part C: bisected alpha_max on the closure scan
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part C: largest alpha with spectral radius <= 1 + 1e-9 on the closure")
print("=" * 78)
print(f"  h/lambda in {{0}} U [1e-6, 1e3] (182 values), 33 x 33 (a_d, nu) rectangle")
print("  G      C_hyp   alpha G_r    alpha G_A_ew  alpha G_A_plain  analytic")
results_C = {}
for G_const in (GUARD_G, float(edge_gam[0])):
    for C_hyp in (0.25, 0.5, 1.0):
        nu2 = (1.18 * C_hyp) ** 2 / 2
        analytic = (2.0 - nu2) / (G_const * C_hyp)
        ar = bisect_alpha(G_relocated, G_const, C_hyp)
        aew = bisect_alpha(G_A_ew, G_const, C_hyp)
        apl = bisect_alpha(G_A_plain, G_const, C_hyp)
        results_C[(G_const, C_hyp)] = (ar, aew, apl, analytic)
        print(
            f"  {G_const:4.2f}   {C_hyp:4.2f}    {ar:9.5f}    {aew:9.5f}     {apl:9.5f}"
            f"      {analytic:9.5f}"
        )
        assert abs(aew - ar) < 1e-6
        assert abs(aew - analytic) < 1e-5
        assert abs(apl - analytic) < 1e-5
a_ew_half = results_C[(GUARD_G, 0.5)][1]
print(
    f"  shipped point: alpha_max(C_hyp = 0.5) = {a_ew_half:.4f} under A_ew;"
    f" shipped 0.5 has {a_ew_half / ALPHA_SHIPPED:.3f}x modal headroom."
)
rho_059 = closure_rho(G_A_ew, GUARD_G, 0.59, 0.5)
rho_050 = closure_rho(G_A_ew, GUARD_G, 0.50, 0.5)
print(
    f"  closure rho at alpha = 0.50: {rho_050:.6f}; at 0.59 (just past the guard"
    f" {old['alpha(C=0.5)']:.4f}): {rho_059:.6f}"
)
assert rho_050 <= 1.0 + 1e-9 and rho_059 > 1.0
print()

# ---------------------------------------------------------------------------
# Part D: exact joint-k lattice bound
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part D: joint-k lattice bound, max_k [alpha C Gamma h + C^2 (K h)^2/2] <= 2")
print("=" * 78)
n_k = 33
g1 = np.linspace(0.0, np.pi, n_k)
kx, ky, kz = np.meshgrid(g1, g1, g1, indexing="ij")
k_all = np.stack([kx.ravel(), ky.ravel(), kz.ravel()], axis=1)[1:]
gam_k, Kh_k = lattice_symbols(k_all)
print(
    f"  {n_k}^3 zone grid: max Gamma h = {gam_k.max():.4f}, max |K| h ="
    f" {Kh_k.max():.4f} (off-axis |K| included)"
)


def lattice_alpha(C_hyp, nu2_scale=1.0):
    """Exact lattice alpha*(C_hyp) of the joint-k Schur condition.

    Parameters
    ----------
    C_hyp : float
        Closure coefficient.
    nu2_scale : float
        Factor on nu^2: 1 for free streaming (f = 1), 1/3 for the isotropic
        M1 limit.

    Returns
    -------
    float
        min_k (2 - C^2 (K h)^2/2) / (C Gamma h).
    """
    budget = 2.0 - nu2_scale * C_hyp**2 * Kh_k**2 / 2.0
    return float(np.min(budget / (C_hyp * gam_k)))


def lattice_C(alpha, nu2_scale=1.0):
    """Exact lattice C*(alpha), by bisection on lattice_alpha.

    Parameters
    ----------
    alpha : float
        Dissipation coefficient ceiling.
    nu2_scale : float
        Factor on nu^2, as in lattice_alpha.

    Returns
    -------
    float
        Largest C_hyp with lattice_alpha(C_hyp) >= alpha.
    """
    lo, hi = 1e-6, 10.0
    while hi - lo > 1e-6:
        mid = 0.5 * (lo + hi)
        if lattice_alpha(mid, nu2_scale) >= alpha:
            lo = mid
        else:
            hi = mid
    return lo


# Stability grid of 2026-09-14 (previous order, alpha pinned, isotropic
# pulse at Z = 1e-4): largest stable and smallest unstable C_hyp.
measured_grid = {0.25: (2.356, 2.531), 0.50: (1.188, 1.276), 0.75: (0.800, 0.840)}
print("  alpha  guard C*  lattice C* (f=1)  lattice C* (f=0)  measured bracket")
for alpha, (m_lo, m_hi) in measured_grid.items():
    cg = guard_C(alpha)
    cl = lattice_C(alpha)
    ci = lattice_C(alpha, 1.0 / 3.0)
    print(
        f"  {alpha:4.2f}   {cg:7.4f}   {cl:7.4f} ({cl / cg:4.2f}x)    {ci:7.4f} ({ci / cg:4.2f}x)"
        f"     ({m_lo:5.3f}, {m_hi:5.3f}] ({m_lo / cg:4.2f}-{m_hi / cg:4.2f}x)"
    )
    assert cl >= cg
    if alpha >= 0.5:
        assert m_lo * 0.95 <= ci <= m_hi * 1.05
print(
    f"  alpha*(C_hyp=0.5): guard {guard_alpha(0.5):.4f}, lattice f=1 {lattice_alpha(0.5):.4f}"
)
assert lattice_alpha(0.5) > 2.0 * guard_alpha(0.5)
print("  The exact lattice symbol puts the boundary inside the measured bracket")
print("  at alpha = 0.5 and 0.75 (isotropic, f = 0): the guard's 2.0x margin there")
print("  is its constant 6.2 = 2 I_W against the lattice maximum 3.31, not")
print("  unmodelled physics. At alpha = 0.25 the measured boundary lies above")
print("  the linear lattice bound, which the linear model does not explain.")
print("  Part B makes this lattice bound order-independent, so it holds for A_ew")
print("  in linear theory; the grid itself ran on the previous order. Unequal")
print("  h across a pair raises the dissipation sum (verify_isrf_dissipation.py")
print("  Part F: 1.22x at h_i/h_j = 1.44, 1.58x at 2.15), which uses part of it.")
print()

# ---------------------------------------------------------------------------
# Part E: transient growth
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part E: transient growth max_n ||G^n||_2 along the closure")
print("=" * 78)


def transient(G_fn, G_const, alpha, C_hyp, n_steps):
    """Transient growth and the step at which it peaks.

    Uses the grid of verify_isrf_dissipation.py Part I: 61 h/lambda values
    log-spaced over [1e-3, 1e3], 9 x 9 (a_d, nu) fractions.

    Parameters
    ----------
    G_fn : callable
        Matrix constructor.
    G_const : float
        Zone constant of the dissipation symbol.
    alpha : float
        Dissipation coefficient ceiling.
    C_hyp : float
        Closure coefficient.
    n_steps : int
        Number of powers probed.

    Returns
    -------
    tuple
        (max norm, argmax n).
    """
    e, a_d, nu = closure_grid(G_const, alpha, C_hyp, np.logspace(-3.0, 3.0, 61), 9)
    M = G_fn(e, a_d, nu).reshape(-1, 2, 2)
    P = np.broadcast_to(np.eye(2, dtype=complex), M.shape).copy()
    best, best_n = 0.0, 0
    for n in range(1, n_steps + 1):
        P = P @ M
        v = float(np.linalg.norm(P, ord=2, axis=(1, 2)).max())
        if v > best:
            best, best_n = v, n
    return best, best_n


design_rows = {
    (6.2, 0.5, 0.5): (1.3487, 1.5403),
    (6.2, 0.59, 0.5): (1.5978, 1.8986),
    (3.15, 1.16, 0.5): (1.3487, 1.6829),
    (6.2, 0.2, 1.0): (1.9128, 2.5349),
}
print("  (G, alpha, C_hyp)    G_r (n*)        A_ew (n*)      ratio  | n <= 1000")
for (G_const, alpha, C_hyp), (ref_r, ref_ew) in design_rows.items():
    tr_r, n_r = transient(G_relocated, G_const, alpha, C_hyp, 200)
    tr_ew, n_ew = transient(G_A_ew, G_const, alpha, C_hyp, 200)
    tr_ew_long, n_ew_long = transient(G_A_ew, G_const, alpha, C_hyp, 1000)
    alpha_lim = (2.0 - (1.18 * C_hyp) ** 2 / 2) / (G_const * C_hyp)
    inside = alpha <= alpha_lim
    print(
        f"  ({G_const:4.2f}, {alpha:4.2f}, {C_hyp:3.1f})  {tr_r:7.4f} ({n_r:3d})"
        f"   {tr_ew:7.4f} ({n_ew:3d})   {tr_ew / tr_r:5.3f}  | {tr_ew_long:7.4f} ({n_ew_long})"
        f"  {'inside' if inside else 'PAST'} the bound {alpha_lim:.4f}"
    )
    assert abs(tr_r - ref_r) < 2e-3 and abs(tr_ew - ref_ew) < 2e-3
    assert 1.10 < tr_ew / tr_r < 1.35
    if inside:
        assert n_ew < 150, "transient peak too close to the n = 200 cap"
        assert abs(tr_ew_long - tr_ew) < 1e-9
    if (G_const, alpha) == (6.2, 0.59):
        assert tr_ew_long > 2.0 * tr_ew
print("  The (6.2, 0.59, 0.5) row lies past the bound: its 200-step value is a")
print("  truncation of an unbounded growth, not a transient bound.")

N_TR = 60
print(f"\n  alpha dependence at C_hyp = 0.5, G = 6.2 (n <= {N_TR}):")
print("  alpha    G_r      A_ew     ratio   n*(A_ew)")
tr_alpha = []
for alpha in (0.0, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, guard_alpha(0.5)):
    tr_r, _ = transient(G_relocated, GUARD_G, alpha, 0.5, N_TR)
    tr_ew, n_ew = transient(G_A_ew, GUARD_G, alpha, 0.5, N_TR)
    tr_alpha.append((alpha, tr_r, tr_ew))
    print(f"  {alpha:6.4f}  {tr_r:7.4f}  {tr_ew:7.4f}  {tr_ew / tr_r:5.3f}   {n_ew}")
    assert n_ew < N_TR - 10
tr_ew_vals = np.array([t[2] for t in tr_alpha])
assert np.all(np.diff(tr_ew_vals) >= -1e-9)
assert abs(tr_alpha[0][2] - tr_alpha[0][1]) < 1e-9
assert abs(tr_alpha[1][2] - tr_alpha[1][1]) < 1e-9

lo, hi = 0.25, 0.5
while hi - lo > 1e-3:
    mid = 0.5 * (lo + hi)
    if transient(G_A_ew, GUARD_G, mid, 0.5, 20)[0] > tr_alpha[0][1] + 1e-6:
        hi = mid
    else:
        lo = mid
alpha_cross = lo
print(
    f"  A_ew transient equals G_r's ({tr_alpha[0][1]:.4f}) for alpha <= {alpha_cross:.3f}"
    " and grows with alpha above it; the excess is one-step (n* = 1),"
    " driven by the dissipation, zero with it off."
)

print("\n  transient at the guard boundary, alpha = (2 - 0.70 C^2)/(6.2 C):")
print("  C_hyp  alpha    G_r      A_ew     ratio")
for C_hyp in (0.25, 0.5, 1.0):
    ag = guard_alpha(C_hyp)
    tr_r, _ = transient(G_relocated, GUARD_G, ag, C_hyp, N_TR)
    tr_ew, n_ew = transient(G_A_ew, GUARD_G, ag, C_hyp, N_TR)
    print(f"  {C_hyp:4.2f}  {ag:6.4f}  {tr_r:7.4f}  {tr_ew:7.4f}  {tr_ew / tr_r:5.3f}")
    assert n_ew < N_TR - 10
print()

# ---------------------------------------------------------------------------
# Part F: e-weighting in thick gas
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part F: dissipation delivered per step in thick gas (null mode nu = 0)")
print("=" * 78)
print("  u' = e (1 - a_d) u under A_ew and G_r; u' = (e - a_d) u without e-weighting.")
print("  Extra decrement from dissipation: A_ew e a_d, unweighted a_d.")
print(
    "  a = C_hyp h/lambda   e        a_d (G=6.2, alpha=0.5, C=0.5)  e a_d    alpha to match"
)
for a in (0.01, 0.1, 0.3, 1.0, 3.0):
    e = np.exp(-a)
    phi = -np.expm1(-a) / a
    a_d = GUARD_G * ALPHA_SHIPPED * C_HYP_SHIPPED * phi
    print(
        f"  {a:5.2f}               {e:6.4f}   {a_d:6.4f}                          "
        f"{e * a_d:6.4f}   {ALPHA_SHIPPED / e:6.4f}"
    )
a_comp = np.log(old["alpha(C=0.5)"] / ALPHA_SHIPPED)
print(
    f"  raising alpha_max to the guard compensates the e factor only for a <="
    f" {a_comp:.3f}; beyond it the guard forbids the compensation. The"
    " previous order G_r already delivered e a_d (same invariants), so this is"
    " no change against the order it replaces."
)
e_t, a_t = 0.3, 0.7
assert abs(G_A_ew(e_t, a_t, 0.0)[0, 0] - G_relocated(e_t, a_t, 0.0)[0, 0]) < 1e-15
print()

# ---------------------------------------------------------------------------
# Part G: eps_1
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part G: eps_1 (ISRF_dissipation_negativity_threshold)")
print("=" * 78)
print(
    f"  eps_1 = {EPS_1_SHIPPED}: saturation scale of alpha_aim = alpha_max x^2 (3 - 2x),"
)
print("  x = min(eps/eps_1, 1), eps = max(-U, 0)/max(<|U|>, -U). It sets where the")
print("  coefficient reaches alpha_max, never how large it may be, so it is absent")
print("  from the linear map: every matrix above depends on alpha only through")
print("  a_d <= G alpha_max C_hyp phi. Under the flux-first order the trigger")
print("  reads U = rho u^n against the kernel mean of u^n (same time level),")
print("  the state the dissipation acts on. No stability re-derivation of")
print("  eps_1 exists; its evidence stays the 2026-09-10 bracket (0.001/0.01/0.1/")
print("  0.5, monotone, 0.01 on the favourable side) and the measured saturation")
print("  three orders past it.")
print()
print("ALL CHECKS PASSED")
