"""Verify WHICH differential operator the Design B `grad(u)` estimator
(design-lw-fuv-design-b.md Section 2.2) must target, and which of two SPH
constructions actually converges to it.

Section 1's flux (first-moment) equation, written for the mass-specific
variables `u = u_V/rho`, `F = F_V/rho` the implementation tracks, is

    DF/Dt = -(D/tau) * (1/rho) grad(rho u) - F/tau,

so the pressure-gradient term is `(1/rho)*grad(rho*u)` = `(1/rho)*grad(u_V)`:
the gradient of the VOLUMETRIC energy density, divided by rho, exactly as
SPH's own momentum equation needs `(1/rho)*grad(P)`, not `grad(P/rho)`.
The original Section 2.2 draft mirrored SPHENIX's `div_v` "difference"
form for `grad(u)`, which converges to plain `grad(u)` = `grad(u_V/rho)`:
the wrong operator whenever `grad(rho) != 0`. The two differ by
`u*grad(ln rho)`, and the difference form's error is not small: for a
UNIFORM radiation energy density `u_V` sitting on a static density
gradient (nothing should happen), it produces a spurious flux
`F = -D grad(u) = +D u_V grad(rho)/rho^2`, pointing toward the dense side,
so the scheme would pump radiation energy into dense gas out of nothing.

This script (a) states the algebra with sympy, (b) runs that exact
exposing configuration (uniform `u_V`, linear `rho(x)`, `F = 0`) on a
regular lattice for both constructions, and (c) repeats
`verify_design_b_div_f_consistency.py`'s neighbour-count convergence
sweep with linear `u(x)`, `rho(x)` fields to identify each construction's
continuum-limit operator. The symmetric construction is
`src/rt/SPHM1RT/rt_gradients.h`'s `radiation_gradient_SPH` `diffmode==1`
branch (the scalar twin of the `div(F)` construction Section 2.2 already
adopted, same shared-value / mirrored-mass-and-sign pattern), whose own
doxygen states it computes `(1/rho)*grad(rho*uin)`.
"""
# =============================================================================
# M1 CLOSURE AUDIT, 2026-09-11: CHECKED, CLOSURE-INDEPENDENT, NO CHANGE.
#
# The question this script settles -- which differential operator the
# gradient estimator converges to, `(1/rho)grad(rho u)` rather than
# `grad(u)` -- is unchanged by the closure. The shipped M1 operator
# contracts exactly the same volumetric field `rho*u` with a per-particle
# tensor before differencing it (`(1/rho)div(D*rho*u)`), and reduces to the
# form tested here, times 1/3, at D = I/3. The spurious-flux failure mode
# demonstrated below (uniform u_V on a density gradient) is likewise
# closure-independent: at F = 0 the closure is exactly isotropic.
# =============================================================================
import numpy as np
import sympy as sp
from scipy.spatial import cKDTree

RNG = np.random.default_rng(20260907)

# ---------------------------------------------------------------------------
# Part 0: the algebra
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part 0: which operator the mass-specific flux equation needs")
print("=" * 78)
x, y, z = sp.symbols("x y z")
rho = sp.Function("rho")(x, y, z)
u = sp.Function("u")(x, y, z)
uV = sp.Symbol("u_V", positive=True)  # a uniform radiation energy density
grad = lambda f: sp.Matrix([sp.diff(f, c_) for c_ in (x, y, z)])
correct = grad(rho * u) / rho
assert sp.simplify(correct - (grad(u) + u * grad(sp.log(rho)))) == sp.zeros(3, 1)
print("  (1/rho) grad(rho u) = grad(u) + u grad(ln rho)      [identity]")
u_uniform_uV = uV / rho
assert sp.simplify(grad(rho * u_uniform_uV) / rho) == sp.zeros(3, 1)
spurious = sp.simplify(grad(u_uniform_uV))
print("  uniform u_V (u = u_V/rho):  (1/rho) grad(rho u) = 0 exactly;")
print(f"                              grad(u) = {spurious.T} != 0")
print("  so the difference form drives F = -D grad(u) = +D u_V grad(rho)/rho^2,")
print("  a spurious flux toward the dense side out of a uniform field.")
print()


# ---------------------------------------------------------------------------
# Kernel (Wendland C2, 3D) and the two pairwise constructions
# ---------------------------------------------------------------------------
def wendland_c2_3d(q):
    q = np.asarray(q, dtype=np.float64)
    norm = 21.0 / (2.0 * np.pi)
    inside = q < 1.0
    w = np.zeros_like(q)
    dwdq = np.zeros_like(q)
    qi = q[inside]
    w[inside] = norm * (1.0 - qi) ** 4 * (4.0 * qi + 1.0)
    dwdq[inside] = norm * (-4.0 * (1.0 - qi) ** 3 * (4.0 * qi + 1.0) + 4.0 * (1.0 - qi) ** 4)
    return w, dwdq


def kernel_dr(r, h):
    """`wi_dr = h^-(dim+1) * dW/dq`, SWIFT's convention (dim=3)."""
    _, dwdq = wendland_c2_3d(r / h)
    return dwdq / h**4


def grad_pair_difference(dx, r, wi_dr, wj_dr, rhoi, rhoj, mi, mj, ui, uj):
    """The ORIGINAL §2.2 form (SPHENIX div_v mirrored):
    grad_i += m_j (u_j - u_i) wi_dr dx/r, finalized by 1/rho_i (folded in
    here); grad_j the mirror with its own kernel and density."""
    r_inv = 1.0 / r
    gi = (mj * (uj - ui) * wi_dr * r_inv / rhoi)[:, None] * dx
    gj = (mi * (ui - uj) * wj_dr * r_inv / rhoj)[:, None] * (-dx)
    return gi, gj


def grad_pair_symmetric(dx, r, wi_dr, wj_dr, rhoi, rhoj, mi, mj, ui, uj):
    """radiation_gradient_SPH diffmode==1: one shared scalar
    (u_i/rho_i wi_dr + u_j/rho_j wj_dr)/r, applied as +m_j dx to i and
    -m_i dx to j. No finalize. Right operator, but NOT the adjoint of the
    diffmode==1 divergence: the staggered scheme built on it is unstable on
    a disordered distribution (verify_design_b_timestepping_stability.py,
    Part F.3). Kept as the rejected candidate."""
    r_inv = 1.0 / r
    shared = (ui / rhoi * wi_dr + uj / rhoj * wj_dr) * r_inv
    gi = (mj * shared)[:, None] * dx
    gj = (-mi * shared)[:, None] * dx
    return gi, gj


def grad_pair_difference_rho_u(dx, r, wi_dr, wj_dr, rhoi, rhoj, mi, mj, ui, uj):
    """radiation_gradient_SPH diffmode==0 (ADOPTED): the difference form
    applied to rho*u with each particle's own kernel and 1/rho^2:
      grad_i += -m_j (rho_i u_i - rho_j u_j) wi_dr / (r rho_i^2) * dx
      grad_j += -m_i (rho_i u_i - rho_j u_j) wj_dr / (r rho_j^2) * dx
    Minus the adjoint of the diffmode==1 divergence in the m*rho inner
    product, so the transport conserves the discrete radiation energy."""
    r_inv = 1.0 / r
    d = rhoi * ui - rhoj * uj
    gi = (-mj * d * wi_dr * r_inv / rhoi**2)[:, None] * dx
    gj = (-mi * d * wj_dr * r_inv / rhoj**2)[:, None] * dx
    return gi, gj


CONSTRUCTIONS = {
    "difference (original §2.2 grad(u))": grad_pair_difference,
    "symmetric diffmode==1 (rejected)": grad_pair_symmetric,
    "difference on rho*u, diffmode==0 (adopted)": grad_pair_difference_rho_u,
}


def build_lattice(n_per_dim, box_l):
    dx_p = box_l / n_per_dim
    c1 = (np.arange(n_per_dim) + 0.5) * dx_p
    xx, yy, zz = np.meshgrid(c1, c1, c1, indexing="ij")
    return np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1), dx_p


def accumulate(pos, rho, mass, u_field, h, fn):
    tree = cKDTree(pos)
    pairs = tree.query_pairs(r=h, output_type="ndarray")
    i_idx, j_idx = pairs[:, 0], pairs[:, 1]
    dx = pos[i_idx] - pos[j_idx]
    r = np.linalg.norm(dx, axis=1)
    w = kernel_dr(r, h)
    gi, gj = fn(dx, r, w, w, rho[i_idx], rho[j_idx], mass[i_idx], mass[j_idx], u_field[i_idx], u_field[j_idx])
    out = np.zeros_like(pos)
    np.add.at(out, i_idx, gi)
    np.add.at(out, j_idx, gj)
    return out


def interior_probe(pos, h, box_l, frac):
    margin = 1.5 * h
    inside = np.all((pos > margin) & (pos < box_l - margin), axis=1)
    cand = np.where(inside)[0]
    target = np.array([frac * box_l, box_l / 2, box_l / 2])
    return cand[np.argmin(np.sum((pos[cand] - target) ** 2, axis=1))]


# ---------------------------------------------------------------------------
# Part 1: the exposing configuration
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part 1: uniform u_V on a static density gradient, F = 0: nothing should happen")
print("=" * 78)
rho0, slope, uV0, box_l, n_per_dim = 1.0, 1.5, 2.0, 1.0, 32
pos, dx_p = build_lattice(n_per_dim, box_l)
rho = rho0 * (1.0 + slope * pos[:, 0])
mass = rho * dx_p**3
u_spec = uV0 / rho  # uniform u_V
h_sweep = [1.3, 1.8, 2.6, 3.8, 5.5]
err_diff, err_sym, err_dru = [], [], []
for hod in h_sweep:
    h = hod * dx_p
    p_idx = interior_probe(pos, h, box_l, 0.5)
    analytic_spurious = -uV0 * rho0 * slope / rho[p_idx] ** 2  # x-component of grad(u)
    gd = accumulate(pos, rho, mass, u_spec, h, grad_pair_difference)[p_idx]
    gs = accumulate(pos, rho, mass, u_spec, h, grad_pair_symmetric)[p_idx]
    gr = accumulate(pos, rho, mass, u_spec, h, grad_pair_difference_rho_u)[p_idx]
    nn = (4.0 / 3.0) * np.pi * hod**3
    print(f"  ~{nn:5.0f} nbrs: difference grad_x = {gd[0]:+.5f} (grad(u)_x analytic {analytic_spurious:+.5f}),"
          f"  symmetric grad_x = {gs[0]:+.2e},  diffmode==0 grad_x = {gr[0]:+.2e}")
    err_diff.append(abs(gd[0] - analytic_spurious))
    err_sym.append(abs(gs[0]))
    err_dru.append(abs(gr[0]))
scale = abs(uV0 * rho0 * slope / rho0**2)
print(f"  symmetric |grad_x| / |spurious scale| stays within {max(err_sym) / scale:.2e};"
      f" diffmode==0 within {max(err_dru) / scale:.2e}; the difference form converges onto the spurious value.")
# u = u_V/rho is not linear here, so the symmetric form's residual is the
# ordinary O(h^2) curvature truncation error (it grows with h, as expected).
# diffmode==0 differences rho*u itself, so a uniform u_V gives exactly zero.
assert max(err_sym) < 0.02 * scale
assert max(err_dru) < 1e-12 * scale
assert err_diff[-1] < 0.05 * scale
print("  PASS: the difference form pumps toward the dense side; the symmetric")
print("  form leaves a uniform u_V alone to better than 1% of that spurious")
print("  signal; diffmode==0 leaves it alone exactly (it differences rho*u).")
print()

# ---------------------------------------------------------------------------
# Part 2: convergence sweep with linear u(x), rho(x)
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part 2: neighbour-count convergence sweep, linear u(x) and rho(x)")
print("=" * 78)
u0, s_u = 1.0, 0.8
u_lin = u0 * (1.0 + s_u * pos[:, 0])
TARGETS = {
    "grad(u)": lambda xp, rp: u0 * s_u,
    "(1/rho) grad(rho u)": lambda xp, rp: u0 * s_u + u0 * (1 + s_u * xp) * rho0 * slope / rp,
}
verdicts = {}
for name, fn in CONSTRUCTIONS.items():
    err = {t: [] for t in TARGETS}
    for hod in h_sweep:
        h = hod * dx_p
        g = accumulate(pos, rho, mass, u_lin, h, fn)
        e_h = {t: [] for t in TARGETS}
        for frac in (0.4, 0.5, 0.6):
            p_idx = interior_probe(pos, h, box_l, frac)
            for t, tf in TARGETS.items():
                e_h[t].append(abs(g[p_idx, 0] - tf(pos[p_idx, 0], rho[p_idx])))
        for t in TARGETS:
            err[t].append(np.mean(e_h[t]))
    print(f"\n  {name}:")
    for hod, *row in zip(h_sweep, *[err[t] for t in TARGETS]):
        nn = (4.0 / 3.0) * np.pi * hod**3
        print(f"    ~{nn:5.0f} nbrs: " + "   ".join(f"{t}={v:.3e}" for t, v in zip(TARGETS, row)))
    last = {t: err[t][-1] for t in TARGETS}
    first = {t: err[t][0] for t in TARGETS}
    best = min(last, key=last.get)
    others = [last[t] for t in TARGETS if t != best]
    if first[best] / max(last[best], 1e-300) >= 5 and last[best] * 10 < min(others):
        verdicts[name] = best
    else:
        verdicts[name] = "AMBIGUOUS"
    print(f"    ==> converges to: {verdicts[name]}")

assert verdicts["difference (original §2.2 grad(u))"] == "grad(u)"
assert verdicts["symmetric diffmode==1 (rejected)"] == "(1/rho) grad(rho u)"
assert verdicts["difference on rho*u, diffmode==0 (adopted)"] == "(1/rho) grad(rho u)"
print()
print("CONFIRMED: the difference form on u targets plain grad(u) (wrong operator")
print("for the mass-specific flux equation); both radiation_gradient_SPH")
print("diffmode==1 and diffmode==0 target (1/rho) grad(rho u), the correct one.")
print("diffmode==0 is adopted: it is minus the adjoint of the diffmode==1")
print("divergence in the m*rho inner product (verify_design_b_timestepping_")
print("stability.py Part F.3), which the staggered scheme's stability on a")
print("disordered distribution requires; diffmode==1 for grad is not, and the")
print("scheme built on it diverges there. ALL CHECKS PASSED")
