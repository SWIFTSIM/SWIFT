"""Verify the time-stepping scheme of Design B (design-lw-fuv-design-b.md
Section 4.5) by von Neumann analysis, and confirm the properties the
design document claims for it.

Three candidate per-step update schemes for the P1-relaxation system

    Du/Dt = -(1/rho) div(rho F) - u/tau + S
    DF/Dt = -(D/tau) (1/rho) grad(rho u) - F/tau,     D = c_hyp*lambda,
                                                      tau = lambda/c_hyp,

are compared on a uniform lattice (uniform rho, c_hyp, lambda, dt), where
the SPH gradient and divergence estimators of Section 2.2 both reduce to
one real, odd Fourier symbol K(k) (grad -> i*K, div -> i*K):

  S0  (the scheme originally written in Section 4.5): explicit FTCS
      transport sub-step for BOTH u and F from the old values, then the
      exact decay exp(-dt/tau) applied to both.
  S1  (S0 with swiftsim_0's exact relaxation, chemistry_flux.h's
      chemistry_part_integrate_flux_source_term, applied to both moments
      but keeping the one-loop, both-from-old-values ordering): each
      moment relaxes exactly toward its frozen-gradient target,
      u -> tau*(S - div F), F -> -D grad u.
  S2  (the scheme adopted by this revision): S1's exact relaxation, with
      the STAGGERED ordering the two existing SPHENIX loops give for free:
      u is finalized from the old F in the density ghost, grad(u) is then
      accumulated from the NEW u in the gradient loop, and F is finalized
      in the extra ghost.

Notation: a = dt/tau, e = exp(-a), X = K*lambda (dimensionless screening
length in symbol units), nu = c_hyp*K*dt (the wave Courant number in
symbol units; note X*(1-e) -> nu as a -> 0).

Part A derives, with sympy, the amplification factor of each scheme and
the exact stability condition of S2. Part B computes the actual Fourier
symbol K(k) of the Wendland C2 kernel-gradient estimator on a 3D cubic
lattice at this project's own resolution_eta, and sweeps lambda/h under
the design's own closure c_hyp = C_CFL*h/dt (Section 4.1) to show S0 and
S1 are unstable in the optically-thin regime while S2 never is. Part C
runs the schemes in the time domain on a 1D periodic chain to confirm the
symbol analysis (growth of S0 at kappa=0, boundedness of S2), then checks
the asymptotic-preserving property: S2's discrete fixed point is the
discrete screened-Poisson (Yukawa) solution with amplitude tau*S,
independent of dt and of c_hyp. Part D checks the frame-choice lag
prediction of Section 1.0.1: a source moving at v_rel through static gas
leaves a steady profile whose centroid trails the source by exactly
v_rel*tau.
"""

import numpy as np
import sympy as sp

# ---------------------------------------------------------------------------
# Part A: symbolic von Neumann analysis
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part A: symbolic amplification factors (uniform lattice, one Fourier mode)")
print("=" * 78)

e, X, nu, a = sp.symbols("e X nu a", positive=True)
K, tau, D, dt, c_hyp = sp.symbols("K tau D dt c_hyp", positive=True)
I = sp.I

# S0: explicit transport from old values, then exact decay of both.
G0 = e * sp.Matrix([[1, -I * dt * K], [-I * c_hyp**2 * dt * K, 1]])
ev0 = list(G0.eigenvals().keys())
g0_abs2 = [sp.simplify(sp.expand(ev * sp.conjugate(ev))) for ev in ev0]
target0 = e**2 * (1 + (c_hyp * K * dt) ** 2)
assert all(sp.simplify(g - target0) == 0 for g in g0_abs2), g0_abs2
print("S0 (FTCS + exact decay):        |g|^2 = e^2 * (1 + nu^2)")
print("   -> stable iff e^2 (1+nu^2) <= 1, i.e. dt/tau >= ln(1+nu^2)/2;")
print("      as tau -> inf (kappa -> 0) e -> 1 and |g| -> sqrt(1+nu^2) > 1:")
print("      unconditionally unstable in the optically-thin limit.")

# S1: exact relaxation of each moment toward its frozen-gradient target,
# both gradients read from the old values (one-loop Jacobi ordering).
G1 = sp.Matrix([[e, -I * tau * (1 - e) * K], [-I * D * (1 - e) * K, e]])
ev1 = list(G1.eigenvals().keys())
g1_abs2 = [sp.simplify(sp.expand(ev * sp.conjugate(ev))) for ev in ev1]
target1 = e**2 + tau * D * K**2 * (1 - e) ** 2  # tau*D = lambda^2
assert all(sp.simplify(g - target1) == 0 for g in g1_abs2), g1_abs2
print("S1 (exact relaxation, Jacobi):  |g|^2 = e^2 + X^2 (1-e)^2")
print("   -> stable iff X^2 <= (1+e)/(1-e) = coth(a/2);")
print("      as a -> 0 this is nu^2 <= 2a, i.e. tau <= 2 dt/nu^2: unstable")
print("      once lambda exceeds ~2 h/(C_CFL (Kh)^2), the thin regime again.")

# S2: exact relaxation, staggered ordering (F sees the NEW u).
G2 = sp.Matrix(
    [[e, -I * tau * (1 - e) * K], [-I * D * (1 - e) * e * K, e - tau * D * K**2 * (1 - e) ** 2]]
)
det2 = sp.simplify(G2.det())
tr2 = sp.simplify(G2.trace())
assert sp.simplify(det2 - e**2) == 0, det2
assert sp.simplify(tr2 - (2 * e - tau * D * K**2 * (1 - e) ** 2)) == 0, tr2
print("S2 (exact relaxation, staggered): det(G) = e^2, tr(G) = 2e - X^2 (1-e)^2")
# Roots of g^2 - T g + e^2 = 0. Complex pair when T^2 < 4 e^2: |g| = e exactly.
# Real pair: both in [-1, 1] iff |T| <= 1 + e^2. Upper bound is automatic
# (T <= 2e <= 1+e^2); the lower bound T >= -(1+e^2) gives the condition.
Xs = sp.symbols("X_s", positive=True)
T_expr = 2 * e - Xs**2 * (1 - e) ** 2
cond = sp.solve(sp.Eq(T_expr, -(1 + e**2)), Xs)
# The two roots are +-(1+e)/(1-e); only the positive one (0 < e < 1) applies.
assert any(sp.simplify(c_ - (1 + e) / (1 - e)) == 0 for c_ in cond), cond
# (1+e)/(1-e) with e = exp(-a) is coth(a/2).
assert sp.simplify((1 + sp.exp(-a)) / (1 - sp.exp(-a)) - sp.coth(a / 2).rewrite(sp.exp)) == 0
print("   -> complex roots: |g| = e exactly (decay at the physical rate, no")
print("      numerical damping and no growth); real roots stay in [-1,1] iff")
print("      X <= (1+e)/(1-e) = coth(a/2).")
print("   Limits: a -> 0 (thin): X <= 2/a  <=>  nu = c_hyp K dt <= 2")
print("           a -> inf (stiff): X <= 1  <=>  K lambda <= 1")
print("   Under the closure c_hyp = C_CFL h/dt:  a = C_CFL h/lambda, so")
print("   X <= coth(a/2) holds for ALL lambda iff C_CFL * max_k(K h) <= 2.")
print()

# ---------------------------------------------------------------------------
# Part B: the actual kernel-gradient symbol, and the closure sweep
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part B: Wendland C2 gradient symbol on a 3D cubic lattice, closure sweep")
print("=" * 78)

ETA = 1.2348          # resolution_eta in every shipped SubgridRadiation example
GAMMA_3D = 1.936492   # kernel_gamma, Wendland C2, 3D (src/kernel_hydro.h)
C_CFL = 0.1           # CFL_condition in every shipped SubgridRadiation example


def wc2_3d_dwdr(r, H):
    """dW/dr of the 3D Wendland C2 kernel with support radius H."""
    q = r / H
    norm = 21.0 / (2.0 * np.pi * H**3)
    inside = q < 1.0
    out = np.zeros_like(r)
    qi = q[inside]
    out[inside] = norm * (-4.0 * (1.0 - qi) ** 3 * (4.0 * qi + 1.0) + 4.0 * (1.0 - qi) ** 4) / H
    return out


def lattice_symbol_3d(k_hat, k_values, dx=1.0):
    """K(k) = sum_j (m/rho) dW/dr(r_j) (dx_j.k_hat / r_j) sin(k dx_j.k_hat)
    over the cubic-lattice neighbours of the origin, m/rho = dx^3.
    On a uniform lattice both the difference form and the symmetric
    (diffmode==1) form reduce to this same odd symbol.
    """
    h = ETA * dx
    H = GAMMA_3D * h
    n = int(np.ceil(H / dx)) + 1
    g = np.arange(-n, n + 1) * dx
    X_, Y_, Z_ = np.meshgrid(g, g, g, indexing="ij")
    pos = np.stack([X_.ravel(), Y_.ravel(), Z_.ravel()], axis=1)
    r = np.linalg.norm(pos, axis=1)
    keep = (r > 0) & (r < H)
    pos, r = pos[keep], r[keep]
    dwdr = wc2_3d_dwdr(r, H)
    proj = pos @ k_hat
    out = np.array([np.sum(dx**3 * dwdr * (proj / r) * np.sin(k * proj)) for k in k_values])
    return out, h


k_values = np.linspace(1e-3, np.pi, 400)
Kh_max = 0.0
for k_hat in (np.array([1.0, 0, 0]), np.array([1.0, 1.0, 0]) / np.sqrt(2), np.ones(3) / np.sqrt(3)):
    Ksym, h = lattice_symbol_3d(k_hat, k_values)
    Kh_max = max(Kh_max, np.max(np.abs(Ksym)) * h)
    print(f"  direction {np.round(k_hat, 3)}: max_k |K(k)| h = {np.max(np.abs(Ksym)) * h:.4f}")
print(f"  ==> (K h)_max = {Kh_max:.4f}; C_CFL (K h)_max = {C_CFL * Kh_max:.4f} (must be <= 2)")
assert C_CFL * Kh_max <= 2.0


def wc2_3d_w(r, H):
    q = r / H
    out = np.zeros_like(r)
    qi = q[q < 1.0]
    out[q < 1.0] = 21.0 / (2.0 * np.pi * H**3) * (1.0 - qi) ** 4 * (4.0 * qi + 1.0)
    return out


# Nyquist content of a kernel-weighted injection on the 3D lattice: the
# share of the deposit sitting in the lattice null modes (K = 0 there),
# which transport cannot move and only absorption removes (Part C.3).
h3 = ETA
H3 = GAMMA_3D * h3
n3 = int(np.ceil(H3)) + 1
g3 = np.arange(-n3, n3 + 1)
X3, Y3, Z3 = np.meshgrid(g3, g3, g3, indexing="ij")
r3 = np.sqrt(X3**2 + Y3**2 + Z3**2)
w3 = wc2_3d_w(r3, H3)
w3 /= w3.sum()
nyq_axis = abs(np.sum(w3 * (-1.0) ** X3))
nyq_corner = abs(np.sum(w3 * (-1.0) ** (X3 + Y3 + Z3)))
print(f"  Nyquist content of a kernel-weighted deposit (3D lattice, eta = {ETA}):"
      f" axis-face mode {nyq_axis:.2e}, corner mode {nyq_corner:.2e} (single-particle deposit: 1)")
print("  (the support is only 2.4 dx at this eta, so a face mode keeps a quarter of")
print("   its amplitude; on an N-particle lattice the null faces are ~3 N^(2/3) of N")
print("   modes, so their spatial weight at the source scales as N^(-1/3))")
assert nyq_axis < 0.5 and nyq_corner < 0.1


def amp_S0(e_, nu_):
    return np.sqrt(e_**2 * (1 + nu_**2))


def amp_S1(e_, X_):
    return np.sqrt(e_**2 + X_**2 * (1 - e_) ** 2)


def amp_S2(e_, X_):
    T = 2 * e_ - X_**2 * (1 - e_) ** 2
    disc = T**2 - 4 * e_**2
    g = np.where(disc < 0, e_, 0.5 * (np.abs(T) + np.sqrt(np.maximum(disc, 0.0))))
    return g


Ksym_x, h = lattice_symbol_3d(np.array([1.0, 0, 0]), k_values)
Kh = np.abs(Ksym_x) * h
lam_over_h = np.logspace(-3, 4, 71)
worst = {"S0": [], "S1": [], "S2": []}
for y in lam_over_h:
    a_ = C_CFL / y            # dt/tau under the closure
    e_ = np.exp(-a_)
    nu_ = C_CFL * Kh          # c_hyp K dt = C_CFL K h
    X_ = Kh * y               # K lambda
    worst["S0"].append(np.max(amp_S0(e_, nu_)))
    worst["S1"].append(np.max(amp_S1(e_, X_)))
    worst["S2"].append(np.max(amp_S2(e_, X_)))
for key in worst:
    worst[key] = np.array(worst[key])

print("\n  max_k |g| under the closure c_hyp = C_CFL h/dt, C_CFL = 0.1:")
print("  lambda/h      S0        S1        S2")
for y, g0, g1, g2 in zip(lam_over_h[::10], worst["S0"][::10], worst["S1"][::10], worst["S2"][::10]):
    print(f"  {y:9.3e}  {g0:8.5f}  {g1:8.5f}  {g2:8.5f}")

first_unstable_S0 = lam_over_h[np.argmax(worst["S0"] > 1 + 1e-12)]
first_unstable_S1 = lam_over_h[np.argmax(worst["S1"] > 1 + 1e-12)]
pred_S0 = 2 * C_CFL / np.log(1 + (C_CFL * Kh_max) ** 2)
print(f"\n  S0 first unstable at lambda/h ~ {first_unstable_S0:.2f} (analytic 2 C_CFL/ln(1+nu_max^2) = {pred_S0:.2f})")
print(f"  S1 first unstable at lambda/h ~ {first_unstable_S1:.2f}")
print(f"  S0 growth per step as lambda -> inf: {worst['S0'][-1]:.5f} (sqrt(1+nu_max^2) = {np.sqrt(1 + (C_CFL * Kh_max) ** 2):.5f})")
print(f"  S2 max |g| over the whole sweep: {np.max(worst['S2']):.12f}")
assert np.max(worst["S2"]) <= 1.0 + 1e-12
assert worst["S0"][-1] > 1.0 and worst["S1"][-1] > 1.0
# The margin: how much larger c_hyp could be before S2's thin-limit bound binds.
print(f"  S2 thin-limit margin: nu_max = {C_CFL * Kh_max:.3f} vs bound 2, i.e. the closure")
print(f"  coefficient could rise from C_CFL = {C_CFL} to {2 / Kh_max:.2f} before S2 destabilizes.")
print()

# ---------------------------------------------------------------------------
# Part C: time-domain check on a 1D periodic chain
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part C: time-domain check, 1D periodic chain")
print("=" * 78)

GAMMA_1D = 1.620185


def wc2_1d_dwdr(r, H):
    """dW/dr of the 1D Wendland C2 kernel, W = (5/(4H)) (1-q)^3 (1+3q)."""
    q = r / H
    inside = q < 1.0
    out = np.zeros_like(r)
    qi = q[inside]
    out[inside] = (5.0 / (4.0 * H)) * (-3.0 * (1.0 - qi) ** 2 * (1.0 + 3.0 * qi) + 3.0 * (1.0 - qi) ** 3) / H
    return out


def chain_operator(n, h, dx=1.0):
    """Antisymmetric matrix A with (grad u)_i = sum_j A_ij u_j for the
    symmetric (diffmode==1) estimator on a uniform periodic chain, m/rho =
    dx. The divergence uses the same matrix. A is exactly antisymmetric, so
    sum_i (div F)_i = 0 (Section 2.2's exact conservation).
    """
    H = GAMMA_1D * h
    A = np.zeros((n, n))
    nmax = int(np.ceil(H / dx))
    for s in range(1, nmax + 1):
        r = np.array([s * dx])
        w = wc2_1d_dwdr(r, H)[0]
        # (u_i + u_j) * dW/dr * sign(x_i - x_j) * (m/rho); the u_i part
        # cancels between +s and -s on a uniform chain, leaving the u_j part.
        for i in range(n):
            A[i, (i - s) % n] += dx * w * (+1.0)  # x_i - x_j = +s dx
            A[i, (i + s) % n] += dx * w * (-1.0)  # x_i - x_j = -s dx
    assert np.allclose(A, -A.T)
    return A


def phi(a_):
    """(1 - exp(-a))/a with the a -> 0 limit; the factor turning tau*(1-e)
    into dt*phi(a) so kappa = 0 (tau = inf) needs no special case."""
    a_ = np.asarray(a_, dtype=float)
    small = a_ < 1e-6
    out = np.empty_like(a_)
    out[small] = 1.0 - 0.5 * a_[small]
    out[~small] = -np.expm1(-a_[~small]) / a_[~small]
    return out


def step_S0(u, F, A, dt, c_hyp, a_, S):
    e_ = np.exp(-a_)
    un = u - dt * (A @ F) + dt * S
    Fn = F - dt * c_hyp**2 * (A @ u)
    return e_ * un, e_ * Fn


def step_S2(u, F, A, dt, c_hyp, a_, S):
    e_, ph = np.exp(-a_), phi(a_)
    un = e_ * u + dt * ph * (S - A @ F)
    Fn = e_ * F - dt * ph * c_hyp**2 * (A @ un)
    return un, Fn


n, dx, h, dt = 256, 1.0, 2.0, 1.0
A = chain_operator(n, h, dx)
c_hyp = C_CFL * h / dt
rng = np.random.default_rng(20260907)
u0 = rng.normal(size=n) * 1e-3
F0 = np.zeros(n)
zero = np.zeros(n)

# C.1: kappa = 0 (tau = inf, a = 0): S0 must grow, S2 must stay bounded.
uS0, FS0 = u0.copy(), F0.copy()
uS2, FS2 = u0.copy(), F0.copy()
nsteps = 5000
for _ in range(nsteps):
    uS0, FS0 = step_S0(uS0, FS0, A, dt, c_hyp, 0.0, zero)
    uS2, FS2 = step_S2(uS2, FS2, A, dt, c_hyp, 0.0, zero)
ratio_S0 = np.max(np.abs(uS0)) / np.max(np.abs(u0))
ratio_S2 = np.max(np.abs(uS2)) / np.max(np.abs(u0))
K1d = np.abs(np.linalg.eigvals(A)).max()
print(f"  C.1 kappa = 0, {nsteps} steps, c_hyp = C_CFL h/dt, nu_max = {c_hyp * K1d * dt:.3f}:")
print(f"      S0 amplitude ratio = {ratio_S0:.3e} (expected ~ (1+nu^2)^(N/2) for the worst mode)")
print(f"      S2 amplitude ratio = {ratio_S2:.3e}")
assert ratio_S0 > 1e3, "S0 did not blow up; the FTCS instability claim would be wrong"
assert ratio_S2 < 10.0, "S2 grew; the staggered-scheme stability claim would be wrong"

# C.2: asymptotic-preserving fixed point, independent of dt/tau and c_hyp.
lam = 5.0 * dx
S = np.zeros(n)
S[n // 2] = 1.0
profiles = []
for c_factor in (0.01, 0.1, 0.5):
    c_ = c_factor * h / dt
    tau_ = lam / c_
    a_ = dt / tau_
    u, F = np.zeros(n), np.zeros(n)
    for it in range(200000):
        un, Fn = step_S2(u, F, A, dt, c_, a_, S)
        if np.max(np.abs(un - u)) < 1e-13 * max(np.max(np.abs(un)), 1e-300) and it > 10:
            u, F = un, Fn
            break
        u, F = un, Fn
    # Direct solve of the discrete screened-Poisson equation u = tau S + lambda^2 A A u
    # (A A is the composed discrete Laplacian; tau*D = lambda^2).
    M = np.eye(n) - lam**2 * (A @ A)
    u_direct = np.linalg.solve(M, tau_ * S)
    rel = np.max(np.abs(u - u_direct)) / np.max(np.abs(u_direct))
    print(f"  C.2 lambda = {lam} dx, c_hyp = {c_factor} h/dt (dt/tau = {a_:.3f}, {it} iterations):")
    print(f"      |u - u_direct|/|u|_max = {rel:.2e}; sum(u)/(tau S_tot) = {np.sum(u) / (tau_ * np.sum(S)):.6f}")
    assert rel < 1e-8
    assert abs(np.sum(u) / (tau_ * np.sum(S)) - 1.0) < 1e-8
    profiles.append(u / tau_)
for p_ in profiles[1:]:
    assert np.max(np.abs(p_ - profiles[0])) / np.max(np.abs(profiles[0])) < 1e-8
x = (np.arange(n) - n // 2) * dx
yukawa_1d = np.sum(S) / (2 * lam) * np.exp(-np.abs(x) / lam)
mask = (np.abs(x) > 2 * h) & (np.abs(x) < 6 * lam)
slope_meas = np.polyfit(np.abs(x[mask]), np.log(profiles[0][mask]), 1)[0]
print(f"      u/tau profiles identical across c_hyp to 1e-8 (fixed point is dt- and c_hyp-independent);")
print(f"      e-folding length of the discrete fixed point: {-1 / slope_meas:.3f} dx (target lambda = {lam} dx)")
print(f"      amplitude at the source vs 1D Yukawa S/(2 lambda): {profiles[0][n // 2] / yukawa_1d[n // 2]:.3f}")
assert abs(-1 / slope_meas - lam) / lam < 0.05


def nyquist_content(u):
    """|u_k(pi)| / |u_k(0)|: the lattice null mode's share of the profile."""
    uk = np.fft.rfft(u)
    return np.abs(uk[-1]) / np.abs(uk[0])


# C.3: the composed Laplacian A@A has a null mode at the lattice Nyquist
# scale (K(pi/dx) = 0), so the fixed point leaves whatever Nyquist content
# the SOURCE has untouched: u_k(pi) = tau S_k(pi) exactly. A one-particle
# delta source is the worst case; the real injection is kernel-weighted.
H1 = GAMMA_1D * h
w = np.maximum(1 - np.abs(x) / H1, 0) ** 3 * (1 + 3 * np.abs(x) / H1)
S_smooth = w / np.sum(w)
tau_ = lam / (C_CFL * h / dt)
u_delta = np.linalg.solve(np.eye(n) - lam**2 * (A @ A), tau_ * S)
u_smooth = np.linalg.solve(np.eye(n) - lam**2 * (A @ A), tau_ * S_smooth)
print(f"  C.3 null-mode diagnostic (lambda = {lam} dx):")
print(f"      delta source:           peak/Yukawa = {u_delta[n // 2] / (tau_ * yukawa_1d[n // 2]):.3f},"
      f" Nyquist content = {nyquist_content(u_delta):.4f}")
print(f"      kernel-weighted source: peak/Yukawa = {u_smooth[n // 2] / (tau_ * yukawa_1d[n // 2]):.3f},"
      f" Nyquist content = {nyquist_content(u_smooth):.4f}")
assert abs(nyquist_content(u_delta) - 1.0) < 1e-10
assert nyquist_content(u_smooth) < 0.05
print()

# ---------------------------------------------------------------------------
# Part D: moving-source lag (frame-choice residual of Section 1.0.1)
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part D: moving source, steady centroid lag = -v_rel * tau")
print("=" * 78)
# Continuum prediction (either the P1 system or its diffusion limit): with a
# source moving at v through static gas, integrating the first moment of the
# u equation in the source frame gives M1 = -v tau M0 exactly, and M0 = S tau.
n, dx, h, dt = 1024, 1.0, 2.0, 1.0
A = chain_operator(n, h, dx)
lam = 8.0 * dx
c_ = C_CFL * h / dt
tau_ = lam / c_
a_ = dt / tau_
for v in (0.05, 0.1):
    u, F = np.zeros(n), np.zeros(n)
    xs = 0.0
    lags = []
    for it in range(int(60 * tau_)):
        S = np.zeros(n)
        # kernel-free linear deposit between the two bracketing particles
        i0 = int(np.floor(xs / dx))
        f = xs / dx - i0
        S[i0 % n] += (1 - f) / dx
        S[(i0 + 1) % n] += f / dx
        u, F = step_S2(u, F, A, dt, c_, a_, S)
        xs += v * dt
        if it > 40 * tau_:
            # centroid relative to the source, on the periodic chain
            rel = ((np.arange(n) * dx - xs + n * dx / 2) % (n * dx)) - n * dx / 2
            lags.append(np.sum(rel * u) / np.sum(u))
    lag = np.mean(lags)
    print(f"  v = {v} dx/dt (v/c_hyp = {v / c_:.2f}): centroid - x_source = {lag:+.3f} dx,"
          f" predicted -v tau = {-v * tau_:+.3f} dx; M0/(S tau) = {np.sum(u) * dx / tau_:.4f}")
    assert abs(lag + v * tau_) < 0.1 * v * tau_ + 0.3 * dx
print()

# ---------------------------------------------------------------------------
# Part E: 0-D transient, exact source fold-in: lag vs cumulative dose
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part E: 0-D declining source, u_{n+1} = e u_n + tau (1-e) S_n")
print("=" * 78)
# A source switched on for T_on then off. The instantaneous field lags its
# quasi-static target tau*S(t) by ~tau, but the time-integrated field (the
# cumulative dose a linear Grackle coupling would receive) equals
# tau * integral(S dt) EXACTLY for any dt/tau: a geometric-series identity
# of the exact fold-in. Nothing is lost, only delayed.
for a_ in (0.01, 1.0, 100.0):
    tau_ = 1.0
    dt_ = a_ * tau_
    e_ = np.exp(-a_)
    n_on = max(int(round(20.0 * tau_ / dt_)), 1)
    n_tot = n_on + max(int(round(40.0 * tau_ / dt_)), 5)
    S_series = np.zeros(n_tot)
    S_series[:n_on] = 1.0
    u = 0.0
    dose = 0.0
    max_lag = 0.0
    for S_n in S_series:
        u = e_ * u + tau_ * (1 - e_) * S_n
        dose += u * dt_
        max_lag = max(max_lag, abs(u - tau_ * S_n))
    print(f"  dt/tau = {a_:6.2f}: dose/(tau int S dt) = {dose / (tau_ * np.sum(S_series) * dt_):.12f},"
          f" max |u - tau S(t)|/(tau S_on) = {max_lag / tau_:.3f}")
    assert abs(dose / (tau_ * np.sum(S_series) * dt_) - 1.0) < 1e-10
print("  (the instantaneous error is the lag; the dose is exact; stiff dt/tau")
print("   has no lag at all, u jumps to tau S each step)")
print()

# ---------------------------------------------------------------------------
# Part F: cross-checks against Design A's three failure modes
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part F: cross-checks against Design A's failure modes (HANDOFF brief)")
print("=" * 78)


def chain_operators_general(xpos, m, rho, hh, box):
    """Per-particle m, rho, h on a periodic chain. Returns
    A_sym : radiation_divergence_SPH diffmode==1 (symmetric),
            (A q)_i = sum_j m_j s_ij (q_i w(r,h_i)/rho_i + q_j w(r,h_j)/rho_j);
    B_dif : radiation_gradient_SPH diffmode==0 (difference on rho*q),
            (B q)_i = (1/rho_i^2) sum_j m_j (rho_j q_j - rho_i q_i) s_ij w(r,h_i).
    B_dif is minus the adjoint of A_sym in the inner product with weight
    m*rho (the discrete volumetric energy), which is what makes the
    staggered scheme energy-conserving on a disordered distribution."""
    n_ = len(xpos)
    A_ = np.zeros((n_, n_))
    B_ = np.zeros((n_, n_))
    for i in range(n_):
        for j in range(n_):
            if i == j:
                continue
            d = xpos[i] - xpos[j]
            d -= box * np.round(d / box)
            r = abs(d)
            s = np.sign(d)
            wi = wc2_1d_dwdr(np.array([r]), GAMMA_1D * hh[i])[0]
            wj = wc2_1d_dwdr(np.array([r]), GAMMA_1D * hh[j])[0]
            if wi == 0.0 and wj == 0.0:
                continue
            A_[i, i] += m[j] * s * wi / rho[i]
            A_[i, j] += m[j] * s * wj / rho[j]
            B_[i, i] += -m[j] * s * wi / rho[i]
            B_[i, j] += m[j] * rho[j] * s * wi / rho[i] ** 2
    return A_, B_


# F.3 (failure 3, injection conservation, and the disordered-distribution
# stability the uniform-lattice analysis of Parts A-C cannot see): the
# FULL scheme, kappa = 0, jittered positions, non-uniform m, rho, h,
# per-particle c_hyp. Two pairings are compared: (div, grad) = (sym, sym),
# which a first draft of this revision adopted, and (sym, dif).
n_, box = 96, 96.0
xpos = np.arange(n_) * 1.0 + rng.uniform(-0.2, 0.2, size=n_)
m = rng.uniform(0.5, 2.0, size=n_)
hh = rng.uniform(1.6, 2.6, size=n_)
rho = m / 1.0
Ag, Bg = chain_operators_general(xpos, m, rho, hh, box)
W = np.diag(m * rho)
assert np.allclose(m @ Ag, 0.0, atol=1e-12)                 # sum_i m_i (div F)_i = 0
assert np.allclose(W @ Bg, -(W @ Ag).T, atol=1e-12)         # B = -W^-1 A^T W
assert not np.allclose(W @ Ag, -(W @ Ag).T, atol=1e-6)      # A is NOT skew in W
print("  F.3 disordered chain, non-uniform m/rho/h: W B_dif = -(W A_sym)^T exactly;")
print("      A_sym is not skew-adjoint to itself.")
S = np.zeros(n_)
S[n_ // 3] = 0.7
S[2 * n_ // 3] = 0.3
c_i = C_CFL * hh / dt
nsteps = 2000
for label, Bop in (("(sym, sym) pairing", Ag), ("(sym, dif) pairing", Bg)):
    u, F = np.zeros(n_), np.zeros(n_)
    energy = []
    for it in range(nsteps):
        un = u + dt * (S - Ag @ F)
        Fn = F - dt * c_i**2 * (Bop @ un)
        u, F = un, Fn
        if it == nsteps // 2:
            S_off = S.copy()
            S = np.zeros(n_)  # switch the source off, then watch the energy
        if it > nsteps // 2:
            energy.append(np.sum(m * rho * (u**2 + F**2 / c_i**2)))
    S = S_off
    budget = np.sum(m * u) / (np.sum(m * S) * dt * (nsteps // 2 + 1))
    energy = np.array(energy)
    excursion = (energy.max() - energy.min()) / energy[0]
    trend = np.polyfit(np.arange(len(energy)), energy / energy[0], 1)[0] * len(energy)
    print(f"      {label}: sum(m u)/injected = {budget:.12f}, max|u| = {np.max(np.abs(u)):.3e},"
          f" energy after source-off: excursion {excursion:.2e}, trend {trend:+.2e}")
    if label.startswith("(sym, sym)"):
        assert np.max(np.abs(u)) > 1e3, "expected the (sym, sym) pairing to blow up"
    else:
        assert abs(budget - 1.0) < 1e-10
        # symplectic Euler conserves a modified energy: bounded O(nu) oscillation, no trend
        assert excursion < 0.1 and abs(trend) < 5e-3
        assert np.max(np.abs(u)) < 1e3
print("      (sym, sym) grows ~10% per step; (sym, dif) conserves sum(m u) to round-off")
print("      and keeps the discrete energy sum m rho (u^2 + F^2/c^2) within a bounded")
print("      O(nu) oscillation with no secular trend (source off). max|u| for (sym, dif)")
print("      is the Nyquist content of the two single-particle sources accumulating on")
print("      those particles with kappa = 0 (Part C.3): confined, linear in time, small")
print("      for kernel-weighted injection.")

# F.2 (failure 2, kappa -> 0): bounded free streaming with continuous
# injection: a 1D source in vacuum must give the plateau u = S/(2 c_hyp)
# behind the front, and the front must not outrun c_hyp t.
n_, dx, h, dt = 2048, 1.0, 2.0, 1.0
A = chain_operator(n_, h, dx)
c_ = C_CFL * h / dt
S = np.zeros(n_)
S[n_ // 2] = 1.0
x = (np.arange(n_) - n_ // 2) * dx
# group speed of the discrete wave: c_hyp * dK/dk, never above c_hyp
ks = 2 * np.pi * np.arange(1, n_ // 2) / n_
Ks = np.array([np.imag((A @ np.exp(1j * kk * np.arange(n_)))[0]) for kk in ks])
print(f"  F.2 kappa = 0: max group speed / c_hyp = {np.max(np.gradient(Ks, ks)):.3f} (1D chain symbol)")
u, F = np.zeros(n_), np.zeros(n_)
for T in (1000, 3000):
    while_steps = T - (1000 if T == 3000 else 0)
    for _ in range(while_steps):
        u, F = step_S2(u, F, A, dt, c_, 0.0, S)
    # median: robust to the dispersive ringing behind the front at early times
    plateau = np.median(u[(np.abs(x) > 20) & (np.abs(x) < 0.5 * c_ * T * dt)])
    fronts = {thr: np.max(np.abs(x[u > thr * plateau])) for thr in (1e-1, 1e-2, 1e-3)}
    print(f"      {T} steps: plateau u = {plateau:.4f} (S/(2 c_hyp) = {1 / (2 * c_):.4f}), c_hyp t = {c_ * T * dt:.0f} dx;"
          f" front at 10%/1%/0.1% of plateau: {fronts[1e-1]:.0f}/{fronts[1e-2]:.0f}/{fronts[1e-3]:.0f} dx;"
          f" sum(u) dx/(S t) = {np.sum(u) * dx / T:.6f}")
    assert abs(plateau * 2 * c_ - 1.0) < 0.05
    assert fronts[1e-1] <= c_ * T * dt + 2 * h
    assert fronts[1e-2] <= c_ * T * dt + 8 * h
    assert abs(np.sum(u) * dx / T - 1.0) < 1e-10
print("      bounded free streaming: correct 1D plateau, energy = S t exactly, the")
print("      10%-level front at c_hyp t; the low-level precursor ahead of it is the")
print("      dispersive tail of a scheme with no numerical damping (its width grows")
print("      sub-linearly with t), not a superluminal mode.")

# F.1 (failure 1, regime of validity): the discrete fixed point's e-folding
# length against the chosen lambda, across lambda/h from the production
# clump regime (lambda << h) to the thin regime. lambda is a chosen input
# here (D = c_hyp lambda), not an emergent one, so it cannot collapse the
# way Design A's did; what CAN happen is that the discrete screened-Poisson
# solution cannot represent a decay shorter than the particle spacing.
n_, dx, h, dt = 1024, 1.0, 2.0, 1.0
A = chain_operator(n_, h, dx)
S = np.zeros(n_)
S[n_ // 2] = 1.0
x = (np.arange(n_) - n_ // 2) * dx
print("  F.1 discrete fixed point vs chosen lambda (h = 2 dx, kernel support 3.2 dx):")
print("      lambda/h   lambda_eff/lambda (far field)   sum(u)/(tau S)")
for lam in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0):
    lam_dx = lam * h
    c_ = C_CFL * h / dt
    tau_ = lam_dx / c_
    u_fp = np.linalg.solve(np.eye(n_) - lam_dx**2 * (A @ A), tau_ * S)
    # e-folding length from the far-field envelope (skip the source region
    # and the null-mode odd-even pattern: use even-indexed particles only)
    sel = (np.abs(x) > 3 * h) & (np.abs(x) < 3 * h + 12 * max(lam_dx, 2 * dx)) & (np.arange(n_) % 2 == 0)
    pos_u = u_fp[sel] > 0
    if np.sum(pos_u) >= 4:
        slope_ = np.polyfit(np.abs(x[sel][pos_u]), np.log(u_fp[sel][pos_u]), 1)[0]
        lam_eff = -1.0 / slope_
    else:
        lam_eff = np.nan
    tot = np.sum(u_fp) / (tau_ * np.sum(S))
    print(f"      {lam:8.2f}   {lam_eff / lam_dx:10.3f}                {tot:.10f}")
    assert abs(tot - 1.0) < 1e-9
    if lam >= 2.0:
        assert abs(lam_eff / lam_dx - 1.0) < 0.05
print("      total energy tau*S exact at every lambda; the e-folding length is")
print("      within 5% of the chosen lambda for lambda >= 2 h and is bounded")
print("      below by the particle spacing for lambda < h (unresolved, not")
print("      collapsed: a clump-regime field stays on the source's own kernel).")

# F.4 heterogeneous coefficients: per-particle tau and c_hyp varying by
# orders of magnitude across the chain (a metallicity gradient plus a
# time-bin spread), continuous injection, must stay bounded. The uniform
# von Neumann analysis does not cover this; this is a smoke test, not a proof.
n_, dx, dt = 512, 1.0, 1.0
xpos = np.arange(n_) * dx + rng.uniform(-0.2, 0.2, size=n_)
hh = rng.uniform(1.6, 2.6, size=n_)
m = rng.uniform(0.5, 2.0, size=n_)
rho = m / dx
A, B = chain_operators_general(xpos, m, rho, hh, float(n_))
S = np.zeros(n_)
S[n_ // 2] = 1.0
kappa_i = 10.0 ** rng.uniform(-3, 0.5, size=n_)  # lambda_i from 0.3 dx to 1000 dx
dt_i = dt * 2.0 ** rng.integers(0, 4, size=n_)   # a 3-bin time-step spread
c_i = C_CFL * hh / dt_i
a_i = c_i * kappa_i * dt
u, F = np.zeros(n_), np.zeros(n_)
peak = []
for it in range(20000):
    e_, ph = np.exp(-a_i), phi(a_i)
    un = e_ * u + dt * ph * (S - A @ F)
    Fn = e_ * F - dt * ph * c_i**2 * (B @ un)
    u, F = un, Fn
    if it % 1000 == 999:
        peak.append(np.max(np.abs(u)))
print(f"  F.4 heterogeneous tau/c_hyp/m/rho/h (lambda 0.3-1000 dx, 3 time bins), 20000 steps:")
print(f"      max|u| every 1000 steps: " + " ".join(f"{p_:.3g}" for p_ in peak))
assert peak[-1] < 2.0 * max(peak[:5]) and np.all(np.isfinite(u))
print("      bounded (no growth over the last 15000 steps); a smoke test only.")
print()
print("ALL CHECKS PASSED")
