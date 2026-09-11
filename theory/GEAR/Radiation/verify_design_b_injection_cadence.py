"""Verify the injection-cadence analysis of the Design B dissipation
design doc (design-lw-fuv-design-b-dissipation.md Section 4.6): the
current once-per-STAR-step lump-sum injection versus a per-gas-sub-step
source RATE applied inside the exact-relaxation `u` update.

Both schemes use the same exact-relaxation integrator (radiation_isrf.c,
`u_new = e*u_prev + dt*phi*(source - div_F)`, `e = exp(-a)`,
`phi = (1-e)/a`). They differ only in WHEN the star's emission enters `u`:

  L (lump, current code, radiation_iact.h:259-305): the star's whole step
    `Delta_t_star*S*phi(Delta_t_star/tau)` is added ONCE, after the gas
    density ghost, gradient loop, extra ghost and cooling of the step the
    star is active in. A gas particle on a 2^k-finer bin runs its other
    sub-steps with no source at all.
  R (rate, proposed): the star writes the mass-specific rate `S` on the
    gas particle; every gas sub-step adds `dt_gas*phi(dt_gas/tau)*S`
    inside the density-ghost update, exactly where `div_F` already enters.

Part A: 0-D (no transport). Closed-form steady states, time averages,
        the value cooling reads, and the per-star-step budget of each
        scheme.
Part B: 1-D periodic chain with transport (the machinery of
        verify_design_b_timestepping_stability.py Part C), kernel-weighted
        source, stiff and thin regimes, gas sub-steps N per star step.
        B.1 measures what cooling reads under each scheme and the L(N=1)
        vs R fixed-point difference; B.2 the kappa = 0 budget; B.3 the
        3-D lattice size of the transport term acting on L's between-lump
        remnant at the Tier-1 corner; B.4 the thin-regime N sweep, where L
        alone (no density contrast) produces the undershoot/overshoot
        class the 1627 log's pinned runs showed.
Part C: 0-D declining source with a star death: dose error of L (dose for
        the elapsed star step, deposited at its end) and R (rate held
        forward, expiring RADIATION_LW_FUV_TAG_LIFETIME_INTERVALS star steps
        after the last touch) against the exact continuous solution.
Part D: several stars on different time bins: a summed rate with
        reset-on-first-touch loses energy; the dose-reservoir form of
        Section 4.6.5 applies every deposited dose exactly once and reduces
        to the constant rate R for one star.

Exit 0 iff every assertion holds.
"""
# =============================================================================
# M1 CLOSURE AUDIT, 2026-09-11: CHECKED, CLOSURE-INDEPENDENT, NO CHANGE.
# Both schemes compared here differ only in WHEN the source enters the
# exact-relaxation `u` update, which the closure does not touch. The
# constant multiplying the source did change (`3*c_hyp/c` -> `c_hyp/c`,
# design-lw-fuv-m1-upgrade.md D2), but it multiplies both schemes
# identically and cancels out of every comparison made below.
# =============================================================================
import numpy as np

# ---------------------------------------------------------------------------
# Shared pieces (same conventions as verify_design_b_timestepping_stability.py)
# ---------------------------------------------------------------------------
GAMMA_1D = 1.620185
C_HYP = 0.5  # LW_FUV_c_hyp_margin default
TAG_LIFETIME_INTERVALS = 2  # RADIATION_LW_FUV_TAG_LIFETIME_INTERVALS


def phi(a_):
    """(1 - exp(-a))/a with the a -> 0 limit, radiation_relaxation_phi_factor."""
    a_ = np.asarray(a_, dtype=float)
    small = a_ < 1e-6
    out = np.empty_like(a_)
    out[small] = 1.0 - 0.5 * a_[small] + a_[small] ** 2 / 6.0
    out[~small] = -np.expm1(-a_[~small]) / a_[~small]
    return out


def wc2_1d_w(r, H):
    """1D Wendland C2 kernel, W = (5/(4H)) (1-q)^3 (1+3q)."""
    q = np.abs(r) / H
    return np.where(q < 1.0, (5.0 / (4.0 * H)) * (1 - q) ** 3 * (1 + 3 * q), 0.0)


def wc2_1d_dwdr(r, H):
    """dW/dr of the 1D Wendland C2 kernel."""
    q = r / H
    inside = q < 1.0
    out = np.zeros_like(r)
    qi = q[inside]
    out[inside] = (
        (5.0 / (4.0 * H))
        * (-3.0 * (1.0 - qi) ** 2 * (1.0 + 3.0 * qi) + 3.0 * (1.0 - qi) ** 3)
        / H
    )
    return out


def chain_operator(n, h, dx=1.0):
    """Antisymmetric kernel-gradient matrix on a uniform periodic chain
    (verify_design_b_timestepping_stability.py Part C)."""
    H = GAMMA_1D * h
    A = np.zeros((n, n))
    nmax = int(np.ceil(H / dx))
    for s in range(1, nmax + 1):
        w = wc2_1d_dwdr(np.array([s * dx]), H)[0]
        for i in range(n):
            A[i, (i - s) % n] += dx * w * (+1.0)
            A[i, (i + s) % n] += dx * w * (-1.0)
    assert np.allclose(A, -A.T)
    return A


# ---------------------------------------------------------------------------
# Part A: 0-D, no transport
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part A: 0-D, no transport; star step Delta = N * gas sub-step delta")
print("=" * 78)
print("  a = delta/tau (gas sub-step in units of tau); S = 1, tau = 1.")
print("  Columns: post-injection steady state (L), cooling-read value at each")
print("  gas sub-step (L: the pre-injection value; R: the only value there is),")
print("  continuous time-average over the star step, and deposit per star step.")
tau, S = 1.0, 1.0
for a in (0.01, 1.0, 3.0, 20.0):
    for N in (1, 2, 4, 16):
        delta = a * tau
        Delta = N * delta
        e_d = np.exp(-delta / tau)
        e_D = np.exp(-Delta / tau)
        # L: iterate to the periodic steady state.
        u = 0.0
        for _ in range(int(60.0 * tau / Delta) + 50):
            for k in range(N):
                u = e_d * u  # decay-only sub-step
                if k == N - 1:
                    u_pre = u
                    u += S * tau * (1.0 - e_D)  # lump, after cooling read
        u_post_L = u
        # Cooling-read values over one period under L: sub-steps 1..N see
        # e_d^k * u_post (the last one being u_pre).
        reads_L = np.array([e_d ** (k + 1) * u_post_L for k in range(N)])
        # R: every sub-step u = e_d u + S tau (1 - e_d); fixed point.
        u = 0.0
        for _ in range(int(60.0 * tau / delta) + 50):
            u = e_d * u + S * tau * (1.0 - e_d)
        u_R = u
        # Continuous-time averages over the star step.
        avg_L = u_post_L * (tau / Delta) * (1.0 - e_D)  # = S tau phi(Delta/tau)
        avg_R = u_R
        deposit_L = S * tau * (1.0 - e_D)  # = S Delta phi(Delta/tau)
        deposit_R = N * S * tau * (1.0 - e_d)  # = S Delta phi(delta/tau)
        # Emitted per star step: S*Delta for both. The exact-relaxation deposit
        # `dt*phi*S` is emission minus the absorption of the fresh emission
        # within that dt (folded into phi, never seen as a field); the rest is
        # absorbed through the gas's own per-sub-step decay, integral(u/tau dt)
        # = avg*Delta/tau, which IS seen as a field. Both sum to S*Delta.
        hidden_L = S * Delta * (1.0 - float(phi(Delta / tau)))
        hidden_R = S * Delta * (1.0 - float(phi(delta / tau)))
        # Discrete field-decay absorption over one period: sum_k u_k (1 - e_d)
        # with u_k the value entering sub-step k.
        decay_L = sum(e_d**k * u_post_L * (1.0 - e_d) for k in range(N))
        decay_R = N * u_R * (1.0 - e_d)
        assert abs(u_post_L / (S * tau) - 1.0) < 1e-10
        assert abs(u_R / (S * tau) - 1.0) < 1e-10
        assert abs(avg_L / (S * tau) - float(phi(Delta / tau))) < 1e-10
        assert abs(deposit_L - S * Delta * float(phi(Delta / tau))) < 1e-12
        assert abs(deposit_R - S * Delta * float(phi(delta / tau))) < 1e-12
        assert abs(decay_L - deposit_L) < 1e-10 and abs(decay_R - deposit_R) < 1e-10
        assert abs(decay_L + hidden_L - S * Delta) < 1e-10
        assert abs(decay_R + hidden_R - S * Delta) < 1e-10
        assert abs(avg_L * Delta / tau - deposit_L) < 1e-10
        print(
            f"  a={a:5.2f} N={N:2d}: L post/(S tau)={u_post_L / (S * tau):.6f}"
            f"  L cooling-read min/(S tau)={reads_L.min() / (S * tau):.4f}"
            f"  L avg/(S tau)={avg_L / (S * tau):.4f}"
            f"  R/(S tau)={u_R / (S * tau):.6f}"
            f"  emission hidden in phi: L {hidden_L / (S * Delta):.4f}, R {hidden_R / (S * Delta):.4f}"
        )
print("  => Both schemes conserve emitted = hidden + field-absorbed, but L hides")
print("     1 - phi(Delta/tau) of the emission inside the injection factor (no")
print("     consumer ever sees it as a field), R only 1 - phi(delta/tau). L's")
print("     post-injection value is exact but is seen only by the snapshot; the")
print("     value cooling reads (pre-injection) is exp(-k delta/tau) of it, and L's")
print("     time-average is S tau phi(Delta/tau) for EVERY N, including N = 1.")
print("     R is the exact solution at every sub-step.")
print()

# ---------------------------------------------------------------------------
# Part B: 1-D chain with transport
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part B: 1-D periodic chain with transport, kernel-weighted source")
print("=" * 78)
n, dx, h = 256, 1.0, 2.0
A = chain_operator(n, h, dx)
x = (np.arange(n) - n // 2) * dx
H1 = GAMMA_1D * h
w_src = wc2_1d_w(x, H1)
w_src /= np.sum(w_src) * dx
S_vec = w_src.copy()  # mass-specific rate per particle, m/rho = dx
delta = 1.0
c_hyp = C_HYP * h / delta


def run_L(lam, N, n_periods):
    """Lump scheme: N decay+transport sub-steps, then the star's lump."""
    tau_ = lam / c_hyp
    a_ = delta / tau_
    e_, ph = np.exp(-a_), float(phi(a_))
    lump = S_vec * (N * delta) * float(phi(N * a_))
    u, F = np.zeros(n), np.zeros(n)
    reads = None
    for p in range(n_periods):
        reads = []
        for k in range(N):
            u_pre = e_ * u - delta * ph * (A @ F)  # density ghost
            F = e_ * F - delta * ph * c_hyp**2 * (A @ u_pre)  # extra ghost from u_pre
            reads.append(u_pre.copy())  # what cooling reads this sub-step
            u = u_pre
            if k == N - 1:
                u = u + lump  # injection, after cooling
    return u, F, np.array(reads)


def run_R(lam, N, n_periods):
    """Rate scheme: every sub-step adds delta*phi*S inside the u update."""
    tau_ = lam / c_hyp
    a_ = delta / tau_
    e_, ph = np.exp(-a_), float(phi(a_))
    u, F = np.zeros(n), np.zeros(n)
    reads = None
    for p in range(n_periods):
        reads = []
        for k in range(N):
            u = e_ * u + delta * ph * (S_vec - A @ F)
            F = e_ * F - delta * ph * c_hyp**2 * (A @ u)
            reads.append(u.copy())
    return u, F, np.array(reads)


def efold(u):
    mask = (np.abs(x) > 1.5 * h) & (u > 1e-12 * u.max())
    if mask.sum() < 4:
        return np.nan
    return -1.0 / np.polyfit(np.abs(x[mask]), np.log(u[mask]), 1)[0]


for h_over_lam in (6.0, 0.6):
    lam = h / h_over_lam
    a_ = C_HYP * h_over_lam
    print(
        f"  h/lambda = {h_over_lam} (a = C_hyp h/lambda = {a_:.2f} per gas sub-step):"
    )
    n_periods = int(max(200, 50 / a_))
    uR, FR, readsR = run_R(lam, 1, 4 * n_periods)
    peak_R = uR[n // 2]
    for N in (1, 2, 4, 8):
        uL, FL, readsL = run_L(lam, N, n_periods)
        uR_N, _, readsR_N = run_R(lam, N, n_periods)
        peak_L = uL[n // 2]
        # Source-kernel window: where the lump lands.
        ker = w_src > 0
        min_read_L = readsL[:, ker].min()
        argmin = np.unravel_index(readsL[:, ker].argmin(), readsL[:, ker].shape)
        min_read_R = readsR_N.min()
        centre_reads_L = readsL[:, n // 2]
        print(
            f"    N={N}: L post-injection peak / R peak = {peak_L / peak_R:.4f};"
            f" L cooling-read min in kernel = {min_read_L / peak_L:+.4f} x L peak"
            f" (sub-step {argmin[0] + 1}/{N});"
            f" L centre cooling-read mean = {centre_reads_L.mean() / peak_L:.4f} x L peak;"
            f" R min anywhere = {min_read_R / peak_R:+.2e} x R peak"
        )
        # R is the same fixed point at every sub-step and non-negative.
        assert np.max(np.abs(uR_N - uR)) / peak_R < 1e-8
        assert min_read_R > -1e-9 * peak_R
        if h_over_lam == 6.0:
            # Stiff: what cooling reads under L is the e-suppressed trough of
            # the sawtooth (Part A), POSITIVE at the centre: F was relaxed from
            # the pre-injection profile, so transport cannot drive the centre
            # negative (the plan-review's round-1 "-40%" and round-2
            # "e(e-0.44) < 0" estimates both used the continuum Laplacian of
            # the kernel, 20x the discrete composed operator's; see B.3).
            e_stiff = np.exp(-a_)
            assert -1e-3 * peak_L < readsL[-1, n // 2] < 2.0 * e_stiff * peak_L
            if N == 1:
                assert readsL[-1, n // 2] > 0.0
            sawtooth_mean = e_stiff * (1 - e_stiff**N) / (N * (1 - e_stiff))
            assert abs(centre_reads_L.mean() / peak_L - sawtooth_mean) < 0.01
    print(
        f"    fixed-point e-folding: L(N=1) post-injection {efold(uL_1 := run_L(lam, 1, n_periods)[0]):.3f} dx,"
        f" R {efold(uR):.3f} dx (lambda = {lam:.3f} dx);"
        f" sum(u) dx / (tau S_tot): L {np.sum(uL_1) * dx / (lam / c_hyp):.5f}, R {np.sum(uR) * dx / (lam / c_hyp):.5f}"
    )
print("  => R's fixed point is what isrf_hyperbolic_propagation_check.py's")
print("     discrete_steady_state_lambda() already iterates (injection inside the")
print("     u update, gradient from the post-injection u); L(N=1) differs from it.")
print()

# B.2: exact conservation under both schemes with kappa = 0 (no sink).
print("  B.2 kappa = 0 budget: sum(u) dx after 8 star steps of N = 4 sub-steps")
lam_inf = 1e12
for name, fn in (("L", run_L), ("R", run_R)):
    u_, _, _ = fn(lam_inf, 4, 8)
    total = np.sum(u_) * dx
    injected = 8 * 4 * delta * np.sum(S_vec) * dx  # phi -> 1
    print(f"      {name}: sum(u) dx / injected = {total / injected:.12f}")
    assert abs(total / injected - 1.0) < 1e-10
print()

# B.3: the size of the transport term acting on the decayed remnant between
# two lumps, in 3-D, on the cubic lattice at eta = 1.2348 with the real
# Wendland C2 and the symmetric pair gradient. Sub-step 2 after a lump
# holds e*(e*I + lam^2 lap_h(I)) at the centre, so the sign there is that of
# e + lam^2 lap_h(w)/w. The continuum kernel Laplacian lap W(0) = -60/H^2
# gives -16 (lam/h)^2 = -0.44 at h/lam = 6; the discrete COMPOSED operator
# div_h(grad_h(.)) of a kernel-scale bump is far weaker.
print("  B.3 3-D cubic lattice, eta = 1.2348: lam^2 lap_h(w)/w at the source centre")
import scipy.sparse as sps

ETA, GAMMA_3D = 1.2348, 1.936492
h3 = ETA
H3 = GAMMA_3D * h3


def wc2_3d_dwdr(r):
    q = r / H3
    return (21 / (2 * np.pi)) / H3**4 * (-20 * q * (1 - q) ** 3) if q < 1 else 0.0


def wc2_3d_w(r):
    q = r / H3
    return (21 / (2 * np.pi)) / H3**3 * (1 - q) ** 4 * (1 + 4 * q) if q < 1 else 0.0


Lbox = 15
g = np.arange(Lbox) - Lbox // 2
X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
pos = np.stack([X, Y, Z], -1).reshape(-1, 3)
n3 = len(pos)
offs = [np.array(o) - 4 for o in np.ndindex(9, 9, 9)]
offs = [o for o in offs if 0 < np.linalg.norm(o) < H3]
rows, cols, vals = [], [], [[], [], []]
ii3 = np.arange(n3)
for o in offs:
    r = np.linalg.norm(o)
    d = wc2_3d_dwdr(r)
    pj = (pos + o + Lbox // 2) % Lbox
    j = np.ravel_multi_index((pj[:, 0], pj[:, 1], pj[:, 2]), (Lbox, Lbox, Lbox))
    rows.append(ii3)
    cols.append(j)
    for c in range(3):
        vals[c].append(np.full(n3, d * o[c] / r))
rows, cols = np.concatenate(rows), np.concatenate(cols)
G3 = [
    sps.csr_matrix((np.concatenate(vals[c]), (rows, cols)), shape=(n3, n3))
    for c in range(3)
]
assert max(abs(G3[c] + G3[c].T).max() for c in range(3)) == 0.0
r0 = np.linalg.norm(pos, axis=1)
wb = np.array([wc2_3d_w(r) for r in r0])
wb /= wb.sum()
c0 = np.argmin(r0)
lap3 = sum(G3[c] @ (G3[c] @ wb) for c in range(3))
ker = wb > 0
for hl in (6.0, 3.0, 1.0):
    lam = h3 / hl
    e_ = np.exp(-C_HYP * hl)
    ratio = lam**2 * lap3[c0] / wb[c0]
    edge = (e_ * wb[ker] + lam**2 * lap3[ker]) / wb[c0]
    print(
        f"      h/lam = {hl}: discrete {ratio:+.4f} vs continuum -60 (lam/H)^2 = {-60 * (lam / H3) ** 2:+.4f};"
        f" e + discrete = {e_ + ratio:+.4f} at the centre;"
        f" min over the kernel of (e w + lam^2 lap_h w)/w_peak = {edge.min():+.1e} at r = {r0[ker][edge.argmin()] / h3:.2f} h"
    )
    if hl == 6.0:
        assert ratio > -0.05 and e_ + ratio > 0.0
        assert edge.min() < 0.0 and edge.min() > -1e-3
print("  => At the Tier-1 corner (h/lam = 6) the centre stays positive at every")
print("     sub-step under L; only the outer kernel edge undershoots, at the")
print("     1e-4-of-peak level, where the local scale is equally suppressed.")
print()

# B.4: the THIN regime, where the causal-reach leg's failure sat (h/lambda ~
# 0.3, 1627 log), swept in sub-steps per star step up to the bin splits the
# pinned runs actually had (~8 at 3e4 K, ~32-64 at 1e6 K). Uniform density,
# uniform c_hyp: the only difference from R is when the emission enters u.
print("  B.4 thin regime, N sub-steps per star step (uniform rho, uniform c_hyp)")
n, h = 512, 2.0
A = chain_operator(n, h, dx)
x = (np.arange(n) - n // 2) * dx
w_src = wc2_1d_w(x, GAMMA_1D * h)
w_src /= np.sum(w_src) * dx
S_vec = w_src.copy()
for h_over_lam in (0.6, 0.3):
    lam = h / h_over_lam
    a_ = C_HYP * h_over_lam
    uR, _, _ = run_R(lam, 1, int(400 / a_))
    peak_R = uR[n // 2]
    for N in (1, 2, 4, 8, 32):
        periods = max(int(400 / (a_ * N)), 40)
        uL, _, readsL = run_L(lam, N, periods)
        min_L = readsL.min() / uL[n // 2]
        print(
            f"      h/lam = {h_over_lam}, N = {N:2d}: L post-injection peak / R peak = {uL[n // 2] / peak_R:.3f};"
            f" min u over the star step / L peak = {min_L:+.4f}"
        )
        if N == 1:
            assert min_L > -1e-12
        if N >= 8:
            assert min_L < -0.05
    assert uR.min() > -1e-12 * peak_R
print("  => With N >= 8 sub-steps per star step the lump scheme alone, with no")
print("     density contrast, undershoots by 6-17% of its own peak and overshoots")
print("     R's peak 2-4x; R has neither. Saturates once N a >~ 2 (the lump")
print("     saturates at S tau). Same regime and magnitude class as the 1627")
print("     log's pinned runs (-26% of plateau, +213% bump at the pinned particle).")
print()

# ---------------------------------------------------------------------------
# Part C: declining source and star death, 0-D, cumulative dose
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part C: 0-D, L(t) declining, star dies; dose error of L and R vs exact")
print("=" * 78)
# Star: Delta_star = 0.1 Myr steps; L(t) = 1 for t < 3 Myr, (t/3)^-1 after,
# 0 from t_death on. Gas: N sub-steps per star step. tau in Myr.
Delta_star = 0.1
t_death = 8.0
t_end = t_death + 1.0


def L_of_t(t):
    if t >= t_death:
        return 0.0
    return 1.0 if t < 3.0 else (t / 3.0) ** -1.0


for tau_ in (0.005, 0.05, 0.5):
    for N in (1, 4):
        delta = Delta_star / N
        e_d = np.exp(-delta / tau_)
        n_star = int(round(t_end / Delta_star))
        # Exact continuous solution, integrated on a fine grid within each sub-step.
        fine = 200
        u_ex = 0.0
        dose_ex = 0.0
        # L scheme
        u_L = 0.0
        dose_L = 0.0
        # R scheme: rate held forward, expiry TAG_LIFETIME_INTERVALS star steps
        # after the last touch by a star with nonzero L.
        u_R = 0.0
        dose_R = 0.0
        rate_R = 0.0
        end_R = -1.0
        for s in range(n_star):
            t0 = s * Delta_star
            # Star active at t0 + Delta_star (end of its step), with L evaluated
            # at the step it just completed (the code's own Delta_t lookback).
            for k in range(N):
                t_sub0 = t0 + k * delta
                # exact
                for f in range(fine):
                    tt = t_sub0 + (f + 0.5) * delta / fine
                    dtt = delta / fine
                    u_ex = u_ex * np.exp(-dtt / tau_) + L_of_t(tt) * tau_ * (
                        1 - np.exp(-dtt / tau_)
                    )
                    dose_ex += u_ex * dtt
                # L: decay-only sub-step; cooling reads it
                u_L = e_d * u_L
                dose_L += u_L * delta
                # R: expiry check at the gas's own step start, then apply rate
                if t_sub0 >= end_R:
                    rate_R = 0.0
                u_R = e_d * u_R + rate_R * tau_ * (1 - e_d)
                dose_R += u_R * delta
            # Star's own step ends here: injection.
            t1 = t0 + Delta_star
            L_step = L_of_t(t0)  # L over the elapsed step, piecewise constant
            if L_step > 0.0:
                u_L += L_step * tau_ * (1 - np.exp(-Delta_star / tau_))
                rate_R = L_step
                end_R = t1 + TAG_LIFETIME_INTERVALS * Delta_star
        dose_true = tau_ * (3.0 + 3.0 * np.log(t_death / 3.0))  # tau * int L dt
        print(
            f"  tau={tau_:5.3f} Myr N={N}: dose_exact/(tau int L) = {dose_ex / dose_true:.5f};"
            f" L: {dose_L / dose_ex:.4f} x exact; R: {dose_R / dose_ex:.4f} x exact"
        )
        # L undercounts by ~phi(Delta/tau) (stiff) from the sawtooth; R overshoots
        # by the forward extrapolation of a declining L plus the post-death tail.
        assert dose_L <= dose_ex * (1.0 + 1e-6)
        assert dose_R >= dose_ex * (1.0 - 1e-6)
        assert dose_R / dose_ex < 1.05
print("  => R's total-dose excess is bounded by the forward lag of a declining L")
print("     (Delta_star/t_age per step) plus the post-death tail of")
print("     TAG_LIFETIME_INTERVALS star steps: a few percent here. L's deficit is")
print("     the sawtooth, order unity when Delta_star >> tau.")
print()
# ---------------------------------------------------------------------------
# Part D: several stars on different bins, 0-D: summed rate vs dose reservoir
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part D: two stars, Delta_A = M * Delta_B, one gas particle on Delta_B/N")
print("=" * 78)
# Emission is accumulated as a pure dose (kappa = 0, e = 1, phi = 1) so the
# time-integrated injection can be compared to what the stars emitted exactly.
for M in (4, 16):
    for N in (1, 4):
        Delta_B = 1.0
        Delta_A = M * Delta_B
        delta = Delta_B / N
        S_A, S_B = 1.0, 0.3
        n_A_steps = 6
        n_sub = int(round(n_A_steps * Delta_A / delta))
        # Scheme "summed rate with reset on first touch of a step": the pure
        # rate form as first drafted; each star writes its rate on its own steps.
        rate = 0.0
        touch_ti = -1
        u_rate = 0.0
        # Scheme "dose reservoir with horizon": each star adds its elapsed-step
        # dose and extends the horizon by its own step; the particle drains
        # f = min(1, delta / t_remaining) of the reservoir per sub-step.
        D = 0.0
        horizon = -1.0
        u_res = 0.0
        emitted = 0.0
        for k in range(n_sub):
            t0 = k * delta
            t1 = t0 + delta
            # gas sub-step: drawdown at its start (drift), then apply.
            if D > 0.0:
                t_rem = horizon - t0
                f = 1.0 if t_rem <= delta else delta / t_rem
                applied = f * D
                D -= applied
                u_res += applied
            u_rate += rate * delta
            # star touches at the END of their own steps (after cooling/kick2).
            touched = False
            for S_s, Delta_s in ((S_A, Delta_A), (S_B, Delta_B)):
                if abs(t1 / Delta_s - round(t1 / Delta_s)) < 1e-9 and t1 > 0:
                    emitted += S_s * Delta_s
                    D += S_s * Delta_s
                    horizon = max(horizon, t1 + Delta_s)
                    if touch_ti != k:
                        rate = 0.0
                        touch_ti = k
                    rate += S_s
                    touched = True
        # Exact bookkeeping: what was emitted up to the last star touch equals
        # what the reservoir has applied plus what it still holds.
        print(
            f"  M={M:2d} N={N}: emitted {emitted:.3f};"
            f" reservoir applied+held {u_res + D:.3f} (held {D:.3f});"
            f" summed-rate applied {u_rate:.3f} ({u_rate / emitted:.3f} x emitted)"
        )
        assert abs((u_res + D) / emitted - 1.0) < 1e-12
        assert u_rate < 0.9 * emitted
    # Single star: the reservoir gives a constant rate equal to S (the R form).
    D, horizon, rates = S_A * Delta_A, Delta_A, []
    for k in range(int(Delta_A / delta)):
        t0 = k * delta
        t_rem = horizon - t0
        f = 1.0 if t_rem <= delta else delta / t_rem
        rates.append(f * D / delta)
        D -= f * D
    assert np.allclose(rates, S_A) and abs(D) < 1e-12
print("  => The summed rate with reset-on-first-touch loses the coarse star's")
print("     emission on the fine star's own steps; the dose reservoir applies every")
print("     deposited dose exactly once, for any bins, and reduces to a constant")
print("     rate S for a single star.")
print()

print("ALL CHECKS PASSED: the two schemes are not energetically equivalent;")
print("R reproduces the exact relaxation solution at every gas sub-step, L only")
print("at the post-injection instant (cooling never reads that instant); in the")
print("stiff regime L's between-lump remnant is e-suppressed and positive at the")
print("centre (the continuum-Laplacian estimate of a negative there was 20x too")
print("large); in the thin regime L with N >= 8 sub-steps per star step")
print("undershoots by 6-17% of its peak with no density contrast at all.")
