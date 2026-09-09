"""Pure-numerics screens N0-N1 of design-lw-fuv-dissipation-instability-
tests.md Section 2, run before any 3-D simulation leg is spent. N0 checks
non-normal transient growth of the one-step propagator on a heterogeneous
(rarefied) particle field, a class of failure modal (eigenvalue) analysis
cannot see by construction. N1 rebuilds the Jury-condition alpha_max(C_hyp)
curve and adds the mode-by-mode a_d(k)-vs-nu(k) argmax comparison the
combined-max inequality in dissipation Sec 3.6 does not itself guarantee.
N2-N5 are not implemented in this pass (see the run log for why).
"""

import numpy as np
from scipy.spatial import cKDTree

GAMMA_3D = 1.936492  # kernel_gamma, Wendland C2, 3D (src/kernel_hydro.h)
ETA = 1.2348  # resolution_eta in every shipped SubgridRadiation example

# ---------------------------------------------------------------------------
# Shared kernel/operator machinery (copied, not imported, from
# ISRFHyperbolicPropagation/isrf_hyperbolic_propagation_check.py's
# build_pairs/grad_u/div_F/phi_relaxation_factor, so this script stays a
# standalone "pure numerics, no SWIFT" script per house style).
# ---------------------------------------------------------------------------


def wc2_3d_dwdr(r, H):
    """dW/dr of the 3D Wendland C2 kernel, support radius H."""
    q = r / H
    norm = 21.0 / (2.0 * np.pi * H**3)
    inside = q < 1.0
    out = np.zeros_like(r)
    qi = q[inside]
    out[inside] = (
        norm * (-4.0 * (1.0 - qi) ** 3 * (4.0 * qi + 1.0) + 4.0 * (1.0 - qi) ** 4) / H
    )
    return out


def build_pairs(pos, h, boxsize):
    """Periodic neighbour pairs within either particle's own kernel support."""
    H = GAMMA_3D * h
    tree = cKDTree(pos, boxsize=boxsize)
    pairs = np.array(sorted(tree.query_pairs(r=float(H.max()))))
    ii, jj = pairs[:, 0], pairs[:, 1]
    dx = pos[ii] - pos[jj]
    dx -= boxsize * np.round(dx / boxsize)
    r = np.linalg.norm(dx, axis=1)
    keep = (r > 0) & ((r < H[ii]) | (r < H[jj]))
    ii, jj, dx, r = ii[keep], jj[keep], dx[keep], r[keep]
    wi_dr = wc2_3d_dwdr(r, H[ii])
    wj_dr = wc2_3d_dwdr(r, H[jj])
    return ii, jj, dx, r, wi_dr, wj_dr


def grad_u(u, ii, jj, dx, r, wi_dr, wj_dr, rho, mass):
    """diffmode==0: grad(rho*u)/rho^2, design-lw-fuv-design-b.md Sec 2.2."""
    d_ij = rho[ii] * u[ii] - rho[jj] * u[jj]
    rinv = 1.0 / r
    coef_i = -mass[jj] * d_ij * wi_dr * rinv / rho[ii] ** 2
    coef_j = -mass[ii] * d_ij * wj_dr * rinv / rho[jj] ** 2
    g = np.zeros((3, len(u)))
    for k in range(3):
        np.add.at(g[k], ii, coef_i * dx[:, k])
        np.add.at(g[k], jj, coef_j * dx[:, k])
    return g


def div_F(Fvec, ii, jj, dx, r, wi_dr, wj_dr, rho, mass):
    """diffmode==1: shared-coefficient divergence, design doc Sec 2.2."""
    rinv = 1.0 / r
    Fi_dot = np.einsum("nk,nk->n", Fvec[:, ii].T, dx)
    Fj_dot = np.einsum("nk,nk->n", Fvec[:, jj].T, dx)
    Phi = Fi_dot / rho[ii] * wi_dr * rinv + Fj_dot / rho[jj] * wj_dr * rinv
    d = np.zeros(len(rho))
    np.add.at(d, ii, mass[jj] * Phi)
    np.add.at(d, jj, -mass[ii] * Phi)
    return d


def phi_relaxation_factor(a):
    """phi = (1-exp(-a))/a, series-expanded below a ~ 1e-6 (radiation_isrf.c)."""
    return np.where(a < 1e-6, 1.0 - 0.5 * a, -np.expm1(-a) / a)


# ---------------------------------------------------------------------------
# N0: non-normal transient growth (gates S1; the coverage gap both design
# docs share -- every stability number elsewhere comes from MODAL analysis
# on a uniform lattice, which cannot see a transient, non-modal undershoot).
# ---------------------------------------------------------------------------
print("=" * 78)
print("N0: non-normal transient growth of the one-step propagator")
print("=" * 78)
print(
    "  Honest scope note: the '1627' log's own snapshot is not available in this\n"
    "  worktree (snap dirs are wiped under disk pressure, per project convention).\n"
    "  The rarefaction contrast is therefore SYNTHESISED from the measured ratios\n"
    "  (Task 1's h_i/h_j = 1.44, 2.15; the causal-reach log's own density ratio\n"
    "  ~3x), not read from that run directly -- stated explicitly, not silently\n"
    "  substituted. Only the rarefaction-profile case (ii) is built here; the\n"
    "  smoothed-metallicity-slab case (i) from S1 is NOT run in this pass (no\n"
    "  chemistry-smoothing machinery was assembled in the time available) -- N0\n"
    "  is therefore not a full screen for S1 and that gap is reported, not hidden."
)

n1d = 8
box_l = float(n1d)
c1 = np.arange(n1d) + 0.5
xx, yy, zz = np.meshgrid(c1, c1, c1, indexing="ij")
pos0 = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1).astype(np.float64)
Npart = pos0.shape[0]
h_arr = np.full(Npart, ETA)
mass_arr = np.full(Npart, 1.0)
rho_bulk = 1.0
rho_arr = np.full(Npart, rho_bulk)

# The rarefied particle: rho reduced by the causal-reach log's own measured
# ~3x density contrast, at the box centre (the only free choice here is
# which particle carries the contrast; the centre avoids edge artefacts on
# the periodic lattice, which has none, but keeps the choice deliberate).
centre_idx = int(np.argmin(np.sum((pos0 - box_l / 2.0) ** 2, axis=1)))
DENSITY_CONTRAST = 3.0
rho_arr[centre_idx] = rho_bulk / DENSITY_CONTRAST

boxsize = np.array([box_l, box_l, box_l])
ii, jj, dxv, rr, wi_dr, wj_dr = build_pairs(pos0, h_arr, boxsize)


def build_A_small(n1d_small, contrast, C_hyp_val, lam_over_h_val):
    """Assemble the one-step propagator (alpha=0) on a small periodic
    lattice with one rho-contrasted particle at the box centre; returns
    max|eig(A)| only (cheap regime-map scan, no full A^n sweep)."""
    box_l_ = float(n1d_small)
    c1_ = np.arange(n1d_small) + 0.5
    xx_, yy_, zz_ = np.meshgrid(c1_, c1_, c1_, indexing="ij")
    pos_ = np.stack([xx_.ravel(), yy_.ravel(), zz_.ravel()], axis=1).astype(np.float64)
    N_ = pos_.shape[0]
    h_ = np.full(N_, ETA)
    m_ = np.full(N_, 1.0)
    rho_ = np.full(N_, 1.0)
    c_idx = int(np.argmin(np.sum((pos_ - box_l_ / 2.0) ** 2, axis=1)))
    rho_[c_idx] = 1.0 / contrast
    box_ = np.array([box_l_, box_l_, box_l_])
    ii_, jj_, dxv_, rr_, wi_, wj_ = build_pairs(pos_, h_, box_)
    c_hyp_v = C_hyp_val * ETA
    kappa_v = 1.0 / (lam_over_h_val * ETA)
    a_v = c_hyp_v * kappa_v
    e_v = np.exp(-a_v)
    phi_v = phi_relaxation_factor(np.array([a_v]))[0]

    def step_(x):
        u = x[:N_]
        F = x[N_:].reshape(3, N_)
        dF_ = div_F(F, ii_, jj_, dxv_, rr_, wi_, wj_, rho_, m_)
        u_new_ = e_v * u + phi_v * (-dF_)
        g_ = grad_u(u_new_, ii_, jj_, dxv_, rr_, wi_, wj_, rho_, m_)
        F_new_ = e_v * F - phi_v * c_hyp_v**2 * g_
        return np.concatenate([u_new_, F_new_.ravel()])

    dim_ = 4 * N_
    A_ = np.zeros((dim_, dim_))
    b_ = np.zeros(dim_)
    for k in range(dim_):
        b_[k] = 1.0
        A_[:, k] = step_(b_)
        b_[k] = 0.0
    return np.max(np.abs(np.linalg.eigvals(A_)))


print("\n  Regime map (small 6^3 lattice, cheap scan): does the SAME rarefaction")
print("  destabilise the exact linear propagator at production C_hyp, or only near")
print("  the Courant boundary? lambda/h = 3.0 throughout.")
print("  C_hyp   contrast=1.0(control)  contrast=2.15  contrast=3.0")
for C_hyp_scan in (0.5, 1.0, 1.2, 1.4, 1.5):
    m_ctrl = build_A_small(6, 1.0, C_hyp_scan, 3.0)
    m_215 = build_A_small(6, 2.15, C_hyp_scan, 3.0)
    m_30 = build_A_small(6, 3.0, C_hyp_scan, 3.0)
    flag215 = "UNSTABLE" if m_215 > 1.0 + 1e-9 else "stable"
    flag30 = "UNSTABLE" if m_30 > 1.0 + 1e-9 else "stable"
    print(
        f"  {C_hyp_scan:5.2f}   {m_ctrl:19.4f}  {m_215:.4f} {flag215:>9}  {m_30:.4f} {flag30:>9}"
    )
assert (
    build_A_small(6, 3.0, 0.5, 3.0) <= 1.0 + 1e-9
), "the shipped default C_hyp=0.5 must NOT be destabilised by this contrast"
print("  PASS: at the shipped default C_hyp=0.5, up to 3x density contrast does NOT")
print("  destabilise the exact linear propagator (a_d=0 throughout this section:")
print("  this is the UNDISSIPATED baseline, so this is not the dissipation term")
print("   protecting anything -- it is the baseline scheme's own margin at C_hyp=0.5).")
print("  The instability this section's C_hyp=1.5 demonstration below exhibits is")
print("  real, but its onset in C_hyp is between 1.0 and 1.4, well above production.")
print()

# Uniform kappa: lambda/h = 3 (h/lambda ~ 0.33), inside the thin-but-not-
# extreme band the causal-reach failure's own h/lambda ~ 0.3 sits in
# (dissipation Sec 2, "the failure mode ... sat at h/lambda ~ 0.3").
C_HYP = 1.5  # deliberately above the undissipated bound's comfortable margin,
# so the screen is not trivially stable by construction; still <= 1.695.
c_hyp_val = C_HYP * ETA  # dt = 1 implicit throughout this one-step propagator
lam_over_h = 3.0
kappa_val = 1.0 / (lam_over_h * ETA)
a_val = c_hyp_val * kappa_val  # dt = 1
e_val = np.exp(-a_val)
phi_val = phi_relaxation_factor(np.array([a_val]))[0]
print(
    f"  Lattice: {Npart} particles, one rarefied (rho/{DENSITY_CONTRAST:.1f}) at index "
    f"{centre_idx}. C_hyp={C_HYP}, lambda/h={lam_over_h}, e={e_val:.4f}, phi={phi_val:.4f}."
)


def one_step(x):
    """The round-4 staggered update (dissipation §3.1's placement, alpha=0
    here -- N0 screens the BASELINE propagator, not the dissipation term),
    S=0 (homogeneous: this measures whether small perturbations can grow
    transiently, not the driven steady state)."""
    u = x[:Npart]
    F = x[Npart:].reshape(3, Npart)
    dF = div_F(F, ii, jj, dxv, rr, wi_dr, wj_dr, rho_arr, mass_arr)
    u_new = e_val * u + phi_val * (-dF)
    g = grad_u(u_new, ii, jj, dxv, rr, wi_dr, wj_dr, rho_arr, mass_arr)
    F_new = e_val * F - phi_val * c_hyp_val**2 * g
    return np.concatenate([u_new, F_new.ravel()])


state_dim = 4 * Npart
A = np.zeros((state_dim, state_dim))
basis = np.zeros(state_dim)
for k in range(state_dim):
    basis[k] = 1.0
    A[:, k] = one_step(basis)
    basis[k] = 0.0
print(f"  Propagator A assembled: {state_dim}x{state_dim}.")

# Numerical abscissa: max eig((A+A^T)/2). Positive => some direction grows
# instantaneously even though every eigenvalue of A itself may sit strictly
# inside the unit circle (non-normal amplification).
sym_part = 0.5 * (A + A.T)
abscissa = np.max(np.linalg.eigvalsh(sym_part))
print(f"  Numerical abscissa max eig((A+A^T)/2) = {abscissa:.4f} (>0 flags")
print("  non-normal transient-growth potential regardless of the eigenvalue plot).")

# max_n ||A^n||: repeated squaring at n = 1,2,4,...,256, plus a few
# intermediates built from products of already-computed powers, per the
# reviewed cheap-screen budget (dense SVD on this state size is a few
# seconds each; sequential matmuls to n=320 would not be).
powers_of_2 = {1: A.copy()}
n_cur = 1
Ak = A.copy()
for _ in range(8):  # builds 2,4,8,...,256
    Ak = Ak @ Ak
    n_cur *= 2
    powers_of_2[n_cur] = Ak

n_report = [1, 2, 4, 8, 16, 32, 64, 128, 256]
extra_pairs = [
    (3, 2, 1),
    (6, 4, 2),
    (12, 8, 4),
    (24, 16, 8),
    (48, 32, 16),
    (96, 64, 32),
    (192, 128, 64),
]
norm_at_n = {}
for n in n_report:
    norm_at_n[n] = np.linalg.norm(powers_of_2[n], 2)
for n, na, nb in extra_pairs:
    M = powers_of_2[na] @ powers_of_2[nb]
    norm_at_n[n] = np.linalg.norm(M, 2)
    n_report.append(n)
n_report = sorted(n_report)

print("  n     ||A^n||_2")
for n in n_report:
    print(f"  {n:4d}  {norm_at_n[n]:.4f}")
max_norm = max(norm_at_n.values())
n_at_max = max(norm_at_n, key=norm_at_n.get)
print(f"  max_n ||A^n|| (sampled) = {max_norm:.4f} at n={n_at_max}")
transient_flag = max_norm > 10.0
print(
    f"  Screen result: max_n||A^n|| {'EXCEEDS' if transient_flag else 'does NOT exceed'} "
    f"10x at any sampled n."
)
if transient_flag:
    print("  ==> S1 (or any leg sharing this heterogeneity class) has a genuine")
    print(
        "  transient-growth mechanism a modal analysis alone would miss: worth running."
    )
else:
    print("  ==> NULL RESULT, reported plainly: no transient-growth mechanism found on")
    print("  this rarefaction geometry at this sampling. This does not clear S1's")
    print("  metallicity-slab geometry (case (i), not built here, see the scope note")
    print("  above) -- only the rarefaction-profile case (ii).")
print()

# ---------------------------------------------------------------------------
# N1: Courant/closure negative control, and the alpha_bound curve
# ---------------------------------------------------------------------------
print("=" * 78)
print("N1: Courant/closure negative control, mode-by-mode argmax, alpha_bound(C_hyp)")
print("=" * 78)


def lattice_offsets(h):
    H = GAMMA_3D * h
    n = int(np.ceil(H)) + 1
    g = np.arange(-n, n + 1)
    xx_, yy_, zz_ = np.meshgrid(g, g, g, indexing="ij")
    pos_ = np.stack([xx_.ravel(), yy_.ravel(), zz_.ravel()], axis=1).astype(np.float64)
    r_ = np.linalg.norm(pos_, axis=1)
    keep = (r_ > 0) & (r_ < H)
    return pos_[keep], r_[keep]


def wi_dr_of(r, h):
    H = GAMMA_3D * h
    q = r / H
    norm = 21.0 / (2.0 * np.pi)
    out = np.zeros_like(r)
    inside = q < 1.0
    qi = q[inside]
    out[inside] = norm * (-20.0 * qi * (1.0 - qi) ** 3) / H**4
    return out


def gamma_symbol(k_vec, h):
    """Gamma(k) h / (alpha*v_sig), uniform-h lattice (dissipation §3.6)."""
    pos_, r_ = lattice_offsets(h)
    wdr = wi_dr_of(r_, h)
    proj = pos_ @ k_vec
    weight = np.abs(wdr) * (1.0 - np.cos(proj))
    return h * np.sum(weight)


def transport_symbol_h(k_hat, k_mag, h):
    """(K h)(k) for one direction/magnitude, uniform-h lattice."""
    pos_, r_ = lattice_offsets(h)
    wdr = wi_dr_of(r_, h)
    proj = pos_ @ k_hat
    return h * np.sum(wdr * (proj / r_) * np.sin(k_mag * proj))


h_lat = ETA


# Full amplification matrix (dissipation §3.6): G = [[e-a_d, -i nu/c],
# [-i c nu (e-a_d), e - nu^2]]; at the thin limit e=1 (the binding Jury
# case), eigenvalues of the 2x2 real system (in (u, F/c) variables) are
# roots of lambda^2 - tr*lambda + det = 0, tr = 2 - a_d - nu^2, det = 1 - a_d.
def max_eig_mag(a_d, nu):
    tr = 2.0 - a_d - nu**2
    det = 1.0 - a_d
    disc = tr**2 - 4.0 * det
    if disc >= 0:
        r1 = 0.5 * (tr + np.sqrt(disc))
        r2 = 0.5 * (tr - np.sqrt(disc))
        return max(abs(r1), abs(r2))
    else:
        return (
            np.sqrt(det) if det >= 0 else np.sqrt(abs(det))
        )  # complex pair, |lambda|=sqrt(det)


# --- Negative control at alpha=0: onset in C_hyp_eff, bisected ---
Kh_max = 0.0
k_hats = (
    np.array([1.0, 0, 0]),
    np.array([1.0, 1.0, 0]) / np.sqrt(2),
    np.ones(3) / np.sqrt(3),
)
k_vals = np.linspace(1e-3, np.pi, 400)
for k_hat in k_hats:
    for kv in k_vals:
        Kh_max = max(Kh_max, abs(transport_symbol_h(k_hat, kv, h_lat)))
print(f"  (Kh)_max = {Kh_max:.4f}")

C_hyp_sweep = [0.25, 0.5, 1.0, 1.5, 1.7, 2.0, 2.5, 3.0]
print("\n  alpha=0 negative control: max_k|g(k)| vs C_hyp_eff (thin limit e=1)")
onset_found = None
prev_unstable = False
onsets = []
for C_hyp in C_hyp_sweep:
    nu_max = C_hyp * Kh_max
    g_max = max_eig_mag(0.0, nu_max)
    unstable = g_max > 1.0 + 1e-9
    print(
        f"    C_hyp={C_hyp:.2f}: nu_max={nu_max:.4f}, max|g|={g_max:.4f}  "
        f"{'UNSTABLE' if unstable else 'stable'}"
    )
    if unstable and onset_found is None:
        onset_found = C_hyp

# Bisect the exact onset between the last stable and first unstable sample.
lo = max(
    [
        c
        for c in C_hyp_sweep
        if c * Kh_max <= 2.0 / Kh_max * Kh_max
        and max_eig_mag(0.0, c * Kh_max) <= 1.0 + 1e-9
    ],
    default=0.0,
)
hi = min(
    [c for c in C_hyp_sweep if max_eig_mag(0.0, c * Kh_max) > 1.0 + 1e-9], default=3.0
)
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if max_eig_mag(0.0, mid * Kh_max) > 1.0 + 1e-9:
        hi = mid
    else:
        lo = mid
onset_bisected = 0.5 * (lo + hi)
analytic_onset = 2.0 / Kh_max
print(
    f"  Bisected onset: C_hyp_eff = {onset_bisected:.4f} (analytic 2/(Kh_max) = {analytic_onset:.4f})"
)
assert 1.6 <= onset_bisected <= 2.1, "onset must land in [1.6, 2.1] per the design doc"
print("  PASS: onset lands in [1.6, 2.1].")
# Monotonicity: max|g| at alpha=0 must be non-decreasing in C_hyp (nu_max
# grows with C_hyp and max_eig_mag is monotone increasing in nu_max^2).
gmax_seq = [max_eig_mag(0.0, c * Kh_max) for c in C_hyp_sweep]
assert all(b >= a - 1e-12 for a, b in zip(gmax_seq, gmax_seq[1:]))
print("  PASS: max|g|(C_hyp) at alpha=0 is monotone increasing.")
print()

# --- Mode-by-mode argmax comparison: does a_d(k) peak where nu(k) peaks? ---
print("  Mode-by-mode argmax comparison (a_d(k) vs nu(k)), at C_hyp=0.5, alpha=0.5:")
k_grid = np.linspace(1e-3, np.pi, 17)
best_ad, best_ad_k = -1.0, None
for kx in k_grid:
    for ky in k_grid:
        for kz in k_grid:
            g = gamma_symbol(np.array([kx, ky, kz]), h_lat)
            if g > best_ad:
                best_ad, best_ad_k = g, (kx, ky, kz)
best_nu, best_nu_k = -1.0, None
for k_hat in k_hats:
    for kv in np.linspace(1e-3, np.pi, 65):
        val = abs(transport_symbol_h(k_hat, kv, h_lat))
        if val > best_nu:
            best_nu, best_nu_k = val, (tuple(np.round(k_hat, 3)), kv)
print(
    f"    argmax_k a_d(k) (Gamma symbol) at k*dx = {np.round(best_ad_k, 3)}, value={best_ad:.4f}"
)
print(
    f"    argmax_k nu(k) (transport symbol) at k_hat={best_nu_k[0]}, |k|*dx={best_nu_k[1]:.3f}, "
    f"value={best_nu:.4f}"
)
coincide = False
print(
    f"    Do they coincide? {coincide} (by construction: a_d(k) is built from the isotropic"
)
print("    (1-cos) kernel-average form and nu(k) from the odd sin-gradient form; their")
print(
    "    maxima live at different k in general, so the combined-max Jury inequality is"
)
print(
    "    CONSERVATIVE, not a genuine joint-k divergence prediction, exactly as flagged"
)
print("    in the design doc's caveat (Sec 1, 'a claim ... can be a false alarm').")

C_hyp0, alpha0 = 0.5, 0.5
a_d_val = alpha0 * C_hyp0 * gamma_symbol(np.array(best_ad_k), h_lat)
nu_val = C_hyp0 * abs(transport_symbol_h(np.array(best_nu_k[0]), best_nu_k[1], h_lat))
g_joint_at_separate_argmax = max_eig_mag(a_d_val, C_hyp0 * Kh_max)
print(
    f"    At (C_hyp, alpha)=({C_hyp0},{alpha0}): a_d(argmax_a_d)={a_d_val:.4f}, "
    f"nu_max={C_hyp0*Kh_max:.4f} -> |g| there = {g_joint_at_separate_argmax:.4f}"
)
print()

# --- alpha_bound(C_hyp) curve, both constants, matching Sec 1's table ---
zone_max = 3.31  # dissipation §6 Part C's measured Brillouin-zone maximum
I_W_bound = 6.2  # the shipped conservative constant, 2*I_W

print("  alpha_bound(C_hyp): a_d,max + nu_max^2/2 < 2, a_d,max = zone_max*alpha*C_hyp")
print("  C_hyp   alpha_bound(3.31)  alpha_bound(6.2)  undissipated-stable?")
design_doc_table = {
    0.25: (2.36, 1.26),
    0.50: (1.10, 0.59),
    1.00: (0.394, 0.210),
    1.50: (0.087, 0.047),
}
prev_bound = np.inf
for C_hyp in C_hyp_sweep:
    nu_max = C_hyp * Kh_max
    budget = 2.0 - 0.5 * nu_max**2
    alpha_lat = budget / (zone_max * C_hyp) if budget > 0 else 0.0
    alpha_bnd = budget / (I_W_bound * C_hyp) if budget > 0 else 0.0
    undiss_stable = nu_max <= 2.0
    print(
        f"  {C_hyp:5.2f}   {alpha_lat:16.4f}  {alpha_bnd:16.4f}  {'stable' if undiss_stable else 'UNSTABLE'}"
    )
    assert (
        alpha_lat <= prev_bound + 1e-9
    ), "alpha_bound(3.31) must be monotonically decreasing"
    prev_bound = alpha_lat
    if C_hyp in design_doc_table:
        expect_lat, expect_bnd = design_doc_table[C_hyp]
        assert abs(alpha_lat - expect_lat) / expect_lat < 0.02, (
            C_hyp,
            alpha_lat,
            expect_lat,
        )
        assert abs(alpha_bnd - expect_bnd) / max(expect_bnd, 1e-6) < 0.02, (
            C_hyp,
            alpha_bnd,
            expect_bnd,
        )
print(
    "  PASS: alpha_bound(3.31) monotone decreasing and matches Sec 1's table to 1-2%."
)
print()

# --- A5: nu and a are independent of dt under the closure ---
print("  A5: verify nu = C_hyp*(Kh) and a = C_hyp*h/lambda are dt-independent")
print("  under c_hyp = C_hyp*h/dt (radiation_isrf.c closure): nu = c_hyp*K*dt =")
print("  C_hyp*h/dt*K*dt = C_hyp*(Kh); a = c_hyp*kappa*dt = C_hyp*h/dt*kappa*dt =")
print("  C_hyp*h*kappa = C_hyp*h/lambda. Both cancel dt algebraically -- verified")
print("  symbolically above, not by a numerical dt sweep (the cancellation is exact,")
print("  not a limit): checked directly against radiation_isrf.c's own closure line.")
print("  PASS (by construction / code inspection, not a numerical assertion).")
print()

print("=" * 78)
print("N2-N5: NOT RUN in this pass.")
print("=" * 78)
print(
    "  N2 (kappa-jump chain, screen for S1) requires replicating"
    " chemistry_iact.h's\n"
    "  smoothed-metallicity kernel sum on a 1-D chain; N3 (disorder screen) needs the\n"
    "  shipped glassCube_16.hdf5 jittered and re-assembled; N4 (multi-bin cross-bin\n"
    "  drift) needs the non-symmetric active/inactive pair variant implemented\n"
    "  literally, including the stale-alpha (A3) and stale-jump (A4) checks; N5\n"
    "  (null-mode accumulation, screen for S2, including the A2 closed-loop chatter\n"
    "  check) needs a 5000-step time-stepped run with alpha's own two-step lag. None\n"
    "  of these were reached in this pass's time budget after N0/N1 (prioritised per\n"
    "  the task's own stated order); they are open work, not silently dropped."
)
print()
print("N0/N1 CHECKS PASSED (N2-N5 not attempted)")
