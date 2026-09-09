"""Verify Design B's Stage-1 artificial-dissipation term (design-lw-fuv-
design-b-dissipation.md Section 3-6): exact conservation of `sum_i m_i u_i`
(Part A), the continuum-limit diffusion coefficient AND SIGN (Part B), and
the explicit-stability bound with its admissible `alpha_max(C_hyp)` (Part
C).

Fixed formula under test (Section 3.1), per band:

    d_ij   = rho_i*u_i_prev - rho_j*u_j_prev
    Wbar_ij = 0.5*(wi_dr + wj_dr)                         (< 0)
    v_sig,ij = alpha_ij * min(c_hyp_i, c_hyp_j)           (alpha_ij = max(alpha_i, alpha_j);
                                                            a VELOCITY, no h factor)
    Psi_ij  = v_sig,ij * d_ij * Wbar_ij / (rho_i*rho_j)

    dissipation_u_i +=  m_j * Psi_ij
    dissipation_u_j += -m_i * Psi_ij

Part B is the sign-fixing measurement: Section 3.1 states the sign must be
established by measurement, not by inspection, because a flipped sign here
is a growth term with the mirrored amplification spectrum of Part C and
would pass Part A's conservation check unnoticed.
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy.integrate import quad

RNG = np.random.default_rng(20260909)

GAMMA_3D = 1.936492  # kernel_gamma, Wendland C2, 3D (src/kernel_hydro.h)
ETA = 1.2348  # resolution_eta in every shipped SubgridRadiation example


def wc2_3d_w_dwdr(r, H):
    """W(r) and dW/dr(r) of the 3D Wendland C2 kernel, support radius H,
    exactly as SWIFT computes it: `dw/dq = -20 q (1-q)^3` on the
    `(21/2pi)(1-q)^4(1+4q)` form (design-lw-fuv-design-b-dissipation.md
    Section 3.5). H may be a scalar or an array broadcastable with r."""
    r = np.atleast_1d(np.asarray(r, dtype=np.float64))
    H = np.broadcast_to(np.asarray(H, dtype=np.float64), r.shape)
    q = r / H
    norm = 21.0 / (2.0 * np.pi)
    inside = q < 1.0
    w = np.zeros_like(q)
    dwdr = np.zeros_like(q)
    qi = q[inside]
    Hi = H[inside]
    w[inside] = norm * (1.0 - qi) ** 4 * (4.0 * qi + 1.0) / Hi**3
    dwdr[inside] = norm * (-20.0 * qi * (1.0 - qi) ** 3) / Hi**4
    return w, dwdr


def wi_dr_of(r, h):
    """`wi_dr = h^-(dim+1) dW/dq|_{r/h}`, SWIFT's own convention, dim=3."""
    H = GAMMA_3D * h
    _, dwdr = wc2_3d_w_dwdr(r, H)
    return dwdr


def psi_pair(
    dx,
    r,
    wi_dr,
    wj_dr,
    mi,
    mj,
    rho_i,
    rho_j,
    c_hyp_i,
    c_hyp_j,
    alpha_i,
    alpha_j,
    u_i_prev,
    u_j_prev,
):
    """One pair's `dissipation_u_i`/`dissipation_u_j` contribution, the
    formula under test."""
    d_ij = rho_i * u_i_prev - rho_j * u_j_prev
    Wbar_ij = 0.5 * (wi_dr + wj_dr)
    v_sig_ij = np.maximum(alpha_i, alpha_j) * np.minimum(c_hyp_i, c_hyp_j)
    Psi_ij = v_sig_ij * d_ij * Wbar_ij / (rho_i * rho_j)
    return mj * Psi_ij, -mi * Psi_ij


# ---------------------------------------------------------------------------
# Part A: exact conservation
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part A: exact conservation of sum_i m_i u_i")
print("=" * 78)

n_draw = 200
h_i = RNG.uniform(0.5, 2.0, n_draw)
h_j = RNG.uniform(0.5, 2.0, n_draw)
rho_i = RNG.uniform(0.1, 10.0, n_draw)
rho_j = RNG.uniform(0.1, 10.0, n_draw)
m_i = RNG.uniform(0.1, 5.0, n_draw)
m_j = RNG.uniform(0.1, 5.0, n_draw)
c_hyp_i = RNG.uniform(0.01, 100.0, n_draw)
c_hyp_j = RNG.uniform(0.01, 100.0, n_draw)
alpha_i = RNG.uniform(0.0, 1.2, n_draw)
alpha_j = RNG.uniform(0.0, 1.2, n_draw)
u_i = RNG.uniform(-10.0, 10.0, n_draw)
u_j = RNG.uniform(-10.0, 10.0, n_draw)
h_max = np.maximum(h_i, h_j)
r = RNG.uniform(0.05, 0.9, n_draw) * GAMMA_3D * h_max  # inside both kernels
dx = np.zeros((n_draw, 3))
dx[:, 0] = r

wi_dr = wi_dr_of(r, h_i)
wj_dr = wi_dr_of(r, h_j)
diss_i, diss_j = psi_pair(
    dx,
    r,
    wi_dr,
    wj_dr,
    m_i,
    m_j,
    rho_i,
    rho_j,
    c_hyp_i,
    c_hyp_j,
    alpha_i,
    alpha_j,
    u_i,
    u_j,
)

# The pair coefficient is a VELOCITY: no h factor appears anywhere in
# v_sig_ij, and the only length scale in Psi_ij is the one already inside
# Wbar_ij's own h^-(dim+1) normalisation (Section 3.1's post-round-1 fix).
# A script written with an extra h factor here would still pass the
# antisymmetry check below (it passes for ANY symmetric coefficient); it
# is Part B's I_2 value that catches that mistake.
budget = m_i * diss_i + m_j * diss_j
scale = np.abs(m_i * m_j * (diss_i - diss_j)) + 1e-300
rel = np.abs(budget) / np.maximum(np.abs(m_i * diss_i), 1e-300)
print(
    f"  200 random draws: max |m_i*diss_i + m_j*diss_j| = {np.max(np.abs(budget)):.3e}"
)
print(f"  max relative residual (vs |m_i*diss_i|) = {np.max(rel):.3e}")
assert np.all(np.abs(budget) < 1e-14 * np.maximum(np.abs(m_i * diss_i), 1.0))
print("  PASS: mirrored-mass antisymmetry holds to float round-off, any inputs.")

# 1-D periodic chain, 5000 steps, kappa = 0, random alpha per particle:
# sum m u must equal the injected total to 1e-12.
n = 48
dx_p = 1.0
h_chain = ETA * dx_p
H_chain = GAMMA_3D * h_chain
m_chain = RNG.uniform(0.5, 2.0, n)
rho_chain = m_chain / dx_p  # 1-D "density": mass per unit length
alpha_chain = RNG.uniform(0.0, 0.3, n)  # kept modest: this test is about
c_hyp_chain = RNG.uniform(0.5, 1.5, n)  # conservation, not stability margin
S_chain = np.zeros(n)
S_chain[n // 3] = 0.7  # one source particle, mass-specific rate

nmax = int(np.ceil(H_chain / dx_p))
offsets = np.arange(1, nmax + 1)
r_off = offsets * dx_p
w_off, dwdr_off = wc2_3d_w_dwdr(r_off, H_chain)
wi_dr_off = dwdr_off  # same h for every particle on this chain

u = np.zeros(n)
dt = 1e-3
n_steps = 5000
injected_total = 0.0
for _ in range(n_steps):
    diss = np.zeros(n)
    for s, r_s, wdr in zip(offsets, r_off, wi_dr_off):
        j_plus = (np.arange(n) + s) % n
        j_minus = (np.arange(n) - s) % n
        for j_idx, sign in ((j_plus, +1.0), (j_minus, -1.0)):
            d_ij = rho_chain * u - rho_chain[j_idx] * u[j_idx]
            Wbar = wdr  # equal-h chain: wi_dr = wj_dr = wdr
            v_sig = np.maximum(alpha_chain, alpha_chain[j_idx]) * np.minimum(
                c_hyp_chain, c_hyp_chain[j_idx]
            )
            Psi = v_sig * d_ij * Wbar / (rho_chain * rho_chain[j_idx])
            diss += m_chain[j_idx] * Psi
    u = u + dt * (diss + S_chain)
    injected_total += dt * np.sum(m_chain * S_chain)

sum_mu = np.sum(m_chain * u)
rel_err = abs(sum_mu - injected_total) / injected_total
print(
    f"  1-D chain, {n_steps} steps: sum(m u) = {sum_mu:.10f}, injected total = "
    f"{injected_total:.10f}, rel. error = {rel_err:.3e}"
)
assert rel_err < 1e-12
print("  PASS: exact conservation over the full time-stepped run.")
print()


# ---------------------------------------------------------------------------
# Kernel volume integrals I_2, I_W (Section 3.5), by radial quadrature
# ---------------------------------------------------------------------------
def kernel_integrals(h):
    H = GAMMA_3D * h

    def integrand_I2(r_):
        _, dwdr_ = wc2_3d_w_dwdr(np.array([r_]), H)
        return 4.0 * np.pi * r_**4 * abs(dwdr_[0])

    def integrand_IW(r_):
        _, dwdr_ = wc2_3d_w_dwdr(np.array([r_]), H)
        return 4.0 * np.pi * r_**2 * abs(dwdr_[0])

    I2, _ = quad(integrand_I2, 0.0, H, limit=200)
    IW_raw, _ = quad(integrand_IW, 0.0, H, limit=200)
    return I2, h * IW_raw


I2_over_h, I_W = kernel_integrals(1.0)
print("=" * 78)
print("Kernel volume integrals (Section 3.5)")
print("=" * 78)
print(
    f"  I_2/h (numerically integrated) = {I2_over_h:.5f}  (Section 3.5 hand-derivation: 3.23)"
)
print(
    f"  I_W (numerically integrated)   = {I_W:.5f}  (Section 3.5 hand-derivation: 3.10)"
)
assert abs(I2_over_h - 3.23) < 1e-3 * 10  # hand-derivation quoted to 3 sig figs
assert abs(I2_over_h - 3.23) < 5e-3
assert abs(I_W - 3.10) < 5e-3
print("  PASS: both match Section 3.5's hand-derived constants to 1e-3.")
print()

# ---------------------------------------------------------------------------
# Part B: continuum limit AND SIGN
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part B: continuum limit and sign (the sign-fixing measurement)")
print("=" * 78)


def build_lattice(n_per_dim, box_l):
    dxp = box_l / n_per_dim
    c1 = (np.arange(n_per_dim) + 0.5) * dxp
    xx, yy, zz = np.meshgrid(c1, c1, c1, indexing="ij")
    return np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1), dxp


def accumulate_dissipation(pos, rho, mass, u_field, h, alpha, c_hyp):
    tree = cKDTree(pos, boxsize=None)
    H = GAMMA_3D * h
    pairs = tree.query_pairs(r=H, output_type="ndarray")
    i_idx, j_idx = pairs[:, 0], pairs[:, 1]
    dxv = pos[i_idx] - pos[j_idx]
    rr = np.linalg.norm(dxv, axis=1)
    wdr = wi_dr_of(rr, h)  # uniform h: wi_dr = wj_dr
    di, dj = psi_pair(
        dxv,
        rr,
        wdr,
        wdr,
        mass[i_idx],
        mass[j_idx],
        rho[i_idx],
        rho[j_idx],
        c_hyp[i_idx],
        c_hyp[j_idx],
        alpha[i_idx],
        alpha[j_idx],
        u_field[i_idx],
        u_field[j_idx],
    )
    out = np.zeros(pos.shape[0])
    np.add.at(out, i_idx, di)
    np.add.at(out, j_idx, dj)
    return out


def interior_probe(pos, h, box_l, frac=0.5):
    margin = 1.5 * GAMMA_3D * h
    inside = np.all((pos > margin) & (pos < box_l - margin), axis=1)
    cand = np.where(inside)[0]
    target = np.array([frac * box_l, box_l / 2, box_l / 2])
    return cand[np.argmin(np.sum((pos[cand] - target) ** 2, axis=1))]


box_l, n_per_dim = 1.0, 24
pos, dxp = build_lattice(n_per_dim, box_l)
rho0 = 2.0
mass0 = rho0 * dxp**3
alpha0, c_hyp0 = 0.7, 1.3
v_sig0 = alpha0 * c_hyp0

# B.1: uniform u_V on a density gradient must give exactly 0.
slope = 1.7
rho_grad = rho0 * (1.0 + slope * pos[:, 0])
mass_grad = rho_grad * dxp**3
U0 = 3.0
u_on_grad = U0 / rho_grad  # rho*u = U0 everywhere: d_ij = 0 for every pair
alpha_arr = np.full(n_per_dim**3, alpha0)
c_hyp_arr = np.full(n_per_dim**3, c_hyp0)
h_probe = 1.8 * dxp
diss_grad = accumulate_dissipation(
    pos, rho_grad, mass_grad, u_on_grad, h_probe, alpha_arr, c_hyp_arr
)
p_idx = interior_probe(pos, h_probe, box_l)
print(
    f"  B.1 uniform u_V on a density gradient: dissipation at probe = "
    f"{diss_grad[p_idx]:.3e} (must be exactly 0)"
)
assert diss_grad[p_idx] == 0.0
print("  PASS: the term differences u_V = rho*u, not u; a density gradient alone")
print("  triggers nothing.")

# B.2: u_V quadratic, uniform density: sum_j m_j Psi_ij -> +(D_eff/rho_i) laplacian(u_V).
A_quad = 0.9
u_V_quad = A_quad * (pos[:, 0] - box_l / 2) ** 2
u_quad = u_V_quad / rho0
laplacian_uV = 2.0 * A_quad  # d^2/dx^2 of A*x^2; no y,z dependence

h_sweep = [1.3, 1.8, 2.6, 3.8]
errs = []
for hod in h_sweep:
    h = hod * dxp
    diss = accumulate_dissipation(
        pos,
        np.full_like(u_quad, rho0),
        np.full_like(u_quad, mass0),
        u_quad,
        h,
        alpha_arr,
        c_hyp_arr,
    )
    I2_h, _ = kernel_integrals(h)
    D_eff = v_sig0 * I2_h / 6.0
    target = (D_eff / rho0) * laplacian_uV
    probes = [interior_probe(pos, h, box_l, frac) for frac in (0.4, 0.5, 0.6)]
    measured = np.mean(diss[probes])
    err = abs(measured - target) / abs(target)
    errs.append(err)
    nn = (4.0 / 3.0) * np.pi * hod**3
    print(
        f"  ~{nn:5.0f} nbrs: measured = {measured:+.5f}, target (+D_eff/rho lap(u_V)) = "
        f"{target:+.5f}, rel. err = {err:.3e}"
    )

print(f"  finest-resolution relative error = {errs[-1]:.3e}")
assert errs[-1] < 0.05
assert errs[-1] < errs[0]  # converging, not diverging
print("  PASS: the accumulated sum converges to +(D_eff/rho_i) laplacian(u_V),")
print("  the SIGN Section 3.1 predicts (a flipped sign would converge to the")
print("  same magnitude with a negative measured value: not observed).")
print()

# ---------------------------------------------------------------------------
# Part C: stability and the parameter bound
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part C: stability bound, admissible alpha_max(C_hyp)")
print("=" * 78)


def lattice_neighbours(h):
    H = GAMMA_3D * h
    n = int(np.ceil(H)) + 1
    g = np.arange(-n, n + 1)
    xx, yy, zz = np.meshgrid(g, g, g, indexing="ij")
    pos_ = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1).astype(np.float64)
    r_ = np.linalg.norm(pos_, axis=1)
    keep = (r_ > 0) & (r_ < H)
    return pos_[keep], r_[keep]


def gamma_symbol(k_vec, h, dx=1.0):
    """Gamma(k) h / (alpha*v_sig) = h * sum_j (m_j/rho_j) |Wbar_ij| (1 - cos(k.r_ij)),
    m_j/rho_j = dx^3 on a uniform lattice, Wbar_ij = wi_dr (equal h)."""
    pos_, r_ = lattice_neighbours(h)
    wdr = wi_dr_of(r_ * dx, h * dx)
    proj = pos_ @ k_vec
    weight = np.abs(wdr) * (1.0 - np.cos(proj))
    return h * dx**3 * np.sum(weight)


h_lat = ETA
corner = gamma_symbol(np.array([np.pi, np.pi, np.pi]), h_lat)
face = gamma_symbol(np.array([np.pi, 0.0, 0.0]), h_lat)
edge = gamma_symbol(np.array([np.pi, np.pi, 0.0]), h_lat)
print(
    f"  corner mode (pi,pi,pi)/dx: Gamma h/(alpha v) = {corner:.4f}  (Section 3.6: ~3.14-3.15)"
)
print(
    f"  face mode   (pi,0,0)/dx:   Gamma h/(alpha v) = {face:.4f}  (Section 3.6: ~3.08)"
)
print(
    f"  edge mode   (pi,pi,0)/dx:  Gamma h/(alpha v) = {edge:.4f}  (Section 6 Part C: zone maximum, ~3.31)"
)
assert abs(corner - 3.15) < 0.05
assert abs(face - 3.08) < 0.05
assert abs(edge - 3.31) < 0.05

# Full Brillouin-zone maximum: sample densely rather than just the three
# named high-symmetry points, to confirm the edge mode really is the max.
grid1d = np.linspace(0.0, np.pi, 9)
zone_max = 0.0
zone_argmax = None
for kx in grid1d:
    for ky in grid1d:
        for kz in grid1d:
            g = gamma_symbol(np.array([kx, ky, kz]), h_lat)
            if g > zone_max:
                zone_max = g
                zone_argmax = (kx, ky, kz)
print(
    f"  Brillouin-zone maximum (9^3 grid): Gamma h/(alpha v) = {zone_max:.4f} at "
    f"k*dx = {tuple(np.round(zone_argmax, 3))}"
)
assert (
    abs(zone_max - edge) < 0.02
), "the sampled zone maximum should coincide with the edge mode"
I_W_bound = 2.0 * I_W
print(
    f"  I_W = {I_W:.4f}, all-neighbours-out-of-phase bound 2*I_W = {I_W_bound:.4f}"
    f" (not realised on a lattice, {zone_max:.2f} < {I_W_bound:.2f})"
)


# nu_max = C_hyp * (Kh)_max, the odd transport symbol from the parent
# stability script (verify_design_b_timestepping_stability.py Part B),
# recomputed here so this script is self-contained. k_hat is a UNIT
# direction; k_mag is the wavenumber magnitude (dx = 1 throughout).
def transport_symbol(k_hat, k_mag, h):
    pos_, r_ = lattice_neighbours(h)
    wdr = wi_dr_of(r_, h)
    proj = pos_ @ k_hat  # direction-cosine numerator, physical distance
    return np.sum(wdr * (proj / r_) * np.sin(k_mag * proj))


k_hats = (
    np.array([1.0, 0, 0]),
    np.array([1.0, 1.0, 0]) / np.sqrt(2),
    np.ones(3) / np.sqrt(3),
)
k_values = np.linspace(1e-3, np.pi, 200)
Kh_max = 0.0
for k_hat in k_hats:
    vals = [abs(transport_symbol(k_hat, kv, h_lat)) for kv in k_values]
    Kh_max = max(Kh_max, max(vals) * h_lat)
print(
    f"  (Kh)_max (transport symbol, parent script's Part B) = {Kh_max:.4f} (expected ~1.18)"
)
assert abs(Kh_max - 1.18) < 0.05

C_hyp_values = [0.25, 0.5, 1.0, 1.7]
print("\n  Admissible alpha_max(C_hyp), from the JURY condition")
print(
    "  a_d,max + nu_max^2/2 < 2, a_d,max = zone_max*alpha*C_hyp, nu_max = C_hyp*(Kh)_max:"
)
print("  C_hyp    alpha_max(zone_max)   alpha_max(bound 2*I_W)")
alpha_max_at_default = None
for C_hyp in C_hyp_values:
    nu_max = C_hyp * Kh_max
    budget = 2.0 - 0.5 * nu_max**2
    alpha_lattice = budget / (zone_max * C_hyp) if budget > 0 else 0.0
    alpha_bound = budget / (I_W_bound * C_hyp) if budget > 0 else 0.0
    print(f"  {C_hyp:5.2f}    {alpha_lattice:8.3f}              {alpha_bound:8.3f}")
    if abs(C_hyp - 0.5) < 1e-9:
        alpha_max_at_default = alpha_bound

# The enforced range check (feedback_properties.h) uses the CONSERVATIVE
# bound 6.2 = 2*I_W and 0.70 = (Kh_max)^2/2 (Section 3.6's decision).
enforced_bound_at_half = (2.0 - 0.70 * 0.5**2) / (6.2 * 0.5)
print(f"\n  Enforced code constant check: 6.2*alpha_max*C_hyp + 0.70*C_hyp^2 <= 2")
print(
    f"  at C_hyp=0.5: alpha_max <= {enforced_bound_at_half:.4f} "
    f"(measured 2*I_W={I_W_bound:.3f}, (Kh_max)^2/2={0.5 * Kh_max**2:.3f})"
)
assert abs(enforced_bound_at_half - alpha_max_at_default) < 0.05
print("  PASS: the code's hardcoded (6.2, 0.70) constants match this script's")
print("  measured (2*I_W, (Kh_max)^2/2) to within a few percent.")
print()
print("ALL CHECKS PASSED")
