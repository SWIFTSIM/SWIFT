"""Conserved-functional diagnostic for the `diffmode==2` anisotropic
gradient loop's skew-adjointness under the M1 pressure-tensor closure.

Background: the M1 upgrade replaced the earlier isotropic scalar gradient
term (which was the EXACT skew-adjoint of the density loop's `div(F)` in
the `m*rho` inner product) with the anisotropic pressure-tensor divergence
(`radiation_gradient_accumulate_band`, `diffmode==2`: a SHARED averaged
kernel derivative `wbar_dr = 0.5*(wi_dr+wj_dr)`, vs. the density loop's own
per-particle `wi_dr`/`wj_dr` in `radiation_divergence_accumulate_band`,
`diffmode==1`, unchanged). `sum(m*u)` itself is untouched (the density loop
is unchanged and exactly antisymmetric); what is at risk is a DIFFERENT,
subtler property: L2/quadratic-energy stability of the linearized coupled
(u, F) system.

This script re-implements both pairwise operators plus the exact-relaxation
`u`/`F` updates (radiation_isrf.c `radiation_end_density_propagation`,
`radiation_end_gradient_propagation`, at kappa=0/source=0 both reduce to
plain forward-Euler forms) and the M1 flux limiter
(`radiation_apply_flux_limiter_band`), directly from source (git commit
`d181a1c27`), and evolves a smooth low-amplitude `(u, F)` seed perturbation
on a FIXED background (pos/rho/mass/h frozen, exactly what the real
per-step density snapshot `rho_prev` does; only `(u, F)` evolve), with
absorption (`kappa`), dissipation (`alpha_max`/`alpha_floor`) and the
source term all off. This isolates the propagation-only linear/quasi-linear
operator that the skew-adjointness concern above is about.

Conserved-functional formula (P1/isotropic-limit convention,
`dF/dt = -c_hyp^2*G` with the `1/3` already folded into `G`):

    E = sum_i m_i*rho_i*(u_i^2/2 + 3*|F_i|^2/(2*c_hyp^2))

Discriminator: at UNIFORM h, `wi_dr == wj_dr` for every pair (same kernel,
same h), so `wbar_dr` reduces to the same per-particle value the old
adjoint pairing used, the residual is identically zero there BY
CONSTRUCTION, and this script's own uniform-h control run below is the
self-check that the re-implementation is faithful (E must be flat to
round-off with A_rho=0). Only a VARIABLE-h configuration can discriminate
the residual at all. `sum(m*u)` is checked at every resolution as an
independent correctness control on the (unchanged, exactly antisymmetric)
density loop.

Verdict rule: drift that SHRINKS WITH RESOLUTION (h -> 0) is an expected,
bounded consistency error, monitor only, no action needed. Drift that
GROWS IN TIME AT FIXED h (accelerating, not just linear accumulation) is a
real non-normal mode, STOP: this script's own authority is not sufficient
to clear it, diagnose the cause before trusting any downstream result.

STATUS: a parallel investigation found
that the isotropic `(u, F)` functional above is NOT a conserved quantity of
the M1 equations at all (even at uniform h, for a genuinely anisotropic
`D`), the residual splits into `R_h` (the kernel-derivative mismatch this
script was built to isolate, removable by a `diffmode==0` gradient-loop
fix) and `R_D` (an irreducible PDE-level anisotropy term, non-zero even at
uniform h, NOT fixable). A low-amplitude run stays near-isotropic
(`D ~ I/3`) throughout, so it under-exercises `R_D` and can read "clean"
while missing the larger effect, exactly what happened on this script's
first corrected run (see IC bug note below). One addition below, per that
investigation's own recommendation (its Section 8): the `D`-weighted
functional `Q_D = sum m*rho*(u^2/2 + F.D^-1.F/(2*c_hyp^2))` (the metric in
which frozen-`D` adjointness genuinely holds) alongside the isotropic `E`.
(A SECOND recommended addition, a uniform-h constant-anisotropic-`D`
static control that isolates `R_D` directly, is NOT implemented here;
the only uniform-h section below is a dynamical run with the real
closure, which stays near-isotropic; still open, not a blocker for the
results below.)

STATUS: the `diffmode==0` fix has LANDED (`radiation_propagation_iact.h`).
The numbers below are post-fix and can be trusted as a verdict on their
own terms.

IC BUG FIXED (found independently while first running this script): the
gradient/divergence pair operates on `u_V = rho*u`, not `u` itself (see
`verify_isrf_dissipation.py`'s own Part B.1: "uniform u_V on a density
gradient must give exactly 0"). An earlier version of this script perturbed
`u` directly with `rho(x)` varying underneath it, which makes `u_V`
strongly non-uniform (inheriting the full density-modulation amplitude,
not the intended small perturbation), an O(1) jump was mistaken for a
"low-amplitude seed", explaining an initial spurious resolution-independent
blowup. Fixed: the background now holds `u_V` uniform (`u = U0/rho`), with
the actual small perturbation applied to `u_V`.
"""

import numpy as np
from scipy.spatial import cKDTree

GAMMA_3D = 1.936492  # kernel_gamma, Wendland C2, 3D (src/kernel_hydro.h)
ETA = 1.2348  # resolution_eta in every shipped SubgridRadiation example


def wc2_3d_w_dwdr(r, H):
    """dW/dr(r) of the 3D Wendland C2 kernel, support radius H, exactly as
    SWIFT computes it (same convention as verify_isrf_dissipation.py)."""
    r = np.atleast_1d(np.asarray(r, dtype=np.float64))
    H = np.broadcast_to(np.asarray(H, dtype=np.float64), r.shape)
    q = r / H
    norm = 21.0 / (2.0 * np.pi)
    inside = q < 1.0
    dwdr = np.zeros_like(q)
    qi = q[inside]
    Hi = H[inside]
    dwdr[inside] = norm * (-20.0 * qi * (1.0 - qi) ** 3) / Hi**4
    return dwdr


def wi_dr_of(r, h):
    """`wi_dr = h^-(dim+1) dW/dq|_{r/h}`, SWIFT's own convention, dim=3."""
    H = GAMMA_3D * h
    return wc2_3d_w_dwdr(r, H)


def m1_closure_f_chi(u, F, c_hyp):
    """`f = min(1, |F|/(c_M*u))` (0 for `u<=0`) and `chi(f)`, matching
    #radiation_get_m1_closure_tensor_band's own zero-flux-guarded formula."""
    F2 = np.sum(F * F, axis=1)
    F_inv = np.where(F2 > 0.0, 1.0 / np.sqrt(np.maximum(F2, 1e-300)), 0.0)
    Fmag = F2 * F_inv
    denom = c_hyp * u
    f = np.where(denom > 0.0, np.minimum(Fmag / np.maximum(denom, 1e-300), 1.0), 0.0)
    sq = 4.0 - 3.0 * f * f
    chi = (3.0 + 4.0 * f * f) / (5.0 + 2.0 * np.sqrt(sq))
    return f, chi


def m1_closure_tensor(u, F, c_hyp):
    """Vectorized #radiation_get_m1_closure_tensor_band (N particles at
    once). Zero-flux guard matching the source's own convention."""
    F2 = np.sum(F * F, axis=1)
    F_inv = np.where(F2 > 0.0, 1.0 / np.sqrt(np.maximum(F2, 1e-300)), 0.0)
    n = F * F_inv[:, None]
    f, chi = m1_closure_f_chi(u, F, c_hyp)
    iso_coeff = 0.5 * (1.0 - chi)
    aniso_coeff = 0.5 * (3.0 * chi - 1.0)
    eye = np.eye(3)[None, :, :]
    outer_nn = n[:, :, None] * n[:, None, :]
    return iso_coeff[:, None, None] * eye + aniso_coeff[:, None, None] * outer_nn


def build_variable_h_lattice(n_per_dim, box_l, rho0, A_rho, n_wave_rho):
    """Uniform-position lattice (fixed spacing), density field varying
    smoothly along x, an SPH particle's own mass is fixed, so this gives
    a smoothly varying h(x) = ETA*(m/rho(x))^(1/3), matching how a real
    disordered/multi-phase distribution's h field looks, without needing
    actual multi-phase relaxation (pos/rho/mass/h are all FROZEN background
    fields here, exactly what the real per-step rho_prev snapshot is)."""
    dxp = box_l / n_per_dim
    c1 = (np.arange(n_per_dim) + 0.5) * dxp
    xx, yy, zz = np.meshgrid(c1, c1, c1, indexing="ij")
    pos = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1)
    rho = rho0 * (1.0 + A_rho * np.sin(2.0 * np.pi * n_wave_rho * pos[:, 0] / box_l))
    mass = np.full(pos.shape[0], rho0 * dxp**3)
    h = ETA * (mass / rho) ** (1.0 / 3.0)
    return pos, dxp, rho, mass, h


def precompute_pairs(pos, box_l, h):
    """Pair list, minimum-image separations, and BOTH particles' own
    kernel-gradient terms (h_i != h_j in general). pos/rho/mass/h never
    change during the run, so this is computed once and reused every step."""
    h_max = np.max(h)
    tree = cKDTree(pos, boxsize=box_l)
    pairs = tree.query_pairs(r=GAMMA_3D * h_max, output_type="ndarray")
    i_idx, j_idx = pairs[:, 0], pairs[:, 1]
    dxv = pos[i_idx] - pos[j_idx]
    dxv -= box_l * np.round(dxv / box_l)
    rr = np.linalg.norm(dxv, axis=1)
    keep = (rr > 0.0) & (rr < GAMMA_3D * np.maximum(h[i_idx], h[j_idx]))
    i_idx, j_idx, dxv, rr = i_idx[keep], j_idx[keep], dxv[keep], rr[keep]
    r_inv = 1.0 / rr
    wi_dr = wi_dr_of(rr, h[i_idx])
    wj_dr = wi_dr_of(rr, h[j_idx])
    wbar_dr = 0.5 * (wi_dr + wj_dr)
    return dict(
        i=i_idx,
        j=j_idx,
        dx=dxv,
        r_inv=r_inv,
        wi_dr=wi_dr,
        wj_dr=wj_dr,
        wbar_dr=wbar_dr,
    )


def divergence_accumulate_varh(pairs, rho, mass, F, N):
    """`div_specific_flux` accumulator, #radiation_divergence_accumulate_band
    (`diffmode==1`, unchanged by the P1-to-M1 upgrade), generalized to
    h_i != h_j (each particle's own wi_dr/wj_dr, not a shared scalar)."""
    i, j, dxv, r_inv = pairs["i"], pairs["j"], pairs["dx"], pairs["r_inv"]
    wi_dr, wj_dr = pairs["wi_dr"], pairs["wj_dr"]
    Fi_dot = np.einsum("ij,ij->i", F[i], dxv)
    Fj_dot = np.einsum("ij,ij->i", F[j], dxv)
    Phi_ij = Fi_dot / rho[i] * wi_dr * r_inv + Fj_dot / rho[j] * wj_dr * r_inv
    out = np.zeros(N)
    np.add.at(out, i, mass[j] * Phi_ij)
    np.add.at(out, j, -mass[i] * Phi_ij)
    return out


def gradient_accumulate_varh(pairs, rho, mass, u, D, N):
    """`grad_u` accumulator, #radiation_gradient_accumulate_band. Updated to
    `diffmode==0` (each particle's own `wi_dr`/`wj_dr` separately, no shared
    average, no grad-h factor), matching the fix that landed in
    `radiation_propagation_iact.h` (working tree, uncommitted, not edited by
    this script): exact
    skew-adjoint of the `diffmode==1` divergence loop in the `D^-1`-weighted
    inner product (not the plain `m*rho` one), for any `h_i != h_j`, whenever
    `D` is locally constant between neighbours."""
    i, j, dxv, r_inv = pairs["i"], pairs["j"], pairs["dx"], pairs["r_inv"]
    wi_dr, wj_dr = pairs["wi_dr"], pairs["wj_dr"]
    rho_i_inv = 1.0 / rho[i]
    rho_j_inv = 1.0 / rho[j]
    Di_dot_dx = np.einsum("nab,nb->na", D[i], dxv)
    Dj_dot_dx = np.einsum("nab,nb->na", D[j], dxv)
    temp_i = Di_dot_dx * (rho[i] * u[i] * r_inv)[:, None]
    temp_j = Dj_dot_dx * (rho[j] * u[j] * r_inv)[:, None]
    diff = temp_i - temp_j
    fac_i = mass[j] * rho_i_inv**2 * wi_dr
    fac_j = mass[i] * rho_j_inv**2 * wj_dr
    out = np.zeros((N, 3))
    np.add.at(out, i, -diff * fac_i[:, None])
    np.add.at(out, j, -diff * fac_j[:, None])
    return out


def m1_closure_D_inverse(D_tensor, f_val, chi_val):
    """`D^-1` for the closure tensor `D = iso*I + aniso*(n outer n)`, using
    the eigen-decomposition directly (D3/investigation Section 4):
    eigenvalue `chi` along `n`, `(1-chi)/2` transverse (2-fold). Well-defined
    at `F=0` (`n=0`, `chi=1/3`): reduces to `3*I` with no special case.
    Degenerates as `f -> 1` (`chi -> 1`, transverse eigenvalue -> 0); floored
    for numerical safety, but this diagnostic must not be trusted near there
    (investigation Section 4's own bound)."""
    N = D_tensor.shape[0]
    eye = np.eye(3)[None, :, :]
    # Recover n-outer-n from D and the iso/aniso coefficients algebraically
    # is fragile; instead rebuild directly from f_val/chi_val's own
    # definitions, matching m1_closure_tensor's construction.
    chi_safe = np.clip(chi_val, 1.0 / 3.0, 1.0 - 1e-4)
    inv_along = 1.0 / chi_safe
    inv_transverse = 2.0 / np.maximum(1.0 - chi_safe, 1e-4)
    outer_nn = (D_tensor - (0.5 * (1.0 - chi_val))[:, None, None] * eye) / np.maximum(
        (0.5 * (3.0 * chi_val - 1.0))[:, None, None], 1e-30
    )
    # At f=0 the aniso coefficient is 0/0-fragile above; outer_nn is exactly
    # 0 there by construction (F=0), so guard explicitly rather than trust
    # the division.
    outer_nn = np.where((chi_val[:, None, None] > 1.0 / 3.0 + 1e-9), outer_nn, 0.0)
    D_inv = (
        inv_transverse[:, None, None] * eye
        + (inv_along - inv_transverse)[:, None, None] * outer_nn
    )
    return D_inv


def apply_flux_limiter(u, c_hyp, F):
    """Vectorized #radiation_apply_flux_limiter_band. Returns the limited F
    and the pre-clamp ratio |F|/(c_hyp*u) (for monitoring whether the
    limiter ever actually engages, it must not, for this diagnostic to
    isolate the LINEAR/quasi-linear propagation operator alone)."""
    F2 = np.sum(F * F, axis=1)
    Fmag = np.sqrt(np.maximum(F2, 0.0))
    active = u > 0.0
    ratio = np.zeros_like(u)
    ratio[active] = Fmag[active] / np.maximum(c_hyp[active] * u[active], 1e-300)
    lim = np.minimum(1.0, np.where(active, 1.0 / np.maximum(ratio, 1e-300), 0.0))
    F_out = F.copy()
    nonzero_F = F2 > 0.0
    scale = np.where(
        active & nonzero_F, np.minimum(1.0, c_hyp * u / np.maximum(Fmag, 1e-300)), 1.0
    )
    F_out[active & nonzero_F] *= scale[active & nonzero_F][:, None]
    F_out[~active] = 0.0
    return F_out, ratio


def run_sim(
    n_per_dim,
    box_l,
    T_end,
    C_courant,
    c_hyp0,
    rho0,
    A_rho,
    n_wave_rho,
    u0,
    eps_u,
    n_wave_u,
    record_every=1,
):
    pos, dxp, rho, mass, h = build_variable_h_lattice(
        n_per_dim, box_l, rho0, A_rho, n_wave_rho
    )
    N = pos.shape[0]
    pairs = precompute_pairs(pos, box_l, h)
    c_hyp = np.full(N, c_hyp0)

    dt = C_courant * h.min() / c_hyp0
    n_steps = int(np.ceil(T_end / dt))

    # Background holds u_V = rho*u UNIFORM (u = U0/rho), perturbation applied
    # to u_V, not to u directly: the gradient/divergence pair operates on
    # u_V, so perturbing u itself with rho(x) varying underneath it would
    # make u_V inherit the full density-modulation amplitude, not the
    # intended small perturbation (see module docstring, "IC BUG FIXED").
    u_V = u0 * (1.0 + eps_u * np.cos(2.0 * np.pi * n_wave_u * pos[:, 0] / box_l))
    u = u_V / rho
    F = np.zeros((N, 3))

    def energy_iso(u_, F_):
        return float(
            np.sum(
                mass
                * rho
                * (0.5 * u_**2 + 3.0 * np.sum(F_ * F_, axis=1) / (2.0 * c_hyp**2))
            )
        )

    def energy_Dinv(u_, F_, D_):
        Dinv = m1_closure_D_inverse(D_, *m1_closure_f_chi(u_, F_, c_hyp))
        Dinv_F = np.einsum("nab,nb->na", Dinv, F_)
        return float(
            np.sum(
                mass
                * rho
                * (0.5 * u_**2 + np.sum(F_ * Dinv_F, axis=1) / (2.0 * c_hyp**2))
            )
        )

    sum_mu0 = float(np.sum(mass * u))
    D0 = m1_closure_tensor(u, F, c_hyp)
    t_hist = [0.0]
    E_hist = [energy_iso(u, F)]
    QD_hist = [energy_Dinv(u, F, D0)]
    sum_mu_hist = [sum_mu0]
    max_f_ratio = 0.0

    for step in range(1, n_steps + 1):
        div_F = divergence_accumulate_varh(pairs, rho, mass, F, N)
        u_new = u - dt * div_F  # kappa=0, source=0: decay=1, phi=1 exactly

        D = m1_closure_tensor(u_new, F, c_hyp)  # built from NEW u, OLD F,
        # matching the real gradient loop's read order (density ghost has
        # already written this step's u; specific_flux is not updated until
        # the gradient loop's own end-of-loop ghost, below).
        grad_u = gradient_accumulate_varh(pairs, rho, mass, u_new, D, N)
        F_new = F - c_hyp[:, None] ** 2 * dt * grad_u  # kappa=0: decay=1, phi=1

        F_lim, ratio = apply_flux_limiter(u_new, c_hyp, F_new)
        max_f_ratio = max(max_f_ratio, float(np.max(ratio)))

        u, F = u_new, F_lim

        if step % record_every == 0 or step == n_steps:
            D_now = m1_closure_tensor(u, F, c_hyp)
            t_hist.append(step * dt)
            E_hist.append(energy_iso(u, F))
            QD_hist.append(energy_Dinv(u, F, D_now))
            sum_mu_hist.append(float(np.sum(mass * u)))

    return dict(
        t=np.array(t_hist),
        E=np.array(E_hist),
        QD=np.array(QD_hist),
        sum_mu=np.array(sum_mu_hist),
        max_f_ratio=max_f_ratio,
        h_min=float(h.min()),
        h_max=float(h.max()),
        dt=dt,
        n_steps=n_steps,
        N=N,
        sum_mu0=sum_mu0,
        E0=E_hist[0],
        QD0=QD_hist[0],
    )


def _drift_stats(label_metric, t, Q, Q0):
    n = len(t)
    i2 = n // 2
    i1a, i1b = i2, i2 + (n - i2) // 2
    i2a, i2b = i1b, n
    slope_first_half = np.polyfit(t[i1a:i1b], Q[i1a:i1b], 1)[0]
    slope_second_half = np.polyfit(t[i2a:i2b], Q[i2a:i2b], 1)[0]
    slope_full = np.polyfit(t[i2:], Q[i2:], 1)[0]
    print(
        f"    [{label_metric}] Q0 = {Q0:.6e}, Q_final = {Q[-1]:.6e}, "
        f"(Q_final-Q0)/Q0 = {(Q[-1]-Q0)/Q0:+.3e}"
    )
    print(
        f"      dQ/dt (second-half fit) = {slope_full:+.4e} "
        f"(normalized {slope_full/Q0:+.4e}/time); first/last quarter of "
        f"second half = {slope_first_half:+.4e} / {slope_second_half:+.4e} "
        f"(ratio {slope_second_half/slope_first_half if slope_first_half != 0 else float('nan'):+.3f}; "
        f"~1 = bounded accumulation, >>1 = accelerating/real instability)"
    )
    return dict(
        slope_full=slope_full,
        slope_full_norm=slope_full / Q0,
        slope_first_half=slope_first_half,
        slope_second_half=slope_second_half,
    )


def report(label, res):
    t, sum_mu = res["t"], res["sum_mu"]
    sum_mu_drift = np.max(np.abs(sum_mu - res["sum_mu0"])) / max(
        abs(res["sum_mu0"]), 1e-300
    )
    print(
        f"  [{label}] N={res['N']}, dt={res['dt']:.4e}, n_steps={res['n_steps']}, "
        f"h in [{res['h_min']:.4e}, {res['h_max']:.4e}]"
    )
    print(f"    sum(m*u): max drift/|sum_mu0| = {sum_mu_drift:.3e} (must be ~roundoff)")
    stats_iso = _drift_stats("isotropic E", t, res["E"], res["E0"])
    stats_Dinv = _drift_stats("D^-1-weighted Q_D", t, res["QD"], res["QD0"])
    print(
        f"    max |F|/(c_hyp*u) reached (limiter engagement check) = "
        f"{res['max_f_ratio']:.3e} (must stay << 1: limiter must not engage)"
    )
    return dict(
        sum_mu_drift=sum_mu_drift,
        iso=stats_iso,
        Dinv=stats_Dinv,
        max_f_ratio=res["max_f_ratio"],
    )


# ---------------------------------------------------------------------------
# Common parameters
# ---------------------------------------------------------------------------
box_l = 1.0
rho0 = 2.0
c_hyp0 = 1.0
u0 = 1.0
eps_u = 2e-3  # smooth, low-amplitude seed perturbation
n_wave_u = 3
n_wave_rho = 2
C_courant = 0.2  # well inside Part C's admissible alpha_max(C_hyp) envelope
T_end = 20.0  # ~20 box light-crossing times (c_hyp0=1, box_l=1): long enough
# for a slowly-accumulating residual to rise clearly above round-off.

print("=" * 78)
print("Control: UNIFORM h (A_rho=0), D^-1 metric must be identically zero")
print("by construction (mode 0 is the exact skew-adjoint of the divergence")
print("loop, in the D^-1 metric, whenever D is locally constant, includes")
print("F=0/uniform h trivially). This is the self-check that this script's")
print("own re-implementation is faithful to the source, not just a")
print("plausible-looking numpy program. NOTE: the isotropic E metric is NOT")
print("expected to be exactly conserved even here, once F grows away from 0")
print("(D drifts from I/3): only Q_D (D^-1 metric) has the machine-precision")
print("guarantee at fixed, locally-constant D.")
print("=" * 78)
res_uniform = run_sim(
    n_per_dim=10,
    box_l=box_l,
    T_end=T_end,
    C_courant=C_courant,
    c_hyp0=c_hyp0,
    rho0=rho0,
    A_rho=0.0,
    n_wave_rho=n_wave_rho,
    u0=u0,
    eps_u=eps_u,
    n_wave_u=n_wave_u,
)
r_uniform = report("uniform h, n=10", res_uniform)
assert r_uniform["sum_mu_drift"] < 1e-8
assert abs(r_uniform["Dinv"]["slope_full_norm"]) < 1e-6, (
    "uniform-h control shows a non-trivial D^-1-metric drift: the "
    "re-implementation is NOT faithful to the source's exact-adjointness "
    "property; do not trust the variable-h results below until this is fixed."
)
print("  PASS: uniform-h control's D^-1 metric is flat to round-off, as required.")
print()

print("=" * 78)
print("Item 4: VARIABLE h, two resolutions (diffmode==0 gradient loop,")
print("matching the working-tree fix in radiation_propagation_iact.h)")
print("=" * 78)
res_lo = run_sim(
    n_per_dim=10,
    box_l=box_l,
    T_end=T_end,
    C_courant=C_courant,
    c_hyp0=c_hyp0,
    rho0=rho0,
    A_rho=0.5,
    n_wave_rho=n_wave_rho,
    u0=u0,
    eps_u=eps_u,
    n_wave_u=n_wave_u,
)
r_lo = report("variable h, n=10 (coarse)", res_lo)
print()
res_hi = run_sim(
    n_per_dim=16,
    box_l=box_l,
    T_end=T_end,
    C_courant=C_courant,
    c_hyp0=c_hyp0,
    rho0=rho0,
    A_rho=0.5,
    n_wave_rho=n_wave_rho,
    u0=u0,
    eps_u=eps_u,
    n_wave_u=n_wave_u,
)
r_hi = report("variable h, n=16 (fine)", res_hi)
print()

assert r_lo["sum_mu_drift"] < 1e-8
assert r_hi["sum_mu_drift"] < 1e-8
print("  sum(m*u) conservation confirmed at round-off at both resolutions")
print("  (isolates the unchanged, exactly antisymmetric density loop).")
print()

for metric in ("iso", "Dinv"):
    label = "isotropic E" if metric == "iso" else "D^-1-weighted Q_D"
    print(
        f"  Resolution comparison, {label} (h_min: {res_lo['h_min']:.4e} -> "
        f"{res_hi['h_min']:.4e}):"
    )
    lo_norm = abs(r_lo[metric]["slope_full_norm"])
    hi_norm = abs(r_hi[metric]["slope_full_norm"])
    print(f"    |dQ/dt|/Q0: coarse = {lo_norm:.4e}, fine = {hi_norm:.4e}")
    shrinks = hi_norm < lo_norm
    print(f"    {'SHRINKS' if shrinks else 'DOES NOT SHRINK'} with resolution.")
    accel_lo = abs(r_lo[metric]["slope_second_half"]) / max(
        abs(r_lo[metric]["slope_first_half"]), 1e-300
    )
    accel_hi = abs(r_hi[metric]["slope_second_half"]) / max(
        abs(r_hi[metric]["slope_first_half"]), 1e-300
    )
    print(
        f"    In-time acceleration ratio (fixed h): coarse = {accel_lo:.3f}, "
        f"fine = {accel_hi:.3f} (>>1 at either resolution = STOP condition)."
    )
    if shrinks and accel_lo < 3.0 and accel_hi < 3.0:
        print(f"    VERDICT ({label}): drift SHRINKS WITH RESOLUTION and does not")
        print("    accelerate in time at fixed h, expected, bounded consistency")
        print("    error. Monitor only, no action needed.")
    else:
        print(f"    VERDICT ({label}): STOP CONDITION, drift does not shrink with")
        print("    resolution and/or accelerates in time at fixed h. Escalate to")
        print("    the operator; do not proceed into Phase 2 on this script's own")
        print("    authority.")
    print()

print("Reminder (investigation Section 4's own bound): the D^-1 metric")
print("degenerates as f -> 1 (transverse eigenvalue of D -> 0), which the")
print("flux limiter explicitly permits. This diagnostic's clean-Q_D result,")
print("if obtained, is evidence for the near-P1/background regime (small f")
print("throughout this run, see the max |F|/(c_hyp*u) figures above) and")
print("must NOT be read as an unconditional stability proof at f -> 1: per")
print("the investigation, that regime's stability must rest on the")
print("dissipation stages and the CFL bound on c_M, not on adjointness.")
print()
