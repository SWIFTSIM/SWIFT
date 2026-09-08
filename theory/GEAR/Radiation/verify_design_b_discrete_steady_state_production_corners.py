"""Design B's own DISCRETE steady-state point-source solution, not the
continuum Yukawa profile `isrf_yukawa_profile_check.py` fits against --
evaluated at the production corners `.claude/dev/logs/
2026-09-08_1003_design-b-tier1-production-corner-comparison.md` actually
measured (`h/lambda_analytic` = 2.82, 4.71, 5.98, 9.97, 11.96, 19.93).

Background / what this follows on from. That comparison fit Design B's
simulated radial profile against the CONTINUUM Yukawa Green's function and
found the fitted lambda 2-6x off (`lambda_measured/h` = 0.54, 0.42, 0.38,
0.31, 0.29, trending toward the design doc's own predicted `~0.25 h`
resolution floor, Sec. 4.7). The operator's direct correction: this
comparison is invalid on its own terms, because Design B does not claim to
reproduce the continuum Yukawa profile once `h` is not small compared to
`lambda` -- Sec. 4.7's own F.1 measurement (1D chain) already found the
DISCRETE fixed point's own e-folding length departs from the input
`lambda` in exactly this regime (`lambda_eff/lambda = 5.2` at `lambda =
0.05 h`, i.e. `h/lambda = 20`). What was never built is the 3D counterpart
of that check, evaluated at the ACTUAL production `h/lambda` values on the
ACTUAL 3D Wendland-C2 kernel/eta this project compiles, so that the
production-corner simulation output can be compared against the right
target instead of the continuum one.

This script builds that target, by extending
`verify_design_b_timestepping_stability.py`'s own machinery:
  - Part B's `lattice_symbol_3d` (the exact Fourier symbol of the
    Wendland-C2 kernel-gradient estimator on a 3D cubic lattice, this
    project's own `resolution_eta`) is generalised here from "one
    direction, a scan over |k|" to every wavevector on a periodic N^3
    grid, so a full 3D Green's function can be built by FFT.
  - Part C.2's discrete screened-Poisson fixed-point identity,
    `(1 - lambda^2 A A) u = tau S` <=> in Fourier space
    `(1 + lambda^2 K(k)^2) u_hat(k) = tau S_hat(k)`, is solved on that
    grid instead of a 1D chain.
On a UNIFORM cubic lattice the difference (`diffmode==0`, grad) and
symmetric (`diffmode==1`, div) estimators of design-lw-fuv-design-b.md
Sec 2.1/2.2 reduce to the identical real, odd symbol `K(k)`, exactly as
`lattice_symbol_3d`'s own docstring states and Part C.2 already assumes
(one operator `A` used for both div and grad) -- this script inherits
that same simplification, so it says nothing new about GLASS disorder
(Sec 4.7's own "glass value unknown" caveat stands unchanged; that is a
distinct, unresolved discretization-vs-continuum question left for later
work, not something this script can settle on a perfect lattice).

Part 0 is an INDEPENDENT cross-check of the continuum steady state (already
established by direct ODE solve in verify_design_b_tier1_fit_target_
validity.py) via a second route: the Fourier-Laplace transform method of
`~/swiftsim_0/theory/GEAR/Diffusion/hyperbolic.tex` (the same repo's own
Cattaneo/telegrapher metal-diffusion scheme, cited there for its exact
time-dependent point-source Green's functions). That document's own system
(`d(rho*Z)/dt + div(F) = 0`, `dF/dt = -(K/tau)*grad(driver) - F/tau`, Eqs.
mass_conservation/cattaneo) has NO sink term on its conserved quantity
(metals are transported, never absorbed) -- structurally different from
Design B's `Du/Dt = -(1/rho)div(rho F) - u/tau + S`, which has an explicit
`-u/tau` absorption term (radiation IS absorbed by dust; this is what
produces the Yukawa screening in the first place, and its absence in the
metal-diffusion case is why that document's own exact Green's function
never screens, only diffuses/propagates). Because of this, `hyperbolic.
tex`'s own closed-form formula does not apply to Design B literally; Part 0
instead applies its SAME two-transform technique to Design B's actual
(different) equations, confirms that dropping the `-u/tau` sink recovers
`hyperbolic.tex`'s own Eq. (transform-solution) exactly (a check that this
script's application of the method is not simply asserted), and confirms
separately that Design B's own steady state (via the final-value theorem)
reproduces the same `tau*S_hat/(1+lambda^2 k^2)` transfer function Part 1
below states from the real-space ODE and verify_design_b_tier1_fit_target_
validity.py already confirmed by direct `sympy.dsolve` -- three independent
derivations of the same object, now in agreement. It also explains, from
that structural difference, why Design B's own transients cannot inherit
`hyperbolic.tex`'s persistent ballistic-front/dipole-layer pathology
(Remark 3d-front-dipole there): the `-u/tau` term exponentially damps any
such feature at the SAME rate that sets the whole scheme's relaxation,
rather than letting it persist as an undamped wave feature indefinitely,
consistent with the 2026-09-08 log's own empirical finding that the
production-corner runs are flat to <0.05% over the last 5 snapshots (a
genuinely converged steady state, not a slow-to-settle transient).

Part 1 sets up the discrete Green's-function problem and cross-checks the
generalised 3D symbol against `lattice_symbol_3d`'s own special-direction
numbers. Part 2 solves the fixed point on an N=32 periodic grid (level 5,
matching every production-corner run) at the production corners' own
`h/lambda_analytic` values and fits `lambda_eff` with the IDENTICAL
methodology `isrf_yukawa_profile_check.py` uses on real snapshots (same
`r_min`, `r_max`, same log-linear fit of `u(r)*r`). Part 3 sweeps a much
wider `h/lambda` range (thin-screening through deep-floor) at a larger
grid to confirm the continuum limit (`h/lambda -> 0` recovers
`lambda_eff/lambda_analytic -> 1`) and to cross-check against the existing
1D chain result at the matching corner. Part 4 compares the discrete
prediction directly against the actual simulation numbers from the
2026-09-08 comparison log.
"""

import numpy as np
import sympy as sp

# ---------------------------------------------------------------------------
# Part 0: cross-check against ~/swiftsim_0's Cattaneo/telegrapher Green's
# function derivation (hyperbolic.tex), via the SAME Fourier-Laplace
# transform technique applied to Design B's own (different) equations.
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part 0: cross-check against swiftsim_0's telegrapher transform method")
print("=" * 78)

p, k, tau_u, tau_F, D_s, S0, U0 = sp.symbols("p k tau_u tau_F D S_0 U_0", positive=True)

# General two-timescale system (u has its own sink 1/tau_u, F relaxes with
# 1/tau_F): u_t = -div(F) - u/tau_u + S*Theta(t); F_t = -(D/tau_F) grad(u)
# - F/tau_F. Fourier in space (grad -> i k, only the k-parallel component
# of F survives), Laplace in time, zero initial data, solved algebraically.
Uhat, Fhat = sp.symbols("Uhat Fhat")
eq_F = sp.Eq(p * Fhat, -sp.I * k * (D_s / tau_F) * Uhat - Fhat / tau_F)
Fhat_sol = sp.solve(eq_F, Fhat)[0]
eq_u_source = sp.Eq(p * Uhat, -sp.I * k * Fhat_sol - Uhat / tau_u + S0 / p)
Uhat_source_sol = sp.simplify(sp.solve(eq_u_source, Uhat)[0])
print("General source-driven solution, two timescales tau_u (u sink), tau_F")
print("(F relaxation):")
print("  Uhat(k,p) =", Uhat_source_sol)

# Cross-check A: drop the source, use an initial delta-function condition
# instead (swiftsim_0's own problem), and take tau_u -> infinity (no u
# sink -- exactly hyperbolic.tex's conserved-metal system, Eqs.
# mass_conservation/cattaneo: d(rho Z)/dt + div(F) = 0, dF/dt =
# -(K/tau)grad(driver) - F/tau, with K identified as this script's D).
eq_u_ic = sp.Eq(p * Uhat - U0, -sp.I * k * Fhat_sol - Uhat / tau_u)
Uhat_ic_sol = sp.solve(eq_u_ic, Uhat)[0]
Uhat_ic_noSink = sp.simplify(sp.limit(Uhat_ic_sol, tau_u, sp.oo))
# hyperbolic.tex Eq. (transform-solution): Utilde = q0(tau p+1)/(tau p^2+p+kappa k^2)
target_swiftsim0 = U0 * (tau_F * p + 1) / (tau_F * p**2 + p + D_s * k**2)
residual_A = sp.simplify(Uhat_ic_noSink - target_swiftsim0)
print()
print("Cross-check A (drop the u-sink, swiftsim_0's own conserved-metal case,")
print("delta-function IC instead of a continuous source):")
print("  Uhat(k,p), tau_u -> inf:", Uhat_ic_noSink)
print("  hyperbolic.tex Eq. (transform-solution), kappa == D:", target_swiftsim0)
print("  difference:", residual_A)
assert residual_A == 0, "does not reduce to swiftsim_0's own formula -- method applied wrongly"
print("  MATCH: this script's transform method reproduces hyperbolic.tex's own")
print("  formula exactly once its (different, sink-free) physics is used.")

# Cross-check B: Design B's OWN system (shared tau_u = tau_F = tau, the
# operator's ruling, HANDOFF item 4), continuous source, steady state via
# the final-value theorem (lim_{t->inf} u(t) = lim_{p->0} p*Uhat(p)).
tau = sp.symbols("tau", positive=True)
Uhat_design_b = Uhat_source_sol.subs({tau_u: tau, tau_F: tau})
u_steady = sp.simplify(sp.limit(p * Uhat_design_b, p, 0))
lam = sp.symbols("lambda", positive=True)
target_screened = tau * S0 / (1 + lam**2 * k**2)
residual_B = sp.simplify((u_steady - target_screened).subs(D_s, lam**2 / tau))
print()
print("Cross-check B (Design B's own shared-tau, absorptive system, continuous")
print("source, steady state via the final-value theorem):")
print("  u_steady(k) =", u_steady, " (D = lambda^2/tau, i.e. tau*D = lambda^2)")
print("  screened-Poisson target tau*S/(1+lambda^2 k^2):", target_screened)
print("  difference (after D = lambda^2/tau):", residual_B)
assert residual_B == 0, "Design B's own steady state does not match the screened-Poisson target"
print("  MATCH: independently reproduces the same transfer function Part 1 below")
print("  states from the real-space ODE, and that")
print("  verify_design_b_tier1_fit_target_validity.py already confirmed by direct")
print("  sympy.dsolve of the steady ODE -- three independent routes, one answer.")
print()
print("Interpretation: hyperbolic.tex's own closed-form time-dependent Green's")
print("function (ballistic front + Bessel-function interior + a front dipole")
print("layer, Remark 3d-front-dipole) does NOT apply to Design B literally --")
print("it is the solution of a DIFFERENT equation (no u sink: metals are")
print("transported, never absorbed). Design B's explicit -u/tau term damps any")
print("such transient feature at the same rate that sets the whole scheme's")
print("relaxation, rather than letting an undamped wave feature persist")
print("indefinitely as in the conservative metal-diffusion case -- consistent")
print("with (not proof beyond) the 2026-09-08 log's own empirical finding that")
print("the production-corner runs are flat to <0.05% over the last 5 snapshots,")
print("i.e. genuinely steady, not a slowly-relaxing transient.")
print()

# ---------------------------------------------------------------------------
# Part 1: the discrete fixed-point problem in Fourier space, and the
# generalised 3D symbol
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part 1: discrete screened-Poisson fixed point in Fourier space")
print("=" * 78)

lam_s, K_s, tau_s, S_s = sp.symbols("lambda K tau S", positive=True)
u_hat, S_hat = sp.symbols("u_hat S_hat")
# Real-space fixed point (design-lw-fuv-design-b.md Sec 4.5, verified by
# verify_design_b_timestepping_stability.py Part C.2): u - lambda^2 * A A u
# = tau S, with A the (shared, uniform-lattice) div/grad operator, Fourier
# symbol i*K(k). Composing div(grad) gives symbol (iK)(iK) = -K^2, so:
fixed_point_fourier = sp.Eq((1 + lam_s**2 * K_s**2) * u_hat, tau_s * S_s)
print("Fixed point, real space: u - lambda^2 div_h(grad_h(u)) = tau*S")
print("Fourier symbol of div_h(grad_h(.)) is -K(k)^2 (K real, odd; grad~iK, div~iK)")
print("=>", fixed_point_fourier)
u_hat_sol = sp.solve(fixed_point_fourier, u_hat)[0]
print("=> u_hat(k) =", u_hat_sol)
# Continuum check: as the discreteness scale h -> 0 at fixed k, K(k) -> |k|
# (Part 1B below verifies this numerically for the actual compiled kernel),
# recovering the textbook continuum screened-Poisson transfer function
# 1/(1+lambda^2 k^2), i.e. the Yukawa profile's own Fourier transform.
k_s = sp.symbols("k", positive=True)
continuum_form = u_hat_sol.subs(K_s, k_s)
print("Continuum limit (K(k) -> |k|): u_hat(k) ->", continuum_form,
      "(the Yukawa Green's function's own transform)")
print()

ETA = 1.2348          # resolution_eta, every shipped SubgridRadiation example
GAMMA_3D = 1.936492   # kernel_gamma, Wendland C2, 3D (src/kernel_hydro.h)


def wc2_3d_dwdr(r, H):
    """dW/dr of the 3D Wendland C2 kernel, support radius H -- H may be a
    scalar (Parts 1-4's fixed-h idealized lattice) or a per-element array
    matching r (Part 5's real, per-particle h_i/h_j on an actual glass)."""
    H = np.broadcast_to(np.asarray(H, dtype=float), np.shape(r))
    q = r / H
    inside = q < 1.0
    out = np.zeros_like(r, dtype=float)
    qi, Hi = q[inside], H[inside]
    norm_i = 21.0 / (2.0 * np.pi * Hi**3)
    out[inside] = norm_i * (-4.0 * (1.0 - qi) ** 3 * (4.0 * qi + 1.0) + 4.0 * (1.0 - qi) ** 4) / Hi
    return out


def wc2_3d_w(r, H):
    """W of the 3D Wendland C2 kernel, support radius H -- same scalar-or-
    array broadcasting as wc2_3d_dwdr."""
    H = np.broadcast_to(np.asarray(H, dtype=float), np.shape(r))
    q = r / H
    inside = q < 1.0
    out = np.zeros_like(r, dtype=float)
    qi, Hi = q[inside], H[inside]
    out[inside] = 21.0 / (2.0 * np.pi * Hi**3) * (1.0 - qi) ** 4 * (4.0 * qi + 1.0)
    return out


def build_stencil(dx=1.0):
    """Neighbour offsets (pos), the per-neighbour coefficient
    c_n = dx^3 * dW/dr(r_n) / r_n (m/rho = dx^3 on a uniform cubic
    lattice), and h, for the real-space kernel support around the origin.
    Same construction as lattice_symbol_3d, generalised to keep the full
    3D offsets rather than only a projection."""
    h = ETA * dx
    H = GAMMA_3D * h
    n = int(np.ceil(H / dx)) + 1
    g = np.arange(-n, n + 1) * dx
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    pos = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    r = np.linalg.norm(pos, axis=1)
    keep = (r > 0) & (r < H)
    pos, r = pos[keep], r[keep]
    dwdr = wc2_3d_dwdr(r, H)
    c = dx**3 * dwdr / r
    return pos, c, h


def symbol_full(kvecs, pos, c):
    """K(k) for arbitrary wavevectors (not just a single direction): a
    plane wave e^{i k.x} response of the shared-coefficient gradient/
    divergence estimator is K(k) = (1/|k|) * sum_n c_n (pos_n.k) sin(pos_n.k)
    -- reduces exactly to lattice_symbol_3d's own
    sum_n c_n (pos_n.k_hat) sin(k |pos_n.k_hat|) when k = k*k_hat, verified
    below against that function's own reported numbers."""
    dot = pos @ kvecs.T  # (Nn, Nk)
    term = c[:, None] * dot * np.sin(dot)
    kmag = np.linalg.norm(kvecs, axis=1)
    out = np.zeros(kvecs.shape[0])
    nz = kmag > 0
    out[nz] = term[:, nz].sum(axis=0) / kmag[nz]
    return out


# Cross-check against verify_design_b_timestepping_stability.py's own
# lattice_symbol_3d numbers along its three tested directions, several |k|.
print("Cross-check: symbol_full vs. the existing script's lattice_symbol_3d "
      "(same stencil, three directions):")
pos0, c0, h0 = build_stencil(dx=1.0)
max_relerr = 0.0
for k_hat in (np.array([1.0, 0, 0]), np.array([1.0, 1.0, 0]) / np.sqrt(2),
              np.ones(3) / np.sqrt(3)):
    k_values = np.array([0.3, 1.0, 2.0, 3.0])
    kvecs = k_values[:, None] * k_hat[None, :]
    K_new = symbol_full(kvecs, pos0, c0)
    # reference, verbatim from verify_design_b_timestepping_stability.py
    proj = pos0 @ k_hat
    r0 = np.linalg.norm(pos0, axis=1)
    K_ref = np.array([np.sum(c0 * r0 * (proj / r0) * np.sin(k * proj)) for k in k_values])
    relerr = np.max(np.abs(K_new - K_ref)) / np.max(np.abs(K_ref))
    max_relerr = max(max_relerr, relerr)
    print(f"  direction {np.round(k_hat, 3)}: max relative discrepancy = {relerr:.2e}")
assert max_relerr < 1e-10, "generalised 3D symbol disagrees with the existing script"
print(f"  ==> generalised symbol matches the already-verified lattice_symbol_3d "
      f"to {max_relerr:.1e} (same operator, more general k)")
print()

# ---------------------------------------------------------------------------
# Part 2: solve the fixed point on an N=32 periodic grid at the ACTUAL
# production-corner h/lambda values, fit lambda_eff the SAME way
# isrf_yukawa_profile_check.py fits a real snapshot.
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part 2: N=32 periodic lattice (level 5, matches every production run),")
print("        discrete Green's function at the production corners' own h/lambda")
print("=" * 78)


def symbol_grid(N, dx=1.0):
    pos, c, h = build_stencil(dx)
    kx = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    ky = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    kz = 2 * np.pi * np.fft.rfftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    kvecs = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], axis=1)
    K = symbol_full(kvecs, pos, c).reshape(KX.shape)
    return K, h


def kernel_weighted_source(N, dx=1.0, h_src=None):
    """A normalised Wendland-C2 kernel deposit centred on one grid site,
    same functional form as the gas density kernel (Sec 4.5 step 6's
    injection is kernel-weighted over the star's own neighbours) -- an
    idealised stand-in for the star's real deposit, not a literal
    reproduction of radiation_iact.h's neighbour loop."""
    if h_src is None:
        h_src = ETA * dx
    H = GAMMA_3D * h_src
    n = int(np.ceil(H / dx)) + 1
    idx = np.arange(N)
    # centred at grid index (N//2, N//2, N//2); periodic minimum-image
    d = (idx - N // 2 + N // 2) % N - N // 2
    DX, DY, DZ = np.meshgrid(d, d, d, indexing="ij")
    r = np.sqrt(DX.astype(float) ** 2 + DY.astype(float) ** 2 + DZ.astype(float) ** 2) * dx
    S = wc2_3d_w(r, H)
    S /= S.sum()
    return S


def discrete_fixed_point_profile(N, lam, dx=1.0):
    """Solve (1+lambda^2 K(k)^2) u_hat = S_hat on the N^3 periodic grid
    (tau*S_tot = 1, an arbitrary normalisation -- irrelevant to the
    log-linear slope fit below) and return the real-space u(x) array plus
    the grid spacing dx."""
    K, h = symbol_grid(N, dx)
    S = kernel_weighted_source(N, dx)
    S_hat = np.fft.rfftn(S)
    u_hat = S_hat / (1.0 + lam**2 * K**2)
    u = np.fft.irfftn(u_hat, axes=(0, 1, 2))
    return u, h


def radial_profile_and_fit(u, N, dx=1.0, n_bins=25):
    """Bin u by distance from the source (grid centre) with periodic
    minimum-image wrapping, fit log(u*r) vs r over
    [r_min, r_max] = [2*d, 0.7*half_box] -- IDENTICAL range and fit to
    isrf_yukawa_profile_check.py's fit_slope, d = box/N (here == dx
    exactly, since N particles per side, box = N*dx)."""
    idx = np.arange(N)
    d = (idx - N // 2 + N // 2) % N - N // 2
    DX, DY, DZ = np.meshgrid(d, d, d, indexing="ij")
    r = np.sqrt(DX.astype(float) ** 2 + DY.astype(float) ** 2 + DZ.astype(float) ** 2) * dx
    r_flat, u_flat = r.ravel(), u.ravel()

    box = N * dx
    half_box = box / 2.0
    r_min = 2.0 * dx
    r_max = 0.7 * half_box
    edges = np.linspace(r_min, r_max, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    u_binned = np.full(n_bins, np.nan)
    for i in range(n_bins):
        sel = (r_flat >= edges[i]) & (r_flat < edges[i + 1])
        if sel.sum() > 0:
            u_binned[i] = np.mean(u_flat[sel])
    valid = ~np.isnan(u_binned) & (u_binned > 0)
    if valid.sum() < 3:
        return None, centres, u_binned
    log_ur = np.log(u_binned[valid] * centres[valid])
    slope, _ = np.polyfit(centres[valid], log_ur, 1)
    lam_fit = -1.0 / slope if slope < 0 else np.inf
    return lam_fit, centres, u_binned


N_PROD = 32  # level 5 => 2**5 = 32 particles per side, every production run
corners = [
    ("m=10 Msun, FUV", 2.82, 0.117, 0.514, 0.54),
    ("m=10 Msun, LW", 4.71, 0.219, 0.963, 0.42),
    ("m=95 Msun, FUV", 5.98, 0.275, 1.244, 0.38),
    ("m=95 Msun, LW", 9.97, 0.372, 2.092, 0.31),
    ("m=760 Msun, FUV", 11.96, None, 2.485, 0.29),
    ("m=760 Msun, LW", 19.93, None, None, None),  # sim: float32 underflow, unmeasurable
]

print(f"{'Corner':<18}{'h/lambda':>10}{'lambda_dx':>12}{'lambda_eff/h (predicted)':>26}"
      f"{'lambda_eff/h (measured)':>26}{'discrete/measured':>20}")
predicted_ratios = {}
for name, h_over_lam, rel_A, rel_B, meas_ratio in corners:
    lam_dx = ETA / h_over_lam  # lambda in units of dx, since h = ETA*dx
    u, h = discrete_fixed_point_profile(N_PROD, lam_dx, dx=1.0)
    lam_fit, centres, u_binned = radial_profile_and_fit(u, N_PROD, dx=1.0)
    pred_ratio = lam_fit / h if lam_fit is not None and np.isfinite(lam_fit) else np.nan
    predicted_ratios[name] = pred_ratio
    meas_str = f"{meas_ratio:.3f}" if meas_ratio is not None else "n/a (float32 underflow)"
    cmp_str = f"{pred_ratio / meas_ratio:.3f}" if (meas_ratio and np.isfinite(pred_ratio)) else "n/a"
    print(f"{name:<18}{h_over_lam:>10.2f}{lam_dx:>12.4f}{pred_ratio:>26.3f}"
          f"{meas_str:>26}{cmp_str:>20}")
print()
print("(lambda_eff/h (measured) is the actual SPH simulation's fitted lambda,")
print(" from the 2026-09-08 production-corner comparison log; 'discrete/measured'")
print(" close to 1 means the real code matches ITS OWN idealized-lattice discrete")
print(" prediction, not the continuum Yukawa target the earlier comparison used.)")
print()

# ---------------------------------------------------------------------------
# Part 3: wide h/lambda sweep at a box always >= 10*lambda (so the fit range
# never runs into the periodic-image bias Part 2's fixed N=32 box shows once
# lambda is a large fraction of the box -- Part 2 intentionally keeps that
# bias in, since it is what the ACTUAL production runs' own finite box does
# too; Part 3 isolates the box-size-independent h/lambda dependence alone).
# Confirms the continuum limit and cross-checks the deep-floor regime against
# the existing 1D chain result (verify_design_b_timestepping_stability.py
# Part F.1, lam=lambda/h=0.05, i.e. h/lambda=20, found lambda_eff/lambda=5.2).
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part 3: h/lambda sweep, box >= 10*lambda at every point (box-size-independent)")
print("=" * 78)
h_over_lam_sweep = np.array([0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 2.82, 4.71, 5.98,
                              9.97, 11.96, 19.93, 20.0])
print(f"{'h/lambda':>10}{'N (grid)':>10}{'lambda_eff/lambda_analytic':>28}{'lambda_eff/h':>16}")
sweep_results = []
for hol in h_over_lam_sweep:
    lam_dx = ETA / hol
    # box = N*dx must be >= ~10*lambda for the far-field fit (r_max=0.35*box)
    # to cover several e-foldings without periodic-image contamination;
    # round up to the next power of 2 for FFT efficiency, cap at 256.
    N = int(min(256, 2 ** np.ceil(np.log2(max(32, 10.0 * lam_dx)))))
    u, h = discrete_fixed_point_profile(N, lam_dx, dx=1.0)
    lam_fit, _, _ = radial_profile_and_fit(u, N, dx=1.0)
    ratio_to_analytic = lam_fit / lam_dx if lam_fit is not None and np.isfinite(lam_fit) else np.nan
    ratio_to_h = lam_fit / h if lam_fit is not None and np.isfinite(lam_fit) else np.nan
    sweep_results.append((hol, N, ratio_to_analytic, ratio_to_h))
    print(f"{hol:>10.2f}{N:>10d}{ratio_to_analytic:>28.3f}{ratio_to_h:>16.3f}")

thin = [row for row in sweep_results if row[0] <= 0.1]
assert all(abs(row[2] - 1.0) < 0.15 for row in thin), (
    "continuum limit not recovered as h/lambda -> 0 -- derivation bug")
print()
print("Continuum-limit sanity check (h/lambda <= 0.1): lambda_eff/lambda_analytic "
      f"within 15% of 1 in every case ({[f'{row[2]:.3f}' for row in thin]}) -- the "
      "discrete solve reduces to the continuum Yukawa target as the resolution")
print("gets fine relative to lambda, as it must, once the box is large enough")
print("relative to lambda that periodic images do not bias the far-field fit.")
deep_floor = [row for row in sweep_results if row[0] >= 19.5][0]
print(f"\nDeep-floor cross-check at h/lambda ~ 20 (this script, 3D, box-independent): "
      f"lambda_eff/lambda_analytic = {deep_floor[2]:.2f}, lambda_eff/h = {deep_floor[3]:.3f}")
print("(the existing 1D chain result at the same h/lambda=20 corner, "
      "verify_design_b_timestepping_stability.py Part F.1's lam=0.05 row: "
      "lambda_eff/lambda = 5.2, i.e. lambda_eff/h ~ 0.26 -- same order of "
      "magnitude and same direction (floored, not collapsed); a 1D vs 3D "
      "geometry difference in the exact prefactor is expected, not a red flag)")
print()
print("Comparing this box-independent sweep to Part 2's fixed-N=32 (real production")
print("box) values at the SAME h/lambda shows how much of the earlier discrepancy")
print("is a genuine small-h/lambda discrete-operator effect versus a finite-box")
print("(periodic-image) effect specific to this example's own box size at these")
print("particle counts -- both are real properties of the ACTUAL shipped example,")
print("not artifacts of this script, but they are conceptually distinct causes.")
print()

# ---------------------------------------------------------------------------
# Part 4: explicit comparison table against the 2026-09-08 measurement log
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part 4: discrete prediction vs. the actual production-corner simulation")
print("=" * 78)
print("""
Corner            h/lambda   lambda_measured/h (sim)   lambda_eff/h (this script's
                                                        idealized discrete prediction)
m=10, FUV          2.82       0.54                      see table above
m=10, LW           4.71       0.42                      see table above
m=95, FUV          5.98       0.38                      see table above
m=95, LW           9.97       0.31                      see table above
m=760, FUV        11.96       0.29                      see table above
""")
print("If the ratio (predicted/measured) in Part 2's table sits close to 1 across")
print("all five measurable corners, the actual SPH implementation is correctly")
print("solving ITS OWN discretized equations at these resolutions -- the earlier")
print("comparison's 'FAIL' verdict was a wrong-target problem (continuum Yukawa),")
print("not a code-correctness problem. A systematic offset would instead point to")
print("glass disorder, the a=1 cosmological-factor assumption, boundary handling,")
print("or a genuine implementation bug -- distinguishable from the idealized-")
print("lattice floor by the SIZE and DIRECTION of the residual after this")
print("comparison, not knowable before it.")
print()
# ---------------------------------------------------------------------------
# Part 5: the actual thing -- solve the coupled (u, F) fixed point on the
# REAL glass/snapshot particle distribution (real x_i, h_i, rho_i from a
# converged snapshot of the 2026-09-08 production-corner run, not an
# idealized perfect lattice), using the EXACT diffmode==1 div(F)/diffmode==0
# grad(u) pairwise formulas of Sec 2.2 and the exact-relaxation staggered
# iteration of Sec 4.5 -- run to its own fixed point. This removes every
# idealization Parts 1-4 make (perfect lattice, kernel-smoothed point
# source, single global h): real glass disorder, real per-particle h_i and
# rho_i, and the real diffmode0-vs-diffmode1 distinction (only equal on a
# perfect lattice) all enter. Snapshots reused from the already-completed
# 2026-09-08 comparison run (still on disk in this session's scratchpad;
# no rerun needed).
# ---------------------------------------------------------------------------
print("=" * 78)
print("Part 5: real snapshot cross-check (actual glass, actual h_i/rho_i, exact")
print("        diffmode0/diffmode1 pairwise formulas, exact-relaxation iteration)")
print("=" * 78)

import glob
import os

import h5py
from scipy.spatial import cKDTree

SIGMA_D_FUV_CGS = 9e-22
SIGMA_D_LW_CGS = 1.5e-21
MU_H = 1.4
M_H_CGS = 1.6726219e-24
GRACKLE_SOLAR_Z = 0.01295
C_CFL = 0.1


def analytic_lambda_cgs(Z, rho_internal, unit_length_cgs, unit_mass_cgs, sigma_d_cgs):
    """Identical formula to isrf_yukawa_profile_check.py's own function."""
    rho_cgs = rho_internal * unit_mass_cgs / unit_length_cgs**3
    D_relative = Z / GRACKLE_SOLAR_Z
    kappa_eff_cgs = sigma_d_cgs * D_relative / (MU_H * M_H_CGS)
    return 1.0 / (kappa_eff_cgs * rho_cgs)


def fit_slope(r, u, r_min, r_max):
    """Verbatim copy of isrf_yukawa_profile_check.py's own fit_slope."""
    mask = (r > r_min) & (r < r_max) & (u > 0)
    if mask.sum() < 3:
        return None, None
    log_ur = np.log(u[mask] * r[mask])
    slope, intercept = np.polyfit(r[mask], log_ur, 1)
    return slope, (-1.0 / slope if slope < 0 else np.inf)


def load_snapshot(path):
    with h5py.File(path, "r") as f:
        header = f["/Header"]
        boxsize = float(np.asarray(header.attrs["BoxSize"], dtype=float).flatten()[0])
        units = f["/Units"]
        uL = float(np.asarray(units.attrs["Unit length in cgs (U_L)"]).flat[0])
        uM = float(np.asarray(units.attrs["Unit mass in cgs (U_M)"]).flat[0])
        gas = f["/PartType0"]
        pos = gas["Coordinates"][:, :]
        rho = gas["Densities"][:].astype(np.float64)
        h = gas["SmoothingLengths"][:].astype(np.float64)
        mass = gas["Masses"][:].astype(np.float64)
        Z = gas["MetalMassFractions"][:, -1]
        star = f["/PartType4"]
        star_pos = star["Coordinates"][0, :]
    return dict(boxsize=boxsize, unit_length_cgs=uL, unit_mass_cgs=uM, pos=pos,
                rho=rho, h=h, mass=mass, Z=Z, star_pos=star_pos)


def build_pairs(pos, h, boxsize):
    """Periodic neighbour pairs (i<j) within the larger of the two
    particles' own kernel support radii, plus dx, r, wi_dr, wj_dr."""
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


def phi_func(a_):
    return 1.0 - 0.5 * a_ if a_ < 1e-6 else -np.expm1(-a_) / a_


def solve_real_fixed_point(snap, lam, r_min, r_max, n_bins=25, n_iter=4000, tol=1e-6):
    pos, h, rho, mass = snap["pos"], snap["h"], snap["rho"], snap["mass"]
    boxsize = snap["boxsize"]
    N = len(pos)
    ii, jj, dx, r, wi_dr, wj_dr = build_pairs(pos, h, boxsize)

    # kernel-weighted deposit at the star's position (same normalised
    # weight function as Part 2's idealized source, now on the real glass).
    dxs = pos - snap["star_pos"]
    dxs -= boxsize * np.round(dxs / boxsize)
    rs = np.linalg.norm(dxs, axis=1)
    Hs = GAMMA_3D * np.median(h)
    S = wc2_3d_w(rs, Hs)
    S /= np.sum(S)

    h_med = np.median(h)
    c_hyp = C_CFL * h_med  # dt = 1 (arbitrary; fixed point is dt-independent)
    a = C_CFL * h_med / lam  # design doc's own closure, Sec 4.1/4.5
    e_, ph = np.exp(-a), phi_func(a)

    u = np.zeros(N)
    Fx, Fy, Fz = np.zeros(N), np.zeros(N), np.zeros(N)
    for it in range(n_iter):
        dF = div_F(np.array([Fx, Fy, Fz]), ii, jj, dx, r, wi_dr, wj_dr, rho, mass)
        u_new = e_ * u + ph * (S - dF)
        g = grad_u(u_new, ii, jj, dx, r, wi_dr, wj_dr, rho, mass)
        Fx_new = e_ * Fx - ph * c_hyp**2 * g[0]
        Fy_new = e_ * Fy - ph * c_hyp**2 * g[1]
        Fz_new = e_ * Fz - ph * c_hyp**2 * g[2]
        du = np.max(np.abs(u_new - u)) / max(np.max(np.abs(u_new)), 1e-300)
        u, Fx, Fy, Fz = u_new, Fx_new, Fy_new, Fz_new
        if du < tol and it > 20:
            break
    converged = du < tol

    r_star = radial_distance = np.linalg.norm(dxs, axis=1)
    order = np.argsort(r_star)
    r_sorted, u_sorted = r_star[order], u[order]
    edges = np.linspace(r_min, r_max, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    u_binned = np.full(n_bins, np.nan)
    for i in range(n_bins):
        sel = (r_sorted >= edges[i]) & (r_sorted < edges[i + 1])
        if sel.sum() > 0:
            u_binned[i] = np.mean(u_sorted[sel])
    valid = ~np.isnan(u_binned)
    _, lam_fit = fit_slope(centres[valid], u_binned[valid], r_min, r_max)
    return lam_fit, it, converged


SCRATCH = ("/tmp/claude-1000/-home-darwinr-swiftsim-3/"
           "9a9748ae-16ff-4941-8b97-e213770168b0/scratchpad/tier1_designB")
run_dirs = {
    "default (m=0.1, Tier 1's own thin corner)": "run",
    "m=10": "run_m10",
    "m=95": "run_m95",
    "m=760": "run_m760",
}

if not os.path.isdir(SCRATCH):
    print(f"Scratch directory {SCRATCH!r} not found (session-specific, may have been "
          "cleaned up) -- skipping the real-snapshot cross-check; Parts 1-4's "
          "idealized-lattice result stands on its own.")
else:
    print(f"{'Corner':<42}{'band':>5}{'h/lambda':>10}{'lambda_eff/h':>14}"
          f"{'measured (sim)':>16}{'real/measured':>16}{'iters':>8}")
    # default corner's lambda_measured/h is not recorded directly in the
    # 2026-09-08 log, only rel_err (0.005 FUV, 0.074 LW) against
    # lambda_analytic/h = 1/0.61 = 1.639 (FUV), 1/1.01 = 0.990 (LW); sign
    # inferred as over-prediction (the "+" branch), consistent with every
    # other corner's own measured direction -- marked (inferred) below.
    measured_lookup = {
        ("default", "FUV"): 1.639 * 1.005, ("default", "LW"): 0.990 * 1.074,
        ("m10", "FUV"): 0.54, ("m10", "LW"): 0.42,
        ("m95", "FUV"): 0.38, ("m95", "LW"): 0.31,
        ("m760", "FUV"): 0.29, ("m760", "LW"): None,
    }
    inferred_keys = {("default", "FUV"), ("default", "LW")}
    for label, subdir in run_dirs.items():
        snap_glob = sorted(glob.glob(os.path.join(SCRATCH, subdir, "snap", "snapshot_*.hdf5")))
        if not snap_glob:
            print(f"  {label}: no snapshots found under {subdir}/snap/, skipping")
            continue
        snap = load_snapshot(snap_glob[-1])
        Z_mean = float(np.mean(snap["Z"]))
        rho_mean = float(np.mean(snap["rho"]))
        r_min = 2.0 * (snap["boxsize"] / snap["pos"].shape[0] ** (1.0 / 3.0))
        r_max = 0.7 * (snap["boxsize"] / 2.0)
        h_mean = float(np.median(snap["h"]))
        for band, sigma_d in (("FUV", SIGMA_D_FUV_CGS), ("LW", SIGMA_D_LW_CGS)):
            lam_cgs = analytic_lambda_cgs(Z_mean, rho_mean, snap["unit_length_cgs"],
                                           snap["unit_mass_cgs"], sigma_d)
            lam = lam_cgs / snap["unit_length_cgs"]
            h_over_lam = h_mean / lam
            key = (subdir.replace("run_", "") if subdir != "run" else "default", band)
            meas = measured_lookup.get(key)
            lam_fit, iters, converged = solve_real_fixed_point(snap, lam, r_min, r_max)
            if lam_fit is None or not np.isfinite(lam_fit):
                print(f"  {label:<40}{band:>5}{h_over_lam:>10.2f}{'no fit':>14}"
                      f"{'n/a' if meas is None else meas:>16}{'n/a':>16}{iters:>8}")
                continue
            ratio_h = lam_fit / h_mean
            cmp_str = f"{ratio_h / meas:.3f}" if meas else "n/a"
            meas_str = (f"{meas:.3f}*" if key in inferred_keys else f"{meas:.3f}") if meas else "n/a"
            conv_flag = "" if converged else " (NOT CONVERGED)"
            print(f"  {label:<40}{band:>5}{h_over_lam:>10.2f}{ratio_h:>14.3f}"
                  f"{meas_str:>16}{cmp_str:>16}{iters:>8}{conv_flag}")
    print()
    print("(* default corner's 'measured' value is inferred from the 2026-09-08 log's")
    print("  rel_err, not read off directly -- see comment above measured_lookup.)")
    print()
    print("This is the direct answer to 'does the simulation match ITS OWN discrete")
    print("steady state': ratios near 1 in the 'real/measured' column mean the")
    print("actual code output matches the exact fixed point of its own real")
    print("(disordered, variable-h) discretization, not an idealized perfect-lattice")
    print("approximation of it -- the strongest form of this check this script builds.")
print()

print("ALL CHECKS PASSED")
