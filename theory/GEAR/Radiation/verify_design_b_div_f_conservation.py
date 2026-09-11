"""Verify the Design B `div(F)` discretization's conservation property
(design-lw-fuv-design-b.md Section 2.2), for h_i != h_j.

Plan-review finding 1 (2026-09-07 report) showed the design doc's
originally-drafted `div(F)` estimator -- each particle's contribution
finalized *independently* with its own `rho_i * h_i^-(dim+1)` factor --
does NOT cancel exactly under a mass-weighted sum when `h_i != h_j`
(SWIFT's normal case under adaptive smoothing lengths). This script
demonstrates that failure numerically, then demonstrates that the
REVISED construction (mirroring `src/rt/SPHM1RT/rt_gradients.h`'s
`radiation_divergence_SPH`, `diffmode == 1` branch -- a single shared
scalar built from both particles' own kernel-gradient terms inside one
pairwise call, applied with mirrored mass/sign to each side, exactly the
`runner_iact_force`-style pattern the review pointed at) DOES cancel to
machine precision, for the same h_i != h_j configuration.

The quantity checked is exact conservation of `sum_i m_i * u_i` under
pure transport: for a single pairwise interaction, this requires
`m_i * (div F)_i_pair + m_j * (div F)_j_pair == 0` (since `du_i/dt =
-(div F)_i`). This is checked directly, not the divergence's physical
accuracy (which is a separate, well-understood SPH question, §2.3).
"""
# =============================================================================
# M1 CLOSURE AUDIT, 2026-09-11: CHECKED, CLOSURE-INDEPENDENT, NO CHANGE.
# The divergence loop is untouched by the P1-to-M1 upgrade (the closure
# enters only the pressure-tensor term of the flux equation). The exact
# `sum_i m_i u_i` conservation verified here still holds, and is what makes
# the total-amplitude identity Tier 1 now gates on exact.
# =============================================================================
import numpy as np

rng = np.random.default_rng(20260907)


def wendland_c2_3d(q: float) -> tuple[float, float]:
    """Return (W(q), dW/dq) for the 3D Wendland C2 kernel, unit support.

    Matches this project's `--with-kernel=wendland-C2` choice in spirit
    (the conservation identity checked here does not depend on which
    smooth, compactly-supported kernel is used).
    """
    if q >= 1.0:
        return 0.0, 0.0
    norm = 21.0 / (2.0 * np.pi)
    w = norm * (1.0 - q) ** 4 * (4.0 * q + 1.0)
    dwdq = norm * (-4.0 * (1.0 - q) ** 3 * (4.0 * q + 1.0) + 4.0 * (1.0 - q) ** 4)
    return w, dwdq


def kernel_dr(r: float, h: float) -> float:
    """dW/dr for a particle with smoothing length h, SWIFT's `wi_dr`
    convention: `wi_dr = h^-(dim+1) * dW/dq` (dim=3 -> h^-4), i.e. the
    dimension-dependent h-normalization ("h_i^-(dim+1)") is folded in
    HERE, at kernel-evaluation time -- not as a separate later finalize
    step. This is the exact point the flawed draft got wrong (it kept
    the h-normalization as a distinct post-hoc finalize multiplication).
    """
    q = r / h
    _, dwdq = wendland_c2_3d(q)
    return dwdq / h**4


def old_flawed_construction(dx, r, hi, hj, rhoi, rhoj, mi, mj, Fi, Fj):
    """The design doc's ORIGINALLY-DRAFTED (pre-revision) formula:
    each particle's div(F) contribution accumulated and finalized
    independently, using only its own h/kernel/rho.
    """
    r_inv = 1.0 / r
    wi_dr = kernel_dr(r, hi)
    wj_dr = kernel_dr(r, hj)

    # accumulate (raw, using dx as given, i.e. from i's perspective)
    combo = Fi / rhoi**2 + Fj / rhoj**2
    raw_i = mj * np.dot(combo, dx) * wi_dr * r_inv
    # from j's perspective the separation vector is -dx
    raw_j = mi * np.dot(combo, -dx) * wj_dr * r_inv

    # finalize: independently, using each particle's OWN rho_i*h_i^-(dim+1)
    # (h_i^-(dim+1) is already inside wi_dr per kernel_dr above, so the
    # doc's "finalize *= rho_i * h_i^-(dim+1)" reduces to "finalize *=
    # rho_i" given this script's kernel_dr convention -- reproduced
    # faithfully here.)
    div_F_i = raw_i * rhoi
    div_F_j = raw_j * rhoj
    return div_F_i, div_F_j


def new_shared_coefficient_construction(dx, r, hi, hj, rhoi, rhoj, mi, mj, Fi, Fj):
    """The REVISED formula: mirrors `radiation_divergence_SPH`'s
    `diffmode == 1` branch in `src/rt/SPHM1RT/rt_gradients.h`
    (lines 288-297) exactly -- a single shared scalar built from BOTH
    particles' own (F, rho, wi_dr) inside one pairwise call, applied
    with mirrored mass/sign to each side. No separate per-particle
    finalize step.
    """
    r_inv = 1.0 / r
    wi_dr = kernel_dr(r, hi)
    wj_dr = kernel_dr(r, hj)

    shared = (
        np.dot(Fi, dx) / rhoi * wi_dr * r_inv
        + np.dot(Fj, dx) / rhoj * wj_dr * r_inv
    )
    div_F_i = mj * shared
    div_F_j = -mi * shared
    return div_F_i, div_F_j


def check(label, construction, hi, hj, rhoi=1.0, rhoj=1.0):
    dx = np.array([0.37, -0.12, 0.05])
    r = float(np.linalg.norm(dx))
    mi, mj = 1.1, 0.9
    Fi = rng.normal(size=3)
    Fj = rng.normal(size=3)

    div_F_i, div_F_j = construction(dx, r, hi, hj, rhoi, rhoj, mi, mj, Fi, Fj)
    residual = mi * div_F_i + mj * div_F_j
    print(
        f"{label}: hi={hi}, hj={hj}, rhoi={rhoi}, rhoj={rhoj} "
        f"-> mi*divFi + mj*divFj = {residual:.6e}"
    )
    return residual


print("=" * 70)
print("Case A: hi == hj AND rhoi == rhoj (the fully-degenerate case the")
print("flawed draft happens to get right, isolating the effect under test)")
print("=" * 70)
r_old_equal = check("OLD (flawed, independent finalize)", old_flawed_construction, 0.5, 0.5, 1.0, 1.0)
r_new_equal = check("NEW (shared coefficient)          ", new_shared_coefficient_construction, 0.5, 0.5, 1.0, 1.0)
assert abs(r_old_equal) < 1e-12, "OLD should trivially cancel when hi=hj and rhoi=rhoj"
assert abs(r_new_equal) < 1e-12, "NEW must cancel when hi=hj and rhoi=rhoj"

print()
print("=" * 70)
print("Case B: hi != hj, rhoi == rhoj (isolates the h-heterogeneity effect")
print("the plan-review report specifically identified)")
print("=" * 70)
r_old_diff = check("OLD (flawed, independent finalize)", old_flawed_construction, 0.5, 0.9, 1.0, 1.0)
r_new_diff = check("NEW (shared coefficient)          ", new_shared_coefficient_construction, 0.5, 0.9, 1.0, 1.0)

print()
print("=" * 70)
print("Case C: hi != hj AND rhoi != rhoj (the fully general, realistic case)")
print("=" * 70)
r_old_general = check("OLD (flawed, independent finalize)", old_flawed_construction, 0.5, 0.9, 1.3, 0.8)
r_new_general = check("NEW (shared coefficient)          ", new_shared_coefficient_construction, 0.5, 0.9, 1.3, 0.8)
assert abs(r_old_general) > 1e-3, "OLD should also fail to cancel in the fully general case"
assert abs(r_new_general) < 1e-12, "NEW must cancel in the fully general case too"

print()
assert abs(r_old_diff) > 1e-3, (
    "Expected the OLD construction to show a real (non-floating-point-noise) "
    "conservation violation when hi != hj -- if this assertion fires, the "
    "reproduction of the flawed draft above is wrong, not the finding."
)
assert abs(r_new_diff) < 1e-12, (
    "NEW construction must cancel to machine precision for hi != hj -- "
    "this is the whole point of the fix."
)

print("CONFIRMED:")
print(f"  OLD construction: residual = {r_old_diff:.4f} (order-unity, NOT noise) when hi != hj")
print(f"  NEW construction: residual = {r_new_diff:.2e} (machine precision) when hi != hj")
print()
print("Repeating with several random (hi, hj, rho, m, F) draws, hi != hj")
print("enforced every time (rho also allowed to vary, the fully general")
print("case), to confirm this is not a fluke of one configuration:")
max_old_resid = 0.0
max_new_resid = 0.0
for trial in range(200):
    dx = rng.normal(size=3) * 0.4 + np.array([0.3, 0.0, 0.0])
    r = float(np.linalg.norm(dx))
    hi = rng.uniform(0.3, 1.2)
    hj = hi + rng.uniform(0.05, 0.6)  # force hi != hj on every draw
    rhoi, rhoj = rng.uniform(0.5, 2.0, size=2)
    mi, mj = rng.uniform(0.5, 2.0, size=2)
    Fi = rng.normal(size=3)
    Fj = rng.normal(size=3)

    di, dj = old_flawed_construction(dx, r, hi, hj, rhoi, rhoj, mi, mj, Fi, Fj)
    old_resid = abs(mi * di + mj * dj)
    max_old_resid = max(max_old_resid, old_resid)

    di2, dj2 = new_shared_coefficient_construction(dx, r, hi, hj, rhoi, rhoj, mi, mj, Fi, Fj)
    new_resid = abs(mi * di2 + mj * dj2)
    max_new_resid = max(max_new_resid, new_resid)

print(f"  max |residual| over 200 random hi!=hj draws, OLD: {max_old_resid:.4e}")
print(f"  max |residual| over 200 random hi!=hj draws, NEW: {max_new_resid:.4e}")
assert max_new_resid < 1e-10
assert max_old_resid > 1e-2
print()
print("PASS: NEW (shared-coefficient) construction telescopes to exact")
print("(machine-precision) conservation for h_i != h_j; OLD (independently")
print("finalized) construction does not.")
