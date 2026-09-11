"""Verify the RSOL-consistent constants of the M1 closure upgrade, and
both steady-state limits they must reproduce.

This is the executable companion to `theory/GEAR/Radiation/02_fuv_isrf.tex`
Section `fuv-steady-limits`; every claim below is one the theory doc works
through by hand, checked here symbolically and exactly (not just to
floating-point tolerance).

The implemented system (theory doc Eqs. `fuv-rsol-zeroth`/
`fuv-lambda-chyp-solve`) replaces the true light speed by a reduced one
everywhere it appears, and rescales the source by `c_hyp/c`:

    du/dt + div(F)                     = -u/tau + (c_hyp/c)*S
    dF/dt + c_hyp**2*div(D(f)*u)       = -F/tau,     tau = lambda/c_hyp

Checks, in order:

1. Free-streaming branch (f = 1, D = n (x) n): `u = L*exp(-r/lambda)/
   (4*pi*c*r**2)` with `F = c_hyp*u*rhat` satisfies BOTH moment equations,
   is self-consistent (`f = 1` throughout), and its amplitude carries the
   TRUE `1/c`, with `c_hyp` cancelling out.
2. Diffusion branch (f -> 0, D = I/3): `u = 3*L*exp(-sqrt(3)*r/lambda)/
   (4*pi*c*lambda*r)` solves the resulting screened-Poisson equation
   exactly, again with `c_hyp` cancelling; its screening length
   `lambda/sqrt(3)` is the TRUE physical P1 value `sqrt(D_phys*tau_phys)`
   with `D_phys = c*lambda/3`, `tau_phys = lambda/c`.
3. The diffusive branch is NOT self-consistent for a single point source:
   substituting it back into `f` gives `f -> 1/sqrt(3)`, not 0. This is
   why M1 selects the free-streaming branch around one star, and why both
   branches have to be derived rather than one assumed.
4. The closure-independent amplitude identity `Integral(u dV) = lambda*L/c`
   (theory doc Eq. `fuv-steady-amplitude`), which Tier 1 gates on.
5. The SUPERSEDED, P1-tuned constants (pressure coefficient `c_hyp**2`
   with no `1/3`, source rescale `3*c_hyp/c`) reproduced and shown to be
   inconsistent: they give a diffusive screening length `lambda` instead
   of the true physical `lambda/sqrt(3)`, and a total amplitude three
   times the physical one. Getting the first of those to look right was
   the whole point of that tuning; it could not also get the second.
"""

import sympy as sp

r, lam, c, c_hyp, L = sp.symbols("r lam c c_hyp L", positive=True)


def radial_div(vec_r):
    """Divergence of a purely radial vector field with component vec_r(r)."""
    return sp.simplify(sp.diff(r**2 * vec_r, r) / r**2)


print("=" * 70)
print("1. Free-streaming branch: f = 1, D = n (x) n")
print("=" * 70)

u_fs = L * sp.exp(-r / lam) / (4 * sp.pi * c * r**2)
F_fs = c_hyp * u_fs  # radial component; |F| = c_hyp*u, i.e. f = 1 exactly

# Zeroth moment, away from the source: div(F) = -u/tau = -(c_hyp/lam)*u.
res_zeroth_fs = sp.simplify(radial_div(F_fs) + (c_hyp / lam) * u_fs)
print("  zeroth-moment residual (should be 0):", res_zeroth_fs)
assert res_zeroth_fs == 0

# First moment. At chi = 1 the pressure tensor is P = u*n(x)n, purely
# radial, so (div P)_r = du/dr + 2*u/r (the transverse components vanish,
# so the -(P_tt + P_pp)/r term of the spherical divergence drops out).
div_P_fs = sp.diff(u_fs, r) + 2 * u_fs / r
res_first_fs = sp.simplify(c_hyp**2 * div_P_fs + (c_hyp / lam) * F_fs)
print("  first-moment residual (should be 0):", res_first_fs)
assert res_first_fs == 0
print("  -> both moments give the same solution; f = |F|/(c_hyp*u) = 1 by")
print("     construction, so the branch is self-consistent.")

# Amplitude: the near-source limit of the zeroth moment fixes it.
flux_through_small_sphere = sp.limit(4 * sp.pi * r**2 * F_fs, r, 0)
print("  lim_{r->0} 4*pi*r^2*|F| =", sp.simplify(flux_through_small_sphere))
assert sp.simplify(flux_through_small_sphere - c_hyp * L / c) == 0
print("  -> equals (c_hyp/c)*L, exactly the injected, rescaled source.")
assert c_hyp not in u_fs.free_symbols
print("  -> c_hyp does not appear in u(r) at all.")

print()
print("=" * 70)
print("2. Diffusion branch: f -> 0, D = I/3")
print("=" * 70)

u_diff = 3 * L * sp.exp(-sp.sqrt(3) * r / lam) / (4 * sp.pi * c * lam * r)

# F = -(c_hyp*lam/3)*grad(u) from the first moment; substituting into the
# zeroth moment gives laplacian(u) - (3/lam^2)*u = -(3*L/(lam*c))*delta.
F_diff = -(c_hyp * lam / 3) * sp.diff(u_diff, r)
res_diff = sp.simplify(radial_div(F_diff) + (c_hyp / lam) * u_diff)
print("  screened-Poisson residual away from the source (should be 0):", res_diff)
assert res_diff == 0
assert c_hyp not in u_diff.free_symbols
print("  -> c_hyp cancels here too.")

lam_eff = lam / sp.sqrt(3)
lam_phys = sp.sqrt((c * lam / 3) * (lam / c))  # sqrt(D_phys*tau_phys)
print(f"  screening length of this branch: lam/sqrt(3) = {lam_eff}")
print(f"  true physical P1 value sqrt(D_phys*tau_phys) = {sp.simplify(lam_phys)}")
assert sp.simplify(lam_eff - lam_phys) == 0
print("  -> the two agree exactly: lambda/sqrt(3) is the physical answer,")
print("     not an artefact of the reduced light speed.")

print()
print("=" * 70)
print("3. The diffusive branch is not self-consistent for a point source")
print("=" * 70)

f_diff = sp.simplify(sp.Abs(F_diff) / (c_hyp * u_diff))
f_far = sp.limit(f_diff, r, sp.oo)
print("  f(r) on the diffusive branch =", f_diff)
print("  lim_{r->inf} f =", f_far, f"= {float(f_far):.4f}")
assert sp.simplify(f_far - 1 / sp.sqrt(3)) == 0
print("  -> f -> 0.577, not 0: a purely absorbing medium never isotropizes")
print("     a point source's field, so M1 selects the free-streaming branch")
print("     there. The diffusive branch is the many-source background one.")

print()
print("=" * 70)
print("4. Closure-independent amplitude identity")
print("=" * 70)

total_fs = sp.integrate(4 * sp.pi * r**2 * u_fs, (r, 0, sp.oo))
print("  Integral(u dV) on the free-streaming branch =", sp.simplify(total_fs))
assert sp.simplify(total_fs - lam * L / c) == 0
total_diff = sp.integrate(4 * sp.pi * r**2 * u_diff, (r, 0, sp.oo))
print("  Integral(u dV) on the diffusive branch      =", sp.simplify(total_diff))
assert sp.simplify(total_diff - lam * L / c) == 0
print("  -> lambda*L/c on BOTH branches, as it must be: the identity uses")
print("     only that D never enters the zeroth moment.")

print()
print("=" * 70)
print("5. The superseded P1-tuned constants, and why they were replaced")
print("=" * 70)

# Old convention: pressure coefficient c_hyp**2 with no 1/3 (equivalently
# D = c_hyp*lam), source rescaled by 3*c_hyp/c.
D_old = c_hyp * lam
tau = lam / c_hyp
lam_eff_old = sp.sqrt(sp.simplify(D_old * tau))
print("  old diffusive screening length sqrt(D*tau) =", lam_eff_old)
assert sp.simplify(lam_eff_old - lam) == 0
print("  -> exactly lambda: correct-looking, but it is the tuning's own")
print("     target, and it disagrees with the true physical lambda/sqrt(3).")

total_old = sp.simplify(tau * (3 * c_hyp / c) * L)
total_new = sp.simplify(tau * (c_hyp / c) * L)
print("  old total amplitude tau*(3*c_hyp/c)*L =", total_old)
print("  new total amplitude tau*(c_hyp/c)*L   =", total_new)
assert sp.simplify(total_old - 3 * lam * L / c) == 0
assert sp.simplify(total_new - lam * L / c) == 0
print("  -> the old convention put three times the physical energy in the")
print("     box. The two changes (the 1/3 moving inside D, and the source")
print("     rescale) are therefore NOT a cancelling pair: the first sets")
print("     the screening length, the second sets the amplitude.")

print()
print("ALL CHECKS PASSED: the RSOL-consistent constants reproduce both")
print("steady-state limits exactly, with c_hyp cancelling out of both, and")
print("the superseded P1 tuning provably cannot do the same.")
