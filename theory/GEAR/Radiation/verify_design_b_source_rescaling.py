"""Verify the Design B source-rescaling fix for the c_hyp-dependent
steady-state amplitude bug (design-lw-fuv-design-b.md Sections 1.1/1.2).

Two independent claims are checked symbolically, exactly (not just to
floating-point tolerance):

1. The ADOPTED fix: rescaling the injection source by
   S_used = S_true * (3*c_hyp/c) makes the c_hyp dependence in the
   steady-state Yukawa amplitude cancel exactly, leaving the amplitude
   implied by the TRUE physical diffusion coefficient D_phys = c*lambda/3.

2. The REJECTED alternative ("split" scheme: reduce only the
   wave-speed/pressure-gradient coefficient in the flux equation to
   c_hyp, keep the flux equation's own damping term at the true rate
   1/tau_phys = c*kappa*rho): this does NOT preserve lambda, giving
   lambda_eff = c_hyp*lambda_true/(sqrt(3)*c), which only equals
   lambda_true when c_hyp = sqrt(3)*c.

Both checks use sympy for exact symbolic algebra; a numeric spot-check
follows each symbolic result as a second, independent confirmation.
"""

import sympy as sp

r, lam, c, c_hyp, S_true = sp.symbols("r lam c c_hyp S_true", positive=True)

# ---------------------------------------------------------------------
# Preliminary: confirm exp(-r/lam)/r solves the homogeneous screened
# (Yukawa/Helmholtz) equation away from the source, i.e. it is the right
# functional form for the Green's function used throughout Sec. 1/6.1.
# ---------------------------------------------------------------------
f = sp.exp(-r / lam) / r
radial_laplacian = sp.simplify((1 / r**2) * sp.diff(r**2 * sp.diff(f, r), r))
residual_homogeneous = sp.simplify(radial_laplacian - f / lam**2)
print("Sanity check: laplacian(f) - f/lam^2 (should be 0 for r>0):")
print("  ", residual_homogeneous)
assert residual_homogeneous == 0, "Yukawa profile does not solve the screened equation"

# ---------------------------------------------------------------------
# Claim 1 (ADOPTED FIX): steady state with D = c_hyp*lam and the
# rescaled source S_used = S_true * 3*c_hyp/c.
# ---------------------------------------------------------------------
D_chosen = c_hyp * lam
S_used = S_true * 3 * c_hyp / c

u_with_fix = S_used / (4 * sp.pi * D_chosen) * f
u_with_fix_simplified = sp.simplify(u_with_fix)

# Physically-expected amplitude: the TRUE diffusion coefficient is
# D_phys = c*lam/3 (standard P1/M1 result, D_phys = c/(3*kappa*rho) with
# lam = 1/(kappa*rho)); the correct amplitude uses D_phys, not D_chosen.
D_phys = c * lam / 3
u_expected = S_true / (4 * sp.pi * D_phys) * f
u_expected_simplified = sp.simplify(u_expected)

residual_fix = sp.simplify(u_with_fix_simplified - u_expected_simplified)
print()
print("Claim 1 (adopted fix): u(r) with S_used, minus the true-physics")
print("prediction 3*S_true/(4*pi*c*lam)*exp(-r/lam)/r (should be exactly 0):")
print("  ", residual_fix)
assert residual_fix == 0, "c_hyp did not cancel exactly"

# Confirm c_hyp has genuinely dropped out of the simplified expression.
print("  u(r) after the fix, fully simplified:", u_with_fix_simplified)
assert c_hyp not in u_with_fix_simplified.free_symbols, "c_hyp still present!"

# Numeric spot-check at two very different c_hyp values: amplitude must
# be IDENTICAL (independent of c_hyp) once the fix is applied.
subs_common = {S_true: 7.0, c: 3e10, lam: 2.5, r: 1.3}
val_slow = float(u_with_fix_simplified.subs({**subs_common, c_hyp: 12.0}))
val_fast = float(u_with_fix_simplified.subs({**subs_common, c_hyp: 4.5e6}))
val_expected = float(u_expected_simplified.subs(subs_common))
print()
print(f"Numeric spot-check: u(c_hyp=12)      = {val_slow:.10e}")
print(f"                    u(c_hyp=4.5e6)  = {val_fast:.10e}")
print(f"                    true prediction = {val_expected:.10e}")
assert abs(val_slow - val_fast) < 1e-12 * abs(val_expected)
assert abs(val_slow - val_expected) < 1e-12 * abs(val_expected)

# ---------------------------------------------------------------------
# Claim 2 (REJECTED alternative): "split" scheme. Flux equation keeps
# the true physical wave-speed-squared/3 coefficient's SPEED replaced by
# c_hyp (i.e. coefficient = c_hyp^2/3, mirroring physical c^2/3 with c ->
# c_hyp only in that one place), but the flux equation's own damping
# term, and the shared u-sink rate, stay at the TRUE physical rate
# 1/tau_phys = c*kappa*rho = c/lam.
#
#   dF/dt + (c_hyp^2/3) grad(u) = -(c/lam) F
#   du/dt + div(F)               = -(c/lam) u + S
#
# Steady state: F = -[(c_hyp^2/3) / (c/lam)] grad(u) = -D_eff * grad(u)
#   with D_eff = c_hyp^2 * lam / (3*c)
# and:            div(F) = -(c/lam)*u + S
#   => laplacian(u) - u/(D_eff * lam/c) = -S/D_eff
#   => lambda_eff^2 = D_eff * (lam/c)
# ---------------------------------------------------------------------
tau_phys = lam / c  # 1/(c*kappa*rho), kappa*rho = 1/lam
D_eff = (c_hyp**2 / 3) * tau_phys
tau_u_true = tau_phys  # "both moments share the same [true] rate"

lam_eff_sq = D_eff * tau_u_true
lam_eff = sp.sqrt(sp.simplify(lam_eff_sq))
lam_eff_simplified = sp.simplify(lam_eff)

target = c_hyp * lam / (sp.sqrt(3) * c)
residual_split = sp.simplify(lam_eff_simplified - target)
print()
print("Claim 2 (rejected alternative): lambda_eff from the split scheme")
print("  lambda_eff =", lam_eff_simplified)
print("  target     =", target)
print("  difference (should be 0):", residual_split)
assert residual_split == 0, "split-scheme lambda_eff does not match the claimed formula"

# Confirm lambda_eff == lambda_true requires c_hyp = sqrt(3)*c exactly.
solutions = sp.solve(sp.Eq(lam_eff_simplified, lam), c_hyp)
print("  c_hyp values making lambda_eff = lam_true:", solutions)
assert solutions == [sp.sqrt(3) * c], "expected c_hyp = sqrt(3)*c as the unique fix point"

print()
print("ALL CHECKS PASSED: source rescaling S_used = S_true*(3*c_hyp/c)")
print("cancels c_hyp exactly; the rejected split alternative provably")
print("does not (lambda_eff = c_hyp*lam/(sqrt(3)*c), fixed only at")
print("c_hyp = sqrt(3)*c, defeating the affordability goal).")
