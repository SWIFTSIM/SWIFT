"""Adjudicate whether Tier 1's pure-Yukawa fit target is the theoretically
correct thing to fit against Design B's own steady state, in the regime
Tier 1 actually exercises (finite absorption, i.e. `kappa > 0`,
`lambda < inf`) -- not the free-streaming (`kappa -> 0`) limit the design
doc's "Known gaps" `~lambda/(3r)` note is about.

Background (design-lw-fuv-design-b.md Sec. 1.2, 6.1, "Known gaps"):
Sec. 1.2 derives that Design B's coupled steady state, for a point source
in a uniform medium, reduces to the screened-Poisson ("modified Helmholtz")
equation `laplacian(u) - u/lambda^2 = -S/D`, and asserts its point-source
Green's function is the Yukawa profile `u(r) = S/(4*pi*D)*exp(-r/lambda)/r`.
`isrf_yukawa_profile_check.py` (Tier 1) fits the SIMULATED radial profile
against exactly that functional form, over `r/lambda` roughly in [1, 9].

Separately, the "Known gaps" section notes the design's Yukawa steady state
"vanishes as lambda -> inf at fixed r, while the true transparent-medium
field is `S_true/(4*pi*r^2*c)`: ratio `3*r/lambda*exp(-r/lambda)`, i.e. the
design under-predicts the steady ISRF in gas with `lambda >> r` by
`~lambda/(3*r)`" -- a real, P1-closure-vs-true-radiative-transfer gap.

The question this script settles: is that `~lambda/(3*r)` note a defect
IN the P1 system's own exact steady-state solution (which would mean
fitting the simulation's output with a pure-Yukawa functional form is
itself biased), or is it a comparison of the P1 system's exact solution
against a DIFFERENT, non-P1 physical target (which would mean the note is
irrelevant to whether Tier 1's fit methodology is sound)?

Four things are established below, symbolically/numerically, not by
assertion:
1. The general steady-state solution of the homogeneous radial ODE
   (`laplacian(u) - u/lambda^2 = 0`, r>0) has an exp(-r/lambda)/r branch
   and an exp(+r/lambda)/r branch; boundedness at r->inf kills the growing
   branch, leaving PURE Yukawa exactly -- for ANY finite lambda>0, at
   every r>0, not merely as an asymptotic approximation.
2. Direct substitution confirms zero ODE residual for the pure-Yukawa
   ansatz at several lambda values spanning many orders of magnitude
   (including values that bracket Tier 1's own measured `h/lambda` corner)
   -- i.e. no r/lambda-dependent shape correction exists in the P1
   system's own exact solution at any finite lambda.
3. The `~lambda/(3*r)` ratio is reproduced exactly as the ratio of the
   P1 system's Yukawa solution to a SEPARATE closed-form expression (the
   kappa=0 transparent-medium field), which does not solve the same ODE
   at all (kappa=0 gives a plain Poisson equation, not a screened one) --
   confirming it is a cross-model comparison, not an internal defect.
4. The ratio is tabulated over Tier 1's actual fit range (r/lambda in
   [1, 9]) and shown to sit on the OPPOSITE side of the "lambda >> r" limit
   the note's own qualifier requires for the "under-predicts" reading to
   even apply -- Tier 1's own regime has r >= lambda throughout, so the
   asymptotic "~lambda/(3*r) under-prediction" language does not describe
   what is happening there in the first place (at r/lambda=1 the design's
   OWN steady state is actually 10% ABOVE the transparent-medium value,
   not below it; only past r/lambda ~ 3 does it fall below, and by
   r/lambda=9 the transparent-medium comparison is moot because the
   medium is manifestly NOT transparent over that distance, kappa*r ~ 9).

A fifth, SEPARATE point (Part 6 below) is flagged, not settled: the above
four points are about the CONTINUUM fit target's functional form only.
Whether the DISCRETE SPH estimator's fitted lambda is itself unbiased at
Tier 1's own corner (lambda comparable to h, per the run this session
actually validated) is a resolution/discretization question the design
doc's own Sec. 6.1 already names separately and this script does not
close -- see Part 6's caveat and named follow-up.

Conclusion printed at the end: Tier 1's pure-Yukawa fit FUNCTIONAL FORM is
the theoretically correct target for validating Design B's own P1
numerics, at any finite kappa; the `~lambda/(3*r)` note is a separate,
already-flagged, deliberately-deferred P1-vs-true-RT closure question (the
cheap-M1 upgrade path), not a shape-fit bias affecting that functional
form. A distinct, NOT-YET-RULED-OUT discretization-bias question (fitted
lambda vs. true lambda at Tier 1's own lambda~h corner) is flagged in Part
6 and left open, per the assigned diagnosis-only scope.
"""

import numpy as np
import sympy as sp

# =====================================================================
# Part 1: exact general solution of the homogeneous steady-state radial
# ODE, spherically symmetric, at FINITE lambda (kappa > 0). This is the
# regime Tier 1 actually runs at -- not the kappa -> 0 limit.
# =====================================================================
print("=" * 78)
print("Part 1: general steady-state solution, finite lambda (kappa > 0)")
print("=" * 78)

r, lam, D, S, c, c_hyp, S_true, w = sp.symbols(
    "r lam D S c c_hyp S_true w", positive=True
)
u = sp.Function("u", positive=True)

# Standard substitution w(r) = r*u(r) turns the radial Laplacian into a
# plain constant-coefficient 2nd-order ODE -- sympy solves this reliably;
# the substitution itself is exact (not an approximation).
w_func = sp.Function("w")
ode_w = sp.Eq(sp.diff(w_func(r), r, 2) - w_func(r) / lam**2, 0)
general_w = sp.dsolve(ode_w, w_func(r))
print("ODE for w=r*u:", ode_w)
print("General solution:", general_w)

C1, C2 = sp.symbols("C1 C2")
# sympy's dsolve returns some linear combination of exp(+-r/lam); extract
# the RHS and confirm it is exactly A*exp(-r/lam) + B*exp(r/lam) for some
# constants (verify by matching, not by assuming the printed form).
rhs = general_w.rhs
expanded = sp.expand(rhs.rewrite(sp.exp))
print("Rewritten in terms of exp():", expanded)

decaying = sp.exp(-r / lam)
growing = sp.exp(r / lam)
assert expanded.has(decaying) or expanded.has(growing), (
    "general solution not expressed in terms of exp(+-r/lam) as expected"
)

print()
print("Physical selection rule: u(r) = w(r)/r must stay BOUNDED as r->inf")
print("(a finite steady-state energy density far from an isolated source,")
print("not a runaway). The exp(+r/lam) branch diverges as r->inf for any")
print("finite lam>0, so its coefficient must vanish. This is a boundary")
print("condition on this exact ODE, not an approximation valid only for")
print("some range of r/lam -- it holds at every r>0.")
print()

# Confirm explicitly: the growing branch violates boundedness for ANY
# lambda in (0, inf), i.e. this selection is not lambda-dependent.
lam_samples = [sp.Rational(1, 1000), sp.Integer(1), sp.Integer(1000), sp.Integer(10**6)]
for lam_val in lam_samples:
    limit_growing = sp.limit((growing / r).subs(lam, lam_val), r, sp.oo)
    assert limit_growing == sp.oo, f"growing branch unexpectedly bounded at lam={lam_val}"
print("Verified: exp(+r/lam)/r diverges as r->inf at every sampled lambda")
print(f"  ({[str(v) for v in lam_samples]}) -- boundedness always excludes it.")

# So the bounded steady-state solution, for ANY finite lambda>0, at EVERY
# r>0, is exactly the single-exponential Yukawa form.
u_general = sp.Symbol("A") * decaying / r
print()
print("=> Bounded steady state, any finite lambda, any r>0:")
print("   u(r) = A * exp(-r/lam) / r      (PURE Yukawa, exact, not asymptotic)")

# =====================================================================
# Part 2: direct residual check -- plug the pure-Yukawa ansatz into the
# FULL inhomogeneous screened-Poisson ODE and confirm zero residual
# (away from r=0) at several lambda values spanning many decades,
# including values that bracket Tier 1's own measured h/lambda corner
# (design doc Sec. 6.1: h/lambda_FUV ~ 0.6, h/lambda_LW ~ 1.0, i.e.
# lambda and h are comparable there -- well inside the "finite lambda"
# regime, nowhere near lambda -> inf).
# =====================================================================
print()
print("=" * 78)
print("Part 2: zero-residual check of the pure-Yukawa ansatz, several lambda")
print("=" * 78)

A = sp.Symbol("A", positive=True)
u_ansatz = A * sp.exp(-r / lam) / r
radial_laplacian = sp.simplify((1 / r**2) * sp.diff(r**2 * sp.diff(u_ansatz, r), r))
residual = sp.simplify(radial_laplacian - u_ansatz / lam**2)
print("laplacian(u_ansatz) - u_ansatz/lam^2 (symbolic, should be 0 for r>0):")
print("  ", residual)
assert residual == 0, "Yukawa ansatz does not solve the screened ODE exactly"

# Numeric spot-check across decades of lambda, evaluated at several r
# values within Tier 1's own fit range for each lambda (r/lam in [1, 9]).
print()
print("Numeric spot-check, residual at r = k*lam for k in [1..9], several lambda:")
max_abs_residual = 0.0
for lam_val in [1e-3, 1.0, 2.6, 100.0, 1e6]:
    for k in range(1, 10):
        r_val = k * lam_val
        subs = {A: 3.7, lam: lam_val, r: r_val}
        res_val = float(residual.subs(subs)) if residual != 0 else 0.0
        max_abs_residual = max(max_abs_residual, abs(res_val))
print(f"  max |residual| over the sweep: {max_abs_residual:.3e} (exactly 0 symbolically)")
assert max_abs_residual == 0.0
print("Confirmed: NO r/lambda-dependent shape correction exists in the P1")
print("system's own exact steady-state solution, at ANY finite lambda.")
print("The fit target the simulation's own numerics should reproduce is")
print("pure Yukawa at every corner Tier 1 could plausibly probe.")

# =====================================================================
# Part 3: the ~lambda/(3r) note is a CROSS-MODEL comparison, not an
# internal defect. Reproduce it exactly, and show it does not solve the
# same ODE the P1 Yukawa profile solves.
# =====================================================================
print()
print("=" * 78)
print("Part 3: the ~lambda/(3r) note is a comparison to a DIFFERENT model")
print("=" * 78)

D_phys = c * lam / 3
u_p1_pointsource = S_true / (4 * sp.pi * D_phys) * sp.exp(-r / lam) / r
u_p1_pointsource = sp.simplify(u_p1_pointsource)
print("P1 exact point-source steady state (design's own Sec. 1.2 result):")
print("  u_P1(r) =", u_p1_pointsource)

u_transparent = S_true / (4 * sp.pi * r**2 * c)
print("True kappa=0 (transparent-medium, free-streaming) field:")
print("  u_transparent(r) =", u_transparent)

# Confirm u_transparent does NOT solve the screened ODE (it solves the
# plain 1/r^2 flux-conservation relation instead) -- i.e. it is not a
# member of the same solution family, it is a genuinely different model.
residual_transparent = sp.simplify(
    (1 / r**2) * sp.diff(r**2 * sp.diff(u_transparent, r), r) - u_transparent / lam**2
)
print()
print("Does u_transparent solve the SAME screened ODE laplacian(u)-u/lam^2=0?")
print("  residual =", residual_transparent, " (nonzero => different model, as expected)")
assert residual_transparent != 0, "transparent-medium field should NOT solve the screened ODE"

ratio = sp.simplify(u_p1_pointsource / u_transparent)
print()
print("Ratio u_P1(r) / u_transparent(r):")
print("  ", ratio)
x = sp.Symbol("x", positive=True)  # x := r/lambda
ratio_of_x = sp.simplify(ratio.subs(r, x * lam))
print("  as a function of x=r/lambda:", ratio_of_x)
assert sp.simplify(ratio_of_x - 3 * x * sp.exp(-x)) == 0, (
    "ratio does not match the design doc's own quoted 3*r/lambda*exp(-r/lambda)"
)
print("Confirmed: matches the design doc's own quoted ratio, 3*r/lambda*exp(-r/lambda).")
print()
print("This ratio compares the P1 system's exact solution (u_P1, which Tier 1")
print("correctly targets) against a DIFFERENT closed-form expression that does")
print("NOT satisfy the same governing equation -- it is not a correction term")
print("that a more careful solve of the P1 system would reveal sitting inside")
print("u_P1 itself. Part 2 already showed u_P1 has zero residual identically.")

# Small-x (r << lambda) asymptotic reading, exactly as the design doc's
# qualifier states it ("lambda >> r"):
small_x_limit = sp.limit(ratio_of_x / x, x, 0)
print()
print(f"Small-x (r<<lambda) limit of ratio/x: {small_x_limit}  "
      "(=> ratio ~ 3x, i.e. u_transparent/u_P1 ~ 1/(3x) = lambda/(3r),")
print("   which is exactly the design doc's own asymptotic reading -- valid")
print("   ONLY as a leading-order approximation for x << 1.")

# =====================================================================
# Part 4: evaluate the ratio across Tier 1's OWN fit range, r/lambda in
# [1, 9], and show it is on the OPPOSITE side of the "lambda >> r"
# (small-x) regime the note's own qualifier requires.
# =====================================================================
print()
print("=" * 78)
print("Part 4: the ratio, evaluated across Tier 1's own fit range r/lambda in [1,9]")
print("=" * 78)

ratio_numeric = sp.lambdify(x, ratio_of_x, "numpy")
xs = np.arange(1, 10)
ratios = ratio_numeric(xs)

print(f"{'r/lambda':>10} | {'u_P1/u_transparent':>20} | {'interpretation':>40}")
print("-" * 78)
for xv, rv in zip(xs, ratios):
    if rv > 1.0:
        note = f"P1 is {100*(rv-1):.1f}% ABOVE transparent-field value"
    else:
        note = f"P1 is {100*(1-rv):.2f}% below (medium optically thick here)"
    print(f"{xv:10d} | {rv:20.6e} | {note:>40}")

assert ratios[0] > 1.0, "expected P1 to exceed the transparent-field value at r/lambda=1"
print()
print("At the LOW end of Tier 1's fit range (r/lambda=1), the design's OWN")
print("exact P1 steady state is *larger* than the transparent-medium value,")
print("not smaller by ~33% as a naive small-x extrapolation of the note would")
print("suggest -- because r/lambda=1 is order-unity optical depth, not the")
print("optically-thin (lambda>>r) regime the note's qualifier requires. At")
print("the HIGH end (r/lambda=9), the medium has kappa*r~9 (manifestly NOT")
print("transparent over that distance) so comparing against the transparent-")
print("medium field there is not even a physically meaningful comparison --")
print("of course a screened field is tiny relative to an unscreened one once")
print("several optical depths separate source and observer; that is exactly")
print("what absorption is supposed to do, not evidence of a P1 defect.")

# =====================================================================
# Part 5: connect explicitly to the free-streaming (kappa->0) pathology.
# Confirm it is a SINGULAR point (lambda -> inf at fixed r), not a smooth
# correction encroaching on finite-lambda solutions from that limit.
# =====================================================================
print()
print("=" * 78)
print("Part 5: free-streaming (kappa->0) pathology -- singular limit, not a")
print("        finite-lambda remnant")
print("=" * 78)

# u_P1 at fixed r as lambda -> inf:
limit_u_p1 = sp.limit(u_p1_pointsource, lam, sp.oo)
print("lim_{lambda->inf} u_P1(r) at FIXED r:", limit_u_p1)
assert limit_u_p1 == 0, "expected the P1 point-source field to vanish as lambda->inf"

print("The transparent-medium field u_transparent has NO lambda dependence")
print("at all (kappa=0 by construction) and is finite and nonzero at every r.")
print()
print("This is the free-streaming pathology: as kappa->0, the diffusion")
print("coefficient D=c_hyp*lambda -> inf, which dilutes the P1 steady-state")
print("energy density to zero at any fixed r -- a SINGULAR limit at exactly")
print("lambda=inf, not a shape defect that develops continuously as lambda")
print("grows large while still finite. Part 2 already confirmed the ODE")
print("residual is EXACTLY zero at lambda=1e6 (a value 1e5-1e6 times Tier")
print("1's own measured lambda scale) -- i.e. right up to the singular point")
print("itself, the P1 system's OWN exact solution is still pure Yukawa, with")
print("no gradual creeping-in of an r/lambda-shaped correction term.")
print()
print("Tier 1's own regime (design doc Sec. 6.1, this session's actual run):")
print("  h/lambda_FUV ~ 0.6, h/lambda_LW ~ 1.0, fit range r/lambda in [1, 9.2]")
print("is nowhere near this singular limit; lambda there is a finite length")
print("comparable to the resolution h, not large compared to the box.")

# =====================================================================
# Part 6: caveat, explicitly out of scope for this script -- a DIFFERENT
# bias mechanism the above analysis does NOT rule out. The design doc's
# own Sec. 6.1 already names it: "for lambda < h the discrete profile's
# e-folding length is floored at ~0.25 h", i.e. once lambda and the SPH
# smoothing length h are comparable, the DISCRETE estimator's own steady
# state deviates from the CONTINUUM Yukawa profile this script analyzed
# -- a resolution effect, unrelated to the P1-vs-transparent-medium
# question Parts 1-5 settle. Tier 1's own measured run sits at
# h/lambda_FUV ~ 0.6 and h/lambda_LW ~ 1.0 (2026-09-07 validation log),
# i.e. lambda ~ h to 1.7*h -- NOT the lambda >> h regime where the
# continuum slope this script derives is guaranteed to be what the
# discrete estimator actually converges to.
# =====================================================================
print()
print("=" * 78)
print("Part 6: caveat -- discretization bias at lambda~h is a SEPARATE,")
print("        NOT-YET-RULED-OUT question")
print("=" * 78)
print(
    """
Parts 1-5 prove the CONTINUUM target is exactly Yukawa at every finite
lambda, with no r/lambda-dependent shape correction. That answers the
task's (a)-vs-(b) question about the fit's FUNCTIONAL FORM.

It does NOT prove the discrete SPH estimator's own steady state converges
to that continuum Yukawa profile's lambda value once lambda and the
smoothing length h are comparable -- the design doc's own Sec. 6.1
explicitly flags a resolution-dependent floor in this regime ("for
lambda < h the discrete profile's e-folding length is floored at
~0.25 h"), and Tier 1's own already-run validation sits at
h/lambda_FUV ~ 0.6, h/lambda_LW ~ 1.0 -- lambda comparable to h, not
lambda >> h. Consistent with (not proof of) a real discretization bias at
this corner: the 2026-09-07 validation run recovered lambda_LW (the more
under-resolved band, h/lambda=1.0) to only 7.7%, about 6x worse than
lambda_FUV (h/lambda=0.6, recovered to 1.2%) -- the direction the
floor argument predicts, though not isolated as the confirmed mechanism
here. Two runs agreeing with each other at the SAME h/lambda corner (as
the matched-nu re-run did, to <1%) would still share this bias if it is
present; that agreement is evidence of reproducibility, not evidence the
recovered lambda equals the true continuum lambda.

This is a distinct, NOT YET RULED OUT question from the one this script
answers, and this script does not attempt to close it. The design doc's
own linear discrete solve (u - lambda^2*div_h(grad_h(u)) = tau*S, the
same machinery verify_design_b_timestepping_stability.py Part C already
implements) is the concrete next check: fit the same pure-Yukawa
estimator against THAT discrete prediction instead of a real run, at
Tier 1's own particle distribution and h/lambda, and read off
lambda_fit/lambda_true as the discretization bias at this corner. Not
built here -- this is a diagnosis of the functional-form question only,
per the assigned scope; the discretization-bias check is a separate,
named follow-up.
"""
)

# =====================================================================
# Conclusion
# =====================================================================
print()
print("=" * 78)
print("CONCLUSION")
print("=" * 78)
print(
    """
(a) is the correct answer to the design doc's own open FUNCTIONAL-FORM
question. Design B's P1-closed system has an EXACT steady-state solution,
at every finite lambda (any nonzero kappa) and every r>0, that is a PURE
Yukawa profile -- not a Yukawa profile plus a lambda/(3r)-type (or any
other) correction term. This was verified two independent ways above: (1)
solving the homogeneous radial ODE directly and showing boundedness kills
the growing exponential branch at every finite lambda (Part 1), and (2)
direct symbolic substitution giving an identically zero ODE residual for
the Yukawa ansatz at lambda values spanning 1e-3 to 1e6, evaluated at
r/lambda = 1..9 in each case (Part 2).

Therefore fitting Tier 1's simulated radial profile against a pure-Yukawa
FUNCTIONAL FORM is the theoretically correct thing to do, not an
approximation of convenience: it is fitting the code's numerical output
against the EXACT analytic solution of the equations Design B actually
solves. No change to isrf_yukawa_profile_check.py's fit functional form
is warranted on these grounds.

The "~lambda/(3r)" note in the design doc's "Known gaps" section is
unrelated to this fit-methodology question. It compares the P1 system's
exact solution against a DIFFERENT closed-form expression -- the true
kappa=0 transparent-medium field -- which does not satisfy the same
governing equation at all (Part 3) and is only the physically relevant
comparison in the optically-thin limit r << lambda, the OPPOSITE end of
Tier 1's own fit range (r/lambda in [1,9], i.e. r >= lambda throughout;
Part 4). Concretely: at r/lambda=1 the design's OWN exact steady state is
10.4% ABOVE the transparent-medium value, not ~33% below it as a naive
lambda/(3r) extrapolation (dropping the exp(-x) factor) would suggest. It
is a P1-closure-vs-true-radiative-transfer fidelity question (already
named in the design doc, already deferred behind the cheap-M1 upgrade
path per the operator's "start with P1" ruling), not a shape-fit bias in
how Tier 1 measures whether the CODE correctly solves the P1 equations.

SEPARATE, NOT-RULED-OUT CAVEAT (Part 6): whether the discrete SPH
estimator's fitted lambda is itself unbiased at Tier 1's own corner
(lambda ~ h to 1.7h) is a resolution/discretization question this script
does not address, already flagged in the design doc's own Sec. 6.1 "floor"
note, and consistent with (not confirmed by) the measured 6x-worse
lambda_LW recovery (h/lambda=1.0) vs lambda_FUV (h/lambda=0.6) in the
2026-09-07 run. Two independent runs agreeing with each other at that same
corner does not rule this out -- it would just mean both share the same
discretization bias, not that the bias is absent. Recommended follow-up
(not built here, diagnosis-only per task scope): fit the same estimator
against the design doc's own discrete linear solve at Tier 1's actual
h/lambda and particle distribution, and read off lambda_fit/lambda_true
as the discretization bias directly.

RECOMMENDATION: no change needed to isrf_yukawa_profile_check.py's
pure-Yukawa FUNCTIONAL FORM -- it targets the exact continuum solution of
the equations being tested, at any finite kappa, with no shape correction
missing. Whether the CURRENT h/lambda corner Tier 1 happens to run at
extracts an unbiased lambda from that (correct) functional form is a
separate, open discretization question, not settled by this script, and
should be checked via the discrete linear-solve comparison named above
before treating the ~1% run-to-run reproducibility as proof of accuracy
rather than proof of reproducibility.
"""
)
