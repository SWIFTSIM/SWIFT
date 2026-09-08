"""Verify the Lagrangian/mass-specific re-derivation of Design B's zeroth-
moment governing equation (design-lw-fuv-design-b.md Section 1, this
session's addendum).

**Why this script exists.** Two prior `/plan-review` rounds established,
empirically (`verify_design_b_div_f_consistency.py`), that the §2.2
`div(F)` SPH estimator computes `(1/rho)*div(rho*F)`, not the plain
`div(F)` that Section 1's governing equation -- copied unmodified from the
theory doc's own **volumetric** RTE-moment convention
(`02_fuv_isrf.tex` Eq. fuv-p1-zeroth/fuv-p1-first) -- literally calls for.
The operator ruled: `u`/`F` are mass-specific throughout this
implementation (matching injection's own `u_i += u_inject/m_i`, and the
shipped `FUVSpecificEnergy`/`LWSpecificEnergy` I/O names), and Section 1
must be re-derived into its correct Lagrangian/mass-specific form rather
than have the discretization chase the doc.

**What this script actually checks, precisely** -- not a tautological
rearrangement, but the one non-obvious, error-prone step: converting an
EULERIAN volumetric PDE (Eq. fuv-p1-zeroth, `d(u_V)/dt` meaning
`partial_t` at a FIXED point in space) into a LAGRANGIAN equation for the
mass-specific quantity `u = u_V/rho`, tracked by an SPH particle that
MOVES with the local fluid velocity `v`. This is genuinely different from
the "standard SPH trick" that converts a volumetric internal-energy
equation into the specific internal-energy equation: that trick works
because the Eulerian energy-density equation is written with an EXPLICIT
mass-advection term, `d(rho*e)/dt + div(rho*e*v) = ...` (the flux already
contains a `rho*e*v` piece by construction, since internal energy is
carried by the gas itself). Eq. fuv-p1-zeroth has NO such term: `F_V` is a
genuinely independent radiative flux, not `u_V*v` plus something else --
radiation is not advected with the gas at leading order. Converting
`partial_t` to `D/Dt` for a field that is NOT advected by construction
therefore genuinely introduces an extra term (a pure kinematic
Reynolds-transport-theorem effect, present regardless of what physics
governs `u_V`'s own dynamics) that the "just divide through by rho"
shortcut silently drops. This script derives that term explicitly with
sympy, rather than asserting it survives or cancels by analogy.

**Result, stated up front**: the dilation piece of that extra term
cancels EXACTLY against the mass-continuity term that appears when
differentiating `u = u_V/rho` itself -- this is the genuinely non-obvious,
previously un-verified step (flagged by round-1 `/plan-review`, sec 4.6,
as "asserted by the standard SPH argument but not rigorously re-derived
term by term"). What is LEFT OVER after that cancellation is a single
term, `(1/rho)*div(u_V*v)`, that the discretization does not compute (no
`v_i`/`v_j` enters `grad(u)`/`div(F)`/the local finalize).

**How to read that term (round 4 of the design doc, superseding this
script's original "missing physics" framing; the algebra below is
unchanged and still correct).** The term is the difference between
writing the REDUCED-speed-of-light P1 system in the simulation-box frame
(Eq. fuv-p1-zeroth as a box-frame field equation, then converted to
Lagrangian coordinates: the term appears) and in the local gas frame (the
RTE moments hold there with comoving quantities; the term does not
appear). Neither is exact: the true P1 system carries O(v/c) mixed-frame
corrections this design never models, and replacing c by c_hyp inflates
whatever is neglected to O(v/c_hyp). In the box frame that is
O(v_box/c_hyp) ~ 30-100 in a galaxy (a star moving faster than c_hyp
relative to the box would trail a Mach cone of its own field, and the
answer would depend on the box's rest frame); in the gas frame it is
O(|v_star - v_gas|/c_hyp), a few km/s over c_hyp, since stars form from
and move with the gas they illuminate. The discretization's gas-frame
choice (drop the term) is therefore the better one, and Step 8's
`v_gas/c_hyp` ratio below is the size of the box-frame price, not of an
error in the gas-frame scheme. What the gas frame costs instead is a
steady centroid lag of exactly `v_rel*tau` for a source moving through
gas (`verify_design_b_timestepping_stability.py`, Part D). GEAR-RT
(`src/rt/GEAR/rt_gradients.h::rt_gradients_predict_drift`) applies
precisely this term at every drift because it chose the box frame, which
is right there since its c_red >> v_box.
"""

import sympy as sp

# ---------------------------------------------------------------------
# Independent variables and fields, all genuinely position/time-dependent
# (sympy Function objects, not scalars) -- required to get partial
# derivatives, chain rule, and product rule right without hand algebra.
# ---------------------------------------------------------------------
t, x, y, z = sp.symbols("t x y z")
coords = (x, y, z)
c, kappa, tau, S = sp.symbols("c kappa tau S", positive=True)
# S is treated as the (already mass-specific) source rate S_V/rho; kept as
# a plain symbol (not a Function) since it never gets differentiated below
# -- only its bookkeeping (S_V = rho*S) matters for this derivation.

rho = sp.Function("rho")(x, y, z, t)
u = sp.Function("u")(x, y, z, t)  # the MASS-SPECIFIC field, u := u_V/rho
vx, vy, vz = (sp.Function(f"v_{c_}")(x, y, z, t) for c_ in "xyz")
v = (vx, vy, vz)
FVx, FVy, FVz = (sp.Function(f"F_V{c_}")(x, y, z, t) for c_ in "xyz")
F_V = (FVx, FVy, FVz)


def div(vec):
    return sum(sp.diff(component, coord) for component, coord in zip(vec, coords))


def grad(scalar):
    return tuple(sp.diff(scalar, coord) for coord in coords)


def dot(a, b):
    return sum(ai * bi for ai, bi in zip(a, b))


def material_derivative(field):
    """D/Dt = partial_t + v . grad, the operator an SPH particle's own
    d(field_i)/dt physically represents (the particle moves at v)."""
    return sp.diff(field, t) + dot(v, grad(field))


u_V = rho * u  # volumetric := mass-specific * density, by definition

print("=" * 78)
print("Step 0: setup")
print("=" * 78)
print("  u_V := rho * u   (u is the mass-specific field being solved for)")
print("  Eulerian governing PDE (Eq. fuv-p1-zeroth, volumetric):")
print("    partial_t(u_V) + div(F_V) = -u_V/tau + rho*S")
print()

# ---------------------------------------------------------------------
# Step 1: the Eulerian PDE, solved for partial_t(u_V).
# ---------------------------------------------------------------------
dudt_V_rhs = -div(F_V) - u_V / tau + rho * S

# ---------------------------------------------------------------------
# Step 2: mass continuity (Eulerian form), solved for partial_t(rho).
# This is the ONLY place fluid-dynamics input (continuity) enters --
# nothing about radiation transport is assumed here.
# ---------------------------------------------------------------------
rho_v = tuple(rho * vi for vi in v)
drhodt_rhs = -div(rho_v)

# ---------------------------------------------------------------------
# Step 3: substitute BOTH into the definition u_V = rho*u, differentiated
# in time, and solve for partial_t(u) -- the only unknown left once the
# two PDEs above have replaced partial_t(u_V) and partial_t(rho).
# ---------------------------------------------------------------------
ut_symbol = sp.diff(u, t)
rhot_symbol = sp.diff(rho, t)

# d(u_V)/dt, expanded via the product rule in terms of ut_symbol/rhot_symbol
dudt_V_expanded = sp.diff(u_V, t)  # = rho*u_t + u*rho_t (sympy expands this)

# Solve "dudt_V_expanded == dudt_V_rhs" for u_t, after eliminating rho_t
# via continuity (substitute FIRST, since u_t appears linearly only once
# rho_t has been replaced).
eq_for_ut = sp.Eq(dudt_V_expanded.subs(rhot_symbol, drhodt_rhs), dudt_V_rhs)
ut_solutions = sp.solve(eq_for_ut, ut_symbol)
assert len(ut_solutions) == 1, "expected a unique solution for partial_t(u)"
ut_derived = sp.simplify(ut_solutions[0])

# ---------------------------------------------------------------------
# Step 4: assemble the material derivative Du/Dt = partial_t(u) + v.grad(u)
# using the DERIVED partial_t(u) -- this is what an SPH particle's own
# d(u_i)/dt update must equal to be physically correct.
# ---------------------------------------------------------------------
Du_Dt_derived = sp.simplify(ut_derived + dot(v, grad(u)))

print("=" * 78)
print("Step 1-4: Du/Dt derived from the Eulerian PDE + mass continuity,")
print("no further assumptions")
print("=" * 78)
print("  Du/Dt =")
sp.pprint(Du_Dt_derived)
print()

# ---------------------------------------------------------------------
# Step 5: THE CLAIM -- this equals the "relative flux" form,
#   Du/Dt = -(1/rho)*div(F_V - u_V*v) - u/tau + S
# i.e. the flux appearing in the correctly Lagrangian-converted equation
# is F_V MINUS the purely-kinematic advective piece u_V*v (a standard
# Reynolds-transport-theorem correction, present for ANY scalar field
# carried across a moving control-volume boundary, independent of the
# field's own physics) -- not F_V unmodified.
# ---------------------------------------------------------------------
G = tuple(F_V[i] - u_V * v[i] for i in range(3))  # "relative-to-fluid" flux
claim_relative_flux = -div(G) / rho - u / tau + S

residual_claim1 = sp.simplify(Du_Dt_derived - claim_relative_flux)
print("=" * 78)
print("Claim 1: Du/Dt == -(1/rho)*div(F_V - u_V*v) - u/tau + S")
print("(residual, should be exactly 0):")
print("=" * 78)
print("  ", residual_claim1)
assert residual_claim1 == 0, "Relative-flux identity failed -- re-derive by hand"
print("PASS")
print()

# ---------------------------------------------------------------------
# Step 6: rewrite in terms of the DISCRETIZATION'S actual tracked
# per-mass flux, F := F_V/rho (the plain, un-reinterpreted substitution --
# deliberately NOT the "comoving-frame flux" reinterpretation, which was
# considered and rejected: the standard mixed-frame RT result is
# F_lab = F_comoving + (4/3)*u_comoving*v at the Eddington closure, not
# u_comoving*v with coefficient 1, so identifying this script's purely
# KINEMATIC coefficient-1 term with the physically-distinct O(v/c)
# comoving-frame correction would be an unverifiable overclaim this
# script cannot check -- see the design doc for the full argument).
# ---------------------------------------------------------------------
F_over_rho = tuple(sp.Function(f"F_{c_}")(x, y, z, t) for c_ in "xyz")
# Substitute F_V = rho*F_over_rho directly (plain definition, no frame claim)
subs_FV_to_F = {F_V[i]: rho * F_over_rho[i] for i in range(3)}

Du_Dt_derived_F = Du_Dt_derived.subs(subs_FV_to_F)
Du_Dt_derived_F = sp.simplify(Du_Dt_derived_F)

discretization_target = -div(tuple(rho * fi for fi in F_over_rho)) / rho - u / tau + S
discretization_target = sp.simplify(discretization_target)

leftover_expected = sp.simplify(div(tuple(u_V * vi for vi in v)) / rho)

print("=" * 78)
print("Claim 2: writing F_V = rho*F (F := the plain per-mass flux the")
print("code actually tracks, NOT a frame reinterpretation), the derived")
print("equation is EXACTLY the discretization's target PLUS one leftover")
print("advective term:")
print("  Du/Dt = -(1/rho)*div(rho*F) - u/tau + S + (1/rho)*div(u_V*v)")
print("=" * 78)
residual_claim2 = sp.simplify(
    Du_Dt_derived_F - (discretization_target + leftover_expected)
)
print("  residual (should be exactly 0):", residual_claim2)
assert residual_claim2 == 0, "Leftover-term decomposition failed"
print("PASS -- the leftover term is exactly (1/rho)*div(u_V*v), isolated")
print("cleanly, not absorbed into anything else.")
print()

# ---------------------------------------------------------------------
# Step 7: confirm the leftover term is NOT identically zero in general
# (i.e. this is a real physical omission, not a term that happens to
# vanish under the equations' own structure) -- check with a concrete,
# generic (non-uniform, non-static) choice of rho, v, u.
# ---------------------------------------------------------------------
concrete_rho = 1.0 + 0.3 * sp.sin(x) * sp.exp(-t)
concrete_u = 2.0 + 0.5 * sp.cos(y) * t
concrete_v = (0.4 * sp.sin(z), 0.1 * x, 0.0)
concrete_u_V = concrete_rho * concrete_u
leftover_concrete = sp.simplify(
    div(tuple(concrete_u_V * vi for vi in concrete_v)) / concrete_rho
)
print("=" * 78)
print("Step 7: leftover term is NOT identically zero for a generic field")
print("(sanity check that this is a real, non-vacuous finding):")
print("=" * 78)
print("  concrete leftover term:", leftover_concrete)
val_at_point = leftover_concrete.subs({x: 0.7, y: 1.1, z: 0.4, t: 0.9})
print("  numeric value at (x,y,z,t)=(0.7,1.1,0.4,0.9):", sp.N(val_at_point))
assert sp.N(val_at_point) != 0, "leftover term vanished at a generic point -- suspicious"
print("PASS: confirmed nonzero at a generic point -- a real term, not an")
print("algebraic artifact that happens to cancel.")
print()

# ---------------------------------------------------------------------
# Step 7b: decompose the leftover term to see WHICH piece SPH's own
# Lagrangian d/dt already provides for free (the "double-counting"
# question this task set out to resolve). div(u_V*v) = u_V*div(v) +
# v.grad(u_V) [product rule]; writing u_V = rho*u and dividing by rho:
#   (1/rho)*div(u_V*v) = u*div(v) + v.grad(u) + u*(v.grad(rho))/rho
# The middle piece, v.grad(u), is the Eulerian advection term an SPH
# particle's own material derivative provides by moving at v. The other
# two are the Reynolds-transport bookkeeping of a box-frame field in a
# gas-frame control volume; the gas-frame scheme drops all three (see the
# module docstring: this is the frame choice, not missing physics).
# ---------------------------------------------------------------------
leftover_full = div(tuple(u_V * vi for vi in v)) / rho
leftover_decomposed = u * div(v) + dot(v, grad(u)) + u * dot(v, grad(rho)) / rho
residual_decomposition = sp.simplify(leftover_full - leftover_decomposed)
print("=" * 78)
print("Step 7b: decompose the leftover term into its three pieces")
print("  (1/rho)*div(u_V*v) == u*div(v) + v.grad(u) + u*(v.grad(rho))/rho")
print("  (residual, should be exactly 0):")
print("=" * 78)
print("  ", residual_decomposition)
assert residual_decomposition == 0, "Leftover-term decomposition identity failed"
print("PASS: v.grad(u) is the piece SPH's own Lagrangian d/dt already")
print("provides (not to be double-counted); u*div(v) and")
print("u*(v.grad(rho))/rho are the box-frame bookkeeping pieces the")
print("gas-frame scheme deliberately does not carry (module docstring).")
print()

# ---------------------------------------------------------------------
# Step 8: the order-of-magnitude estimate quantifying how large the
# leftover term is relative to the transport term the discretization
# DOES compute, using the design's own steady-state scaling
# (F ~ -D*grad(u), D = c_hyp*lambda, so |F_V| ~ rho*c_hyp*u/L for some
# length scale L ~ lambda, i.e. |F_V| ~ c_hyp*u_V to order of magnitude).
# This is a scaling argument, not a symbolic identity -- reported as such.
# ---------------------------------------------------------------------
c_hyp_sym, v_gas_sym = sp.symbols("c_hyp v_gas", positive=True)
u_V_scale, L_scale = sp.symbols("u_V_scale L", positive=True)
F_V_scale_est = c_hyp_sym * u_V_scale  # |F_V| ~ c_hyp * u_V, from F ~ D*grad(u)/... ~ c_hyp*lambda*(u/lambda)
div_FV_scale_est = F_V_scale_est / L_scale
div_uVv_scale_est = u_V_scale * v_gas_sym / L_scale
ratio = sp.simplify(div_uVv_scale_est / div_FV_scale_est)
print("=" * 78)
print("Step 8: order-of-magnitude ratio, leftover term vs. the transport")
print("term the discretization already computes")
print("=" * 78)
print("  |div(u_V*v)| / |div(F_V)|  ~  ", ratio, " = v_gas / c_hyp")
print("  Here v_gas is the gas velocity RELATIVE TO THE FRAME the equation")
print("  is written in. In the box frame (galaxy orbital speeds, 100-200")
print("  km/s) this is 30-100 for c_hyp ~ c_s: the size of the error a")
print("  box-frame reduced-speed-of-light scheme carries in its neglected")
print("  O(v/c_hyp) mixed-frame terms, and why this design does NOT add the")
print("  term. In the gas frame the analogous ratio is v_rel/c_hyp with")
print("  v_rel the star-gas relative speed, a few km/s over c_hyp; that is")
print("  the residual the gas-frame scheme actually carries (design doc")
print("  1.0.1, validation legs 6.5/6.6).")
print()

# ---------------------------------------------------------------------
# Step 9: the §1.1/§1.2 steady-state reduction -- confirm the amplitude/
# screening-length results built on plain div(F) survive UNCHANGED under
# a LOCALLY UNIFORM density assumption (the setup §6.1/§6.2's own
# validation legs actually use -- a uniform box), and exhibit the extra
# drift term that appears once rho is allowed to vary.
# ---------------------------------------------------------------------
rho_const = sp.Symbol("rho0", positive=True)
D_sym = sp.Symbol("D", positive=True)
u_generic = sp.Function("u")(x, y, z)
F_generic = tuple(-D_sym * sp.diff(u_generic, coord) for coord in coords)  # F = -D*grad(u)

# General (non-uniform rho) reduction: (1/rho)*div(rho*F) = div(F) + F.grad(ln(rho))
rho_field_generic = sp.Function("rho")(x, y, z)
lhs_general = sp.diff(rho_field_generic * F_generic[0], x) / rho_field_generic \
    + sp.diff(rho_field_generic * F_generic[1], y) / rho_field_generic \
    + sp.diff(rho_field_generic * F_generic[2], z) / rho_field_generic
lhs_general = sp.simplify(lhs_general)
div_F_generic = div(F_generic)
grad_ln_rho = tuple(sp.diff(sp.log(rho_field_generic), coord) for coord in coords)
drift_term = dot(F_generic, grad_ln_rho)
rhs_general = sp.simplify(div_F_generic + drift_term)
residual_general = sp.simplify(lhs_general - rhs_general)
print("=" * 78)
print("Step 9a: general (non-uniform rho) reduction")
print("  (1/rho)*div(rho*F) == div(F) + F.grad(ln(rho))  (residual, should be 0):")
print("=" * 78)
print("  ", residual_general)
assert residual_general == 0

# Uniform-density special case: rebuild with rho as a genuine constant
# (derivative = 0, not a substitution into an already-differentiated
# expression), confirm the drift term vanishes and the operator reduces
# to plain div(F) exactly.
lhs_uniform = sp.diff(rho_const * F_generic[0], x) / rho_const \
    + sp.diff(rho_const * F_generic[1], y) / rho_const \
    + sp.diff(rho_const * F_generic[2], z) / rho_const
lhs_uniform = sp.simplify(lhs_uniform)
residual_uniform = sp.simplify(lhs_uniform - div_F_generic)
print()
print("Step 9b: uniform-density special case (rho = rho0, a constant):")
print("  (1/rho0)*div(rho0*F) == div(F)  (residual, should be 0):")
print("  ", residual_uniform)
assert residual_uniform == 0
print("PASS: the §1.1/§1.2 Yukawa/amplitude derivations (which assume a")
print("locally uniform ambient density around the point source -- exactly")
print("the setup §6.1/§6.2's own validation boxes use) are UNCHANGED by")
print("the mass-specific reformulation: (1/rho)*div(rho*F) reduces")
print("identically to plain div(F) there. Away from that assumption (a")
print("real density gradient), an extra drift term F.grad(ln rho) appears")
print("-- not a new problem introduced by this rework, since §1.2's own")
print("derivation already implicitly assumed a locally uniform medium for")
print("the point-source Green's function to apply at all.")
print()

print("=" * 78)
print("ALL CHECKS PASSED")
print("=" * 78)
print("Summary:")
print("  1. Du/Dt = -(1/rho)*div(F_V - u_V*v) - u/tau + S   [exact, general]")
print("  2. Writing F_V = rho*F (plain per-mass flux, no frame claim):")
print("     Du/Dt = -(1/rho)*div(rho*F) - u/tau + S + (1/rho)*div(u_V*v)")
print("     The discretization computes everything except the LAST term.")
print("  3. That leftover term is real (nonzero for a generic field) and")
print("     scales as v_gas/c_hyp relative to the transport term already")
print("     computed. It is the price of writing the reduced-speed-of-light")
print("     system in the BOX frame (v_gas = velocity relative to the box,")
print("     30-100 c_hyp in a galaxy); the gas-frame scheme drops it and")
print("     carries an O(v_rel/c_hyp) residual instead (design doc 1.0.1),")
print("     measured as a v_rel*tau centroid lag.")
print("  4. Sections 1.1/1.2's amplitude/screening-length results are")
print("     UNCHANGED under a locally uniform ambient density (exactly")
print("     the validation setup already used); a drift term")
print("     F.grad(ln rho) appears only once density is allowed to vary,")
print("     which the existing amplitude derivation never assumed away.")
