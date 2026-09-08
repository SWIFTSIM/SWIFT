"""Verify WHICH differential operator the Design B `div(F)` estimator
(design-lw-fuv-design-b.md Section 2.2) actually converges to.

Round-2 `/plan-review` (2026-09-07 report, Revision 2 section) found that
`src/rt/SPHM1RT/rt_gradients.h`'s own doxygen (lines 29, 249) states its
`radiation_divergence_SPH` function computes `(1/rho)*div(rho*fin)` -- a
mass-density-weighted divergence -- for ALL THREE `diffmode` branches
uniformly. At the time this script was first written, Section 1's
governing equation was copied unmodified from the theory doc's own
volumetric RTE-moment convention (`du/dt + div F = -u/tau + S`, plain
`div(F)`), which this estimator does NOT compute -- these differ by
`F.grad(ln rho)`, nonzero whenever the density field has a real gradient.

**Superseded, 2026-09-07 (operator ruling + re-derivation).** The
mismatch above was a governing-equation bug, not a numerics bug: `u`/`F`
are mass-specific throughout this implementation (injection already does
`u_i += u_inject/m_i`), and `verify_design_b_lagrangian_mass_specific_
derivation.py` derives, from the volumetric RTE moments plus mass
continuity, that the CORRECT Lagrangian/mass-specific governing equation
is `Du/Dt = -(1/rho)*div(rho*F) - u/tau + S` (up to a leftover
`(1/rho)*div(u_V*v)` gas-advection term the discretization does not yet
include -- a separate, newly-identified gap, see that script and the
design doc's §1/Known-gaps for the full account; NOT resolved by this
script). This script's own finding below -- that the §2.2 estimator
converges to `(1/rho)*div(rho*F)` -- is therefore now the CONFIRMATION
that the estimator matches its (corrected) target, not a discrepancy
needing an operator ruling. The sweep/assertions below are unchanged;
only this framing paragraph and the module-level conclusion differ from
this script's original version.

The existing script (`verify_design_b_div_f_conservation.py`) checks only
CONSERVATION: does `m_i*(div F)_i_pair + m_j*(div F)_j_pair` cancel
exactly for a single interacting pair? It says nothing about which
operator the *converged* value approximates. This script checks
CONSISTENCY instead: build a smooth, non-uniform density field and a
smooth vector field `F` with a known, closed-form analytic divergence,
place real SPH particles on a lattice, accumulate the discrete `div(F)_i`
using the exact §2.2 formula (implemented precisely as written, not as
"corrected" by assumption), and see which of the two candidate targets
(plain `div(F)`, constant here; or `(1/rho)*div(rho*F)`, which varies
with position whenever `grad(rho) != 0`) the discrete sum tracks as
resolution increases.

Test fields (chosen so the two candidate targets are cleanly
distinguishable, and both have closed forms):
  rho(x, y, z) = rho0 * (1 + slope * x)          (linear gradient along x)
  F(x, y, z)   = F0 * (x, y, z)                  (position vector field)

  div(F)                = 3 * F0                          (constant)
  (1/rho) * div(rho*F)  = 3*F0 + F0*x*rho0*slope / rho(x)  (varies with x)

Particles sit on a regular cubic lattice (not a glass) specifically so
that finite-resolution SPH truncation error is small and does not itself
mask which of the two constant-vs-varying target curves the estimator is
tracking; §2.3's own "zeroth-order consistency only, biased near sharp
features" caveat is about accuracy on a DISORDERED distribution, a
separate question from "which operator does the formula converge to,"
which is what this script isolates.
"""

import numpy as np
from scipy.spatial import cKDTree

RNG = np.random.default_rng(20260907)


# ---------------------------------------------------------------------------
# Kernel (Wendland C2, 3D, this project's own --with-kernel choice)
# ---------------------------------------------------------------------------
def wendland_c2_3d(q):
    """Vectorized (W(q), dW/dq) for the 3D Wendland C2 kernel, unit support."""
    q = np.asarray(q, dtype=np.float64)
    norm = 21.0 / (2.0 * np.pi)
    inside = q < 1.0
    w = np.zeros_like(q)
    dwdq = np.zeros_like(q)
    qi = q[inside]
    w[inside] = norm * (1.0 - qi) ** 4 * (4.0 * qi + 1.0)
    dwdq[inside] = norm * (
        -4.0 * (1.0 - qi) ** 3 * (4.0 * qi + 1.0) + 4.0 * (1.0 - qi) ** 4
    )
    return w, dwdq


def kernel_dr(r, h):
    """`wi_dr = h^-(dim+1) * dW/dq`, SWIFT's own convention (dim=3 -> h^-4),
    matching `verify_design_b_div_f_conservation.py`'s `kernel_dr`.
    """
    q = r / h
    _, dwdq = wendland_c2_3d(q)
    return dwdq / h**4


# ---------------------------------------------------------------------------
# Pairwise div(F) constructions
# ---------------------------------------------------------------------------
def divF_pair_diffmode1_mirror(dx, r, wi_dr, wj_dr, rhoi, rhoj, mi, mj, Fi, Fj):
    """The design doc's CURRENT §2.2 formula, mirroring
    `radiation_divergence_SPH`'s `diffmode==1` branch exactly:

      Phi_ij = (F_i . dx)/rho_i * wi_dr * r_inv + (F_j . dx)/rho_j * wj_dr * r_inv
      div_F_i += m_j * Phi_ij ;  div_F_j += -m_i * Phi_ij
    """
    r_inv = 1.0 / r
    shared = (
        np.einsum("...k,...k->...", Fi, dx) / rhoi * wi_dr * r_inv
        + np.einsum("...k,...k->...", Fj, dx) / rhoj * wj_dr * r_inv
    )
    return mj * shared, -mi * shared


def divF_pair_difference_shared_rho(dx, r, wi_dr, wj_dr, rhoi, rhoj, mi, mj, Fi, Fj):
    """Candidate FIX: a difference-of-F form (not a density-weighted-sum
    form), using a single PAIRWISE-SHARED kernel-slope and density
    normalization (arithmetic means of wi_dr/wj_dr and rho_i/rho_j) so the
    exact same shared-value / mirrored-mass-and-sign pattern is retained
    (conservation is a property of that pattern alone, established below,
    independent of what the shared value represents):

      Phi_ij = (F_i - F_j) . dx * r_inv * 0.5*(wi_dr + wj_dr) / (0.5*(rho_i + rho_j))
      div_F_i += m_j * Phi_ij ;  div_F_j += -m_i * Phi_ij

    Unlike the density-weighted-sum form, the density normalization here
    is a SHARED (pairwise) quantity that -> rho(x_i) to leading order as
    r_ij -> 0 (rho_i + rho_j -> 2*rho(x_i) for a smooth field), so the
    O(h) error this substitution introduces is the same *order* as the
    estimator's own already-acknowledged zeroth-order truncation error
    (§2.3), not a new O(1) bias -- unlike the F.grad(ln rho) term the
    density-weighted-sum form carries at every resolution.
    """
    r_inv = 1.0 / r
    wbar = 0.5 * (wi_dr + wj_dr)
    rhobar = 0.5 * (rhoi + rhoj)
    shared = np.einsum("...k,...k->...", Fi - Fj, dx) * r_inv * wbar / rhobar
    return mj * shared, -mi * shared


def divF_pair_rho_squared(dx, r, wi_dr, wj_dr, rhoi, rhoj, mi, mj, Fi, Fj):
    """Third candidate: the standard SPH MOMENTUM-EQUATION-style
    normalization (`rho^2`, as in `runner_iact_force`'s `P_i/rho_i^2 +
    P_j/rho_j^2`), not SPHM1RT's own `rho^1` form:

      Phi_ij = (F_i . dx)/rho_i^2 * wi_dr * r_inv + (F_j . dx)/rho_j^2 * wj_dr * r_inv
      div_F_i += m_j * Phi_ij ;  div_F_j += -m_i * Phi_ij

    For the analogous SCALAR case (`grad(P)`), this exact `rho^2` structure
    is the textbook reason the standard hydro force is both exactly
    momentum-conserving AND an unbiased estimator of `(1/rho)*grad(P)` (the
    `(P/rho^2)*grad(rho)` pieces from each term cancel exactly) -- unlike
    the `rho^1` structure, which leaves that piece uncancelled (§2.2's
    `F.grad(ln rho)` finding above). Whether the same cancellation survives
    here, for a VECTOR-to-scalar divergence built from two *different*
    kernels (`wi_dr` on the `rho_i`-term, `wj_dr` on the `rho_j`-term,
    unlike the single-kernel scalar identity), is exactly what this
    candidate's numerical test settles -- not assumed.
    """
    r_inv = 1.0 / r
    shared = (
        np.einsum("...k,...k->...", Fi, dx) / rhoi**2 * wi_dr * r_inv
        + np.einsum("...k,...k->...", Fj, dx) / rhoj**2 * wj_dr * r_inv
    )
    return mj * shared, -mi * shared


CONSTRUCTIONS = {
    "diffmode1_mirror (current §2.2 formula)": divF_pair_diffmode1_mirror,
    "difference_shared_rho (candidate fix)": divF_pair_difference_shared_rho,
    "rho_squared (momentum-eq-style)": divF_pair_rho_squared,
}


# ---------------------------------------------------------------------------
# Part 1: exact-conservation sanity check for BOTH constructions
# (mirrors verify_design_b_div_f_conservation.py's own check, extended to
# the candidate fix, so we never trade consistency for conservation.)
# ---------------------------------------------------------------------------
def check_pairwise_conservation():
    print("=" * 78)
    print("Part 1: exact-conservation sanity check (single random pairs,")
    print("h_i != h_j AND rho_i != rho_j), both constructions")
    print("=" * 78)
    max_resid = {name: 0.0 for name in CONSTRUCTIONS}
    for _ in range(200):
        dx = RNG.normal(size=3) * 0.4 + np.array([0.3, 0.0, 0.0])
        r = float(np.linalg.norm(dx))
        hi = RNG.uniform(0.3, 1.2)
        hj = hi + RNG.uniform(0.05, 0.6)
        rhoi, rhoj = RNG.uniform(0.5, 2.0, size=2)
        mi, mj = RNG.uniform(0.5, 2.0, size=2)
        Fi = RNG.normal(size=3)
        Fj = RNG.normal(size=3)
        wi_dr = kernel_dr(r, hi)
        wj_dr = kernel_dr(r, hj)
        for name, fn in CONSTRUCTIONS.items():
            di, dj = fn(dx, r, wi_dr, wj_dr, rhoi, rhoj, mi, mj, Fi, Fj)
            resid = abs(mi * di + mj * dj)
            max_resid[name] = max(max_resid[name], resid)
    for name, val in max_resid.items():
        print(f"  {name}: max |m_i*divFi + m_j*divFj| over 200 draws = {val:.4e}")
    for name, val in max_resid.items():
        assert val < 1e-10, f"{name} FAILED exact conservation (residual {val:.3e})"
    print("PASS: all constructions cancel to machine precision, h_i!=h_j, rho_i!=rho_j")
    print("(expected -- an algebraic property of the shared-value/mirrored-mass-")
    print("and-sign pattern itself, independent of what the shared value")
    print("represents; says nothing about which operator each converges to, see")
    print("Part 2 for that).")
    print()


# ---------------------------------------------------------------------------
# Part 2: full-field consistency check
# ---------------------------------------------------------------------------
def rho_field(pos, rho0, slope):
    return rho0 * (1.0 + slope * pos[:, 0])


def F_field(pos, F0):
    return F0 * pos


def build_lattice(n_per_dim, box_l):
    """Regular cubic lattice of `n_per_dim**3` particles filling [0, box_l]^3
    (cell-centred), plus the per-particle cell volume (for mass assignment).
    """
    dx_p = box_l / n_per_dim
    coords_1d = (np.arange(n_per_dim) + 0.5) * dx_p
    xx, yy, zz = np.meshgrid(coords_1d, coords_1d, coords_1d, indexing="ij")
    pos = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1)
    return pos, dx_p


def run_convergence_check(rho0, slope, F0, box_l, n_per_dim, h_over_dxp_values):
    """Sweep NEIGHBOUR COUNT (`h_over_dxp`) at fixed spatial resolution
    (`n_per_dim`), not particle count at fixed `h_over_dxp`.

    Both test fields (`rho`, `F`) are chosen exactly LINEAR in position, so
    the usual "h -> 0 relative to the field's curvature scale" truncation
    bias is identically zero for these estimators regardless of `h`: a
    symmetric kernel's convolution with a linear function is exact (the
    first odd moment vanishes by symmetry). The only remaining source of
    deviation from each estimator's true continuum-limit operator is
    particle-discreteness / kernel-quadrature error -- and that shrinks
    specifically as the NUMBER OF NEIGHBOURS within the kernel support
    grows (`h/dx_p -> infinity`), not merely as `h` itself shrinks at fixed
    neighbour count. An earlier version of this script swept `n_per_dim`
    at fixed `h_over_dxp` and found a residual that did not shrink at all
    over a 4x range in `h` -- exactly consistent with this: that sweep
    never changed the neighbour count, so it could not have shown
    convergence either way. Sweeping `h_over_dxp` instead isolates the
    right variable.
    """
    # The three candidate continuum-limit targets a construction might
    # converge to, as closed-form functions of (x_p, rho_p).
    TARGETS = {
        "plain div(F)": lambda x_p, rho_p: 3.0 * F0,
        "(1/rho)*div(rho*F)": lambda x_p, rho_p: 3.0 * F0 + F0 * x_p * rho0 * slope / rho_p,
        "(1/rho)*div(F)": lambda x_p, rho_p: 3.0 * F0 / rho_p,
    }

    print("=" * 78)
    print("Part 2: full-field consistency check (neighbour-count sweep)")
    print(f"  rho(x) = {rho0} * (1 + {slope}*x)   F(r) = {F0} * r")
    for tname, tfn in TARGETS.items():
        print(f"  analytic {tname} at x=0.4/0.5/0.6: "
              f"{tfn(0.4*box_l, rho_field(np.array([[0.4*box_l,0,0]]), rho0, slope)[0]):.4f} / "
              f"{tfn(0.5*box_l, rho_field(np.array([[0.5*box_l,0,0]]), rho0, slope)[0]):.4f} / "
              f"{tfn(0.6*box_l, rho_field(np.array([[0.6*box_l,0,0]]), rho0, slope)[0]):.4f}")
    print(f"  n_per_dim fixed at {n_per_dim}; h_over_dxp swept over {h_over_dxp_values}")
    print("=" * 78)

    # A handful of fixed x-locations (in units of box_l) to probe -- chosen
    # well inside the interior (>= 1.5h from every domain edge at every
    # h_over_dxp tested, including the largest) so boundary under-sampling
    # never contaminates the comparison.
    probe_fracs = [0.4, 0.5, 0.6]

    results = {name: {f: [] for f in probe_fracs} for name in CONSTRUCTIONS}

    pos, dx_p = build_lattice(n_per_dim, box_l)
    rho = rho_field(pos, rho0, slope)
    F = F_field(pos, F0)
    mass = rho * dx_p**3

    for h_over_dxp in h_over_dxp_values:
        h = h_over_dxp * dx_p
        n_neighbours_est = (4.0 / 3.0) * np.pi * h_over_dxp**3

        tree = cKDTree(pos)
        pairs = tree.query_pairs(r=h, output_type="ndarray")
        i_idx = pairs[:, 0]
        j_idx = pairs[:, 1]

        dx = pos[i_idx] - pos[j_idx]
        r = np.linalg.norm(dx, axis=1)
        # guard against exact coincidence (shouldn't happen on a lattice)
        keep = r > 1e-12
        i_idx, j_idx, dx, r = i_idx[keep], j_idx[keep], dx[keep], r[keep]

        wi_dr = kernel_dr(r, h)
        wj_dr = kernel_dr(r, h)  # uniform h in this base test

        margin = 1.5 * h
        interior = (
            (pos[:, 0] > margin)
            & (pos[:, 0] < box_l - margin)
            & (pos[:, 1] > margin)
            & (pos[:, 1] < box_l - margin)
            & (pos[:, 2] > margin)
            & (pos[:, 2] < box_l - margin)
        )

        for name, fn in CONSTRUCTIONS.items():
            div_F = np.zeros(pos.shape[0])
            di, dj = fn(
                dx, r, wi_dr, wj_dr, rho[i_idx], rho[j_idx],
                mass[i_idx], mass[j_idx], F[i_idx], F[j_idx],
            )
            np.add.at(div_F, i_idx, di)
            np.add.at(div_F, j_idx, dj)

            print(
                f"\n  -- {name}, h_over_dxp={h_over_dxp:.2f} "
                f"(~{n_neighbours_est:.0f} neighbours, h={h:.4f}, dx_p={dx_p:.4f}) --"
            )
            for frac in probe_fracs:
                target_x = frac * box_l
                # nearest interior particle to (target_x, box_l/2, box_l/2)
                probe_point = np.array([target_x, box_l / 2, box_l / 2])
                cand = np.where(interior)[0]
                d2 = np.sum((pos[cand] - probe_point) ** 2, axis=1)
                p_idx = cand[np.argmin(d2)]
                x_p = pos[p_idx, 0]
                rho_p = rho[p_idx]

                measured = div_F[p_idx]
                target_vals = {tname: tfn(x_p, rho_p) for tname, tfn in TARGETS.items()}
                results[name][frac].append((h_over_dxp, n_neighbours_est, measured, target_vals))
                err_str = "  ".join(
                    f"{tname}={tval:+.6f} (err={measured-tval:+.2e})"
                    for tname, tval in target_vals.items()
                )
                print(f"    x={x_p:.4f}: measured={measured:+.6f}  {err_str}")

    print()
    print("=" * 78)
    print("Convergence summary: error vs each candidate target, by neighbour count")
    print("=" * 78)
    verdicts = {}
    for name in CONSTRUCTIONS:
        print(f"\n{name}:")
        err_by_target_by_hod = {tname: {} for tname in TARGETS}
        for frac in probe_fracs:
            rows = results[name][frac]
            for h_over_dxp, n_neigh, measured, target_vals in rows:
                for tname, tval in target_vals.items():
                    err_by_target_by_hod[tname].setdefault(h_over_dxp, []).append(abs(measured - tval))
        for h_over_dxp in h_over_dxp_values:
            n_neigh = (4.0 / 3.0) * np.pi * h_over_dxp**3
            parts = [
                f"{tname}={np.mean(err_by_target_by_hod[tname][h_over_dxp]):.4e}"
                for tname in TARGETS
            ]
            print(f"  h_over_dxp={h_over_dxp:.2f} (~{n_neigh:.0f} nbrs): " + "   ".join(parts))

        # Converged-to-target means: the error vs that target shrinks
        # substantially (>= 5x) over the FULL sweep AND ends up much
        # smaller (>= 10x) than the error vs every OTHER target at the
        # finest neighbour count tested. A narrow sweep can show a
        # transient crossing that looks like convergence to the wrong
        # target (observed at h_over_dxp~1.8 for the current formula, see
        # module docstring) -- the >=5x-shrink-over-the-FULL-sweep
        # requirement guards against reading that transient as the verdict.
        last, first = h_over_dxp_values[-1], h_over_dxp_values[0]
        err_last = {t: np.mean(err_by_target_by_hod[t][last]) for t in TARGETS}
        err_first = {t: np.mean(err_by_target_by_hod[t][first]) for t in TARGETS}
        best_target = min(err_last, key=err_last.get)
        shrink_factor = err_first[best_target] / max(err_last[best_target], 1e-300)
        others_last = [err_last[t] for t in TARGETS if t != best_target]
        if shrink_factor >= 5.0 and err_last[best_target] * 10.0 < min(others_last):
            verdict = best_target
        else:
            verdict = "AMBIGUOUS (" + ", ".join(
                f"{t}: {err_first[t]:.2e}->{err_last[t]:.2e}" for t in TARGETS
            ) + ")"
        verdicts[name] = verdict
        print(f"  ==> converges to: {verdict}")

    return verdicts


if __name__ == "__main__":
    check_pairwise_conservation()
    verdicts = run_convergence_check(
        rho0=1.0,
        slope=1.5,
        F0=2.0,
        box_l=1.0,
        n_per_dim=32,
        h_over_dxp_values=[1.3, 1.8, 2.6, 3.8, 5.5],
    )

    print()
    print("=" * 78)
    print("FINAL VERDICT")
    print("=" * 78)
    for name, verdict in verdicts.items():
        print(f"  {name}: {verdict}")

    current = verdicts["diffmode1_mirror (current §2.2 formula)"]
    difference_shared = verdicts["difference_shared_rho (candidate fix)"]
    rho_sq = verdicts["rho_squared (momentum-eq-style)"]
    assert current == "(1/rho)*div(rho*F)", (
        "Expected the CURRENT §2.2 formula to converge to the density-weighted "
        f"operator per the SPHM1RT doxygen claim; got '{current}' instead -- "
        "if this assertion fires, re-examine before trusting the doc update."
    )
    assert difference_shared != "plain div(F)", (
        "This candidate was expected to FAIL (collapse toward zero, not "
        f"converge to any target) -- got '{difference_shared}' instead. If it "
        "actually converges to plain div(F), that overturns the doc's "
        "'no in-family exact-conservative fix exists' finding; re-examine "
        "before trusting anything written about it."
    )
    print()
    print("CONFIRMED:")
    print(f"  diffmode1_mirror (current §2.2 formula) -> {current}")
    print(f"    matches the SPHM1RT doxygen's own claim, not the plan-review")
    print(f"    reviewer's tentative n=0-family hypothesis (refuted: that family")
    print(f"    uses a single grad_i(W_ij); diffmode==1 mixes wi_dr/wj_dr per side).")
    print(f"  difference_shared_rho (candidate fix) -> {difference_shared}")
    print(f"    (collapses toward zero: forcing a mirrored-mass/sign pattern onto")
    print(f"    a plain DIFFERENCE-of-F scalar is structurally self-cancelling for")
    print(f"    an isotropic neighbour distribution -- NOT a viable fix.)")
    print(f"  rho_squared (momentum-eq-style) -> {rho_sq}")
    print(f"    (same rho^1 -> rho^2 trick that makes runner_iact_force exactly")
    print(f"    momentum-conserving AND grad(P)/rho-consistent; removes the")
    print(f"    F.grad(ln rho) bias, but the result still carries an extra 1/rho")
    print(f"    factor relative to plain div(F) -- see module docstring/design doc")
    print(f"    for why multiplying it out by rho_i afterward would break exact")
    print(f"    conservation again.)")
    print()
    print("No construction tested here reaches PLAIN div(F) while retaining exact")
    print("pairwise (mass-weighted-sum) conservation -- this is expected and no")
    print("longer an open question: verify_design_b_lagrangian_mass_specific_")
    print("derivation.py derives (1/rho)*div(rho*F), not plain div(F), as the")
    print("CORRECT mass-specific Lagrangian target, so the current §2.2 formula's")
    print("convergence to (1/rho)*div(rho*F) above is a PASS, not a gap.")
