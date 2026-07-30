/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include <config.h>

#if defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) || \
    defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)

/* Local headers. */
#include "swift.h"

/* System headers. */
#include <stdio.h>
#include <string.h>

/**
 * @file testChemistryFluxAntisymmetry.c
 * @brief Gate test for the chemistry diffusion flux-exchange dispatcher
 * (runner_iact_chemistry_flux_exchange(), chemistry_iact.h).
 *
 * Four things are exercised, each a distinct arm of the dispatcher's decision
 * table; they are reported separately and must not be conflated:
 *
 *  - GATE (both_updatable_here=1, the "sym" call sites): a single evaluation
 *    updates both particles; checks that this evaluation is internally
 *    conservative (Fi == -Fj) and that it does not depend on which physical
 *    particle a caller happens to pass first (relevant since several call
 *    sites place the cell-j particle first).
 *  - MINDT ORDER-INDEPENDENCE: one particle active, one inactive. Directly
 *    exercises chemistry_compute_pair_fluxes()'s mindt computation
 *    (the three sites fixed this session) for the two argument bindings of
 *    the same physical pair, independent of the dispatcher's arm selection.
 *  - LOCAL MIXED-BAND TIE-BREAK (both_updatable_here=0, both local): the
 *    dispatcher's own two-sided visitation is simulated by calling it once
 *    per argument order; exactly one of the two must act (the position
 *    winner), the other must skip, and ties must skip on both sides.
 *  - CROSS-RANK (both_updatable_here=0, only the first argument local): the
 *    canonical-order computation must always apply only to the local slot,
 *    regardless of which side of the pair happens to be "first" here, and a
 *    tie must skip on both sides.
 */

/* ------------------------------------------------------------------------ */
/* Test scaffolding                                                         */
/* ------------------------------------------------------------------------ */

enum test_geometry { GEOM_WELL_BEHAVED, GEOM_ILL_CONDITIONED };

enum test_limiter_regime {
  LIMITER_INERT,
  LIMITER_BITING,
  LIMITER_NEAR_NOISE_GATE,
  LIMITER_NEAR_CAPACITY_POS,
  LIMITER_NEAR_CAPACITY_NEG
};

static const char *geom_name(enum test_geometry g) {
  return g == GEOM_WELL_BEHAVED ? "well-behaved" : "ill-conditioned";
}

static const char *limiter_name(enum test_limiter_regime l) {
  switch (l) {
    case LIMITER_INERT:
      return "inert";
    case LIMITER_BITING:
      return "biting";
    case LIMITER_NEAR_NOISE_GATE:
      return "near-noise-gate";
    case LIMITER_NEAR_CAPACITY_POS:
      return "near-capacity(+)";
    case LIMITER_NEAR_CAPACITY_NEG:
      return "near-capacity(-)";
  }
  return "?";
}

/* The metal specie index used throughout: only element 0 is driven with a
 * non-trivial state, the others are populated with distinct, non-degenerate
 * values so a bug that only shows up for m > 0 (e.g. an index accidentally
 * hardcoded to 0 in a fix) would still be caught by the max-over-elements
 * residual computed below. */
#define TEST_METAL 0

/**
 * @brief Fill a #part with a self-consistent, non-degenerate baseline state.
 *
 * Position defaults to the origin; callers that exercise the position
 * tie-break (flux_exchange_precedes()) must set p->x[0] explicitly
 * afterwards.
 * Left at zero here: velocity gradients, flux gradients (except where
 * run_config() drives them non-zero for the hyperbolic MUSCL term), and the
 * cell slope-limiter bounds are wide open. run_config() drives a genuine
 * jump in metal_mass between pi and pj (the mixed-depth bug's trigger) plus
 * a modest metal-density gradient (gradients.rhoZ, required for the
 * parabolic build to produce any flux at all -- see run_config()).
 */
static void init_part(struct part *p, long long id, float h, float mass,
                      float rho, double metal_mass_frac) {
  bzero(p, sizeof(struct part));

  p->id = id;
  p->x[0] = 0.0;
  p->x[1] = 0.0;
  p->x[2] = 0.0;
  p->v[0] = 0.0;
  p->v[1] = 0.0;
  p->v[2] = 0.0;
  p->h = h;
  p->mass = mass;
  p->rho = rho;
  p->u = 1.0f;
  p->time_bin = 1;

  /* Geometry: well-behaved identity-like matrix by default; per-test
   * configuration below may force the ill-conditioned fallback via wcorr. */
  p->geometry.volume = 1.0f;
  p->geometry.matrix_E[0][0] = 1.0f;
  p->geometry.matrix_E[1][1] = 1.0f;
  p->geometry.matrix_E[2][2] = 1.0f;
  p->geometry.wcorr = 1.0f;
  p->geometry.condition_number = 1.0f;

  struct chemistry_part_data *ch = &p->chemistry_data;
  ch->kappa = 1.0;
#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  ch->tau = 1.0;
#endif
  ch->limiter.maxr = 1.0f;
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    /* Wide-open cell limiter bounds: with zero gradients the cell limiter
     * takes its gradtrue == 0 fast path (alpha = 1) regardless, but keep
     * these sane in case a future change makes gradients non-zero. */
    ch->limiter.rhoZ[m][0] = -1e10;
    ch->limiter.rhoZ[m][1] = 1e10;
    ch->limiter.Z[m][0] = -1e10;
    ch->limiter.Z[m][1] = 1e10;
    ch->fct_theta[m] = 1.0; /* FCT inert unless overridden by the test */
    /* Distinct, non-degenerate metal masses across elements so an
     * accidentally-hardcoded element index would still be caught. */
    ch->metal_mass[m] = mass * (metal_mass_frac + 0.001 * m);
    if (ch->metal_mass[m] < 0.0) ch->metal_mass[m] = 0.0;
    if (ch->metal_mass[m] > 0.999 * mass) ch->metal_mass[m] = 0.999 * mass;
  }
  /* TEST_METAL gets exactly the requested fraction (no per-element jitter),
   * since that is the element the residual is computed on. */
  ch->metal_mass[TEST_METAL] = mass * metal_mass_frac;
}

static void set_geometry(struct part *p, enum test_geometry g) {
  if (g == GEOM_WELL_BEHAVED) {
    p->geometry.wcorr = 1.0f;
  } else {
    /* fvpm_part_geometry_well_behaved() returns
     * (p->geometry.wcorr > const_gizmo_min_wcorr); forcing it below that
     * threshold on either particle routes chemistry.c's area computation
     * through the SPH-face-area fallback, which is the branch that
     * actually fires in dense fragmented regions (see the plan, step 1). */
    p->geometry.wcorr = 0.0f;
  }
}

/* Reset only the fields that accumulate across a call: the flux
 * contributions and the FCT prep-pass outflow sums. Leaves the particle's
 * physical state (masses, geometry, kappa, tau, position, ...) untouched so
 * repeated calls on the same template see the identical input state. */
static void reset_flux_accumulators(struct part *p) {
  struct chemistry_part_data *ch = &p->chemistry_data;
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    ch->flux.metal_mass[m] = 0.0;
    ch->fct_sum_out[m] = 0.0;
#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
    ch->flux.flux[m][0] = 0.0;
    ch->flux.flux[m][1] = 0.0;
    ch->flux.flux[m][2] = 0.0;
#endif
  }
}

/* ------------------------------------------------------------------------ */
/* Result bookkeeping                                                       */
/* ------------------------------------------------------------------------ */

struct residual_result {
  double signed_residual; /* (Fi + Fj) / max(|Fi|, |Fj|) */
  double Fi;
  double Fj;
};

/**
 * @brief Evaluate the antisymmetry residual of the GATE (both_updatable_here)
 * arm for a single configured pair.
 *
 * Calls runner_iact_chemistry_flux_exchange(pi, pj, both_updatable_here=1)
 * and separately (pj, pi, -dx, both_updatable_here=1) -- both calls take the
 * "today's sym path, no reordering" arm (see the dispatcher's decision
 * table), so this measures the same thing the pre-fix gate test did: how
 * much the flux computation's floating-point result depends on which
 * physical particle a caller happens to bind to the function's first
 * parameter.
 */
static struct residual_result measure_gate_residual(
    struct part *pi_template, struct part *pj_template, const float dx[3],
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data) {
  struct part pi, pj;
  memcpy(&pi, pi_template, sizeof(struct part));
  memcpy(&pj, pj_template, sizeof(struct part));

  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float time_base = 1e-5;
  const integertime_t ti_current = 1;

  reset_flux_accumulators(&pi);
  reset_flux_accumulators(&pj);
  runner_iact_chemistry_flux_exchange(
      r2, dx, pi.h, pj.h, &pi, &pj, /*both_updatable_here=*/1,
      /*local_first=*/1, /*local_second=*/1, /*a=*/1.f,
      /*H=*/0.f, time_base, ti_current, cosmo, with_cosmology, chem_data);
  const double Fi = pi.chemistry_data.flux.metal_mass[TEST_METAL];

  memcpy(&pi, pi_template, sizeof(struct part));
  memcpy(&pj, pj_template, sizeof(struct part));
  const float mdx[3] = {-dx[0], -dx[1], -dx[2]};
  reset_flux_accumulators(&pi);
  reset_flux_accumulators(&pj);
  runner_iact_chemistry_flux_exchange(
      r2, mdx, pj.h, pi.h, &pj, &pi, /*both_updatable_here=*/1,
      /*local_first=*/1, /*local_second=*/1, /*a=*/1.f,
      /*H=*/0.f, time_base, ti_current, cosmo, with_cosmology, chem_data);
  const double Fj = pj.chemistry_data.flux.metal_mass[TEST_METAL];

  struct residual_result res;
  res.Fi = Fi;
  res.Fj = Fj;
  const double denom = max(fabs(Fi), fabs(Fj));
  res.signed_residual = (denom > 0.0) ? (Fi + Fj) / denom : 0.0;
  return res;
}

/**
 * @brief Evaluate order-independence of
 * chemistry_compute_pair_fluxes()'s mindt for a mixed-activity
 * pair (pi active, pj inactive): calls it with (pi, pj, dx) and separately
 * with (pj, pi, -dx) and compares mindt. Bypasses the dispatcher entirely --
 * a real call site never presents an inactive particle as the "first"
 * (entitled) argument, so this targets the three mindt sites directly
 * rather than the dispatcher's arm selection.
 */
static double measure_mindt_residual(
    struct part *pi_template, struct part *pj_template, const float dx[3],
    const struct cosmology *cosmo,
    const struct chemistry_global_data *chem_data) {
  struct part pi, pj;
  memcpy(&pi, pi_template, sizeof(struct part));
  memcpy(&pj, pj_template, sizeof(struct part));
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  double totflux[GEAR_CHEMISTRY_ELEMENT_COUNT][4];
  float mindt_fwd = 0.f, vmax = 0.f;
#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  float delxbar_i = 0.f, delxbar_j = 0.f;
  chemistry_compute_pair_fluxes(r2, dx, pi.h, pj.h, &pi, &pj, chem_data, cosmo,
                                totflux, &mindt_fwd, &vmax, &delxbar_i,
                                &delxbar_j);
#else
  chemistry_compute_pair_fluxes(r2, dx, pi.h, pj.h, &pi, &pj, chem_data, cosmo,
                                totflux, &mindt_fwd, &vmax);
#endif

  const float mdx[3] = {-dx[0], -dx[1], -dx[2]};
  float mindt_bwd = 0.f;
#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  chemistry_compute_pair_fluxes(r2, mdx, pj.h, pi.h, &pj, &pi, chem_data, cosmo,
                                totflux, &mindt_bwd, &vmax, &delxbar_i,
                                &delxbar_j);
#else
  chemistry_compute_pair_fluxes(r2, mdx, pj.h, pi.h, &pj, &pi, chem_data, cosmo,
                                totflux, &mindt_bwd, &vmax);
#endif

  const double denom = max(fabs(mindt_fwd), fabs(mindt_bwd));
  return (denom > 0.0) ? (mindt_fwd - mindt_bwd) / denom : 0.0;
}

/* ------------------------------------------------------------------------ */
/* Configuration builder                                                    */
/* ------------------------------------------------------------------------ */

static void build_chem_data(struct chemistry_global_data *chem_data,
                            double noise_gate) {
  bzero(chem_data, sizeof(struct chemistry_global_data));
  chem_data->diffusion_coefficient = 1.0f;
  chem_data->diffusion_mode = isotropic_constant;
#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  chem_data->tau = 1.0;
  chem_data->relaxation_time_mode = constant_mode;
  chem_data->hyperbolic_limiter_scope = limiter_all_components;
#endif
  chem_data->hll_riemann_solver_psi = 0.1f;
  chem_data->hll_riemann_solver_epsilon = 0.0f;
  chem_data->C_CFL_chemistry = 0.5f;
  chem_data->flux_limiter_noise_gate = noise_gate;
  chem_data->flux_limiter_safety = 0.5;
  chem_data->flux_limiter_sink_stability = 0.25;
  chem_data->flux_limiter_startup = 0.1;
}

/**
 * @brief Build one pair configuration (particles at the origin; dx supplies
 * the synthetic separation) and measure the GATE residual and the mindt
 * order-independence residual.
 *
 * @param mixed_activity if non-zero, pi is left active (flux.dt > 0) and pj
 * inactive (flux.dt < 0); the mindt check is only meaningful in this case
 * (the GATE arm doesn't touch flux.dt at all).
 */
static void run_config(const char *label, enum test_geometry geom,
                       float h_ratio, enum test_limiter_regime limiter_regime,
                       int hyp_flux_nonzero, int fct_enabled,
                       int mixed_activity, int with_cosmology,
                       struct residual_result *gate_out, double *mindt_out) {
  const float hi = 1.0f;
  const float hj = hi * h_ratio;

  /* Baseline masses tuned so the limiter regimes land where intended; see
   * the comment for each branch. Values found by inspection of
   * chemistry_limit_metal_mass_flux()'s constraints (chemistry_flux.h). */
  float mass_i = 1.0f, mass_j = 1.0f;
  double frac_i = 0.5, frac_j = 0.1;
  double noise_gate = 1e-15; /* production default */

  switch (limiter_regime) {
    case LIMITER_INERT:
      mass_i = 100.0f;
      mass_j = 100.0f;
      frac_i = 0.20;
      frac_j = 0.19;
      break;
    case LIMITER_BITING:
      mass_i = 1.0f;
      mass_j = 1.0f;
      frac_i = 0.99;
      frac_j = 1e-6;
      break;
    case LIMITER_NEAR_NOISE_GATE:
      noise_gate = 0.35;
      mass_i = 10.0f;
      mass_j = 10.0f;
      frac_i = 0.50;
      frac_j = 0.499999;
      break;
    case LIMITER_NEAR_CAPACITY_POS:
      mass_i = 1.0f;
      mass_j = 1.0f;
      frac_i = 0.999999;
      frac_j = 0.999999999;
      break;
    case LIMITER_NEAR_CAPACITY_NEG:
      mass_i = 1.0f;
      mass_j = 1.0f;
      frac_i = 0.999999999;
      frac_j = 0.999999;
      break;
  }

  struct part pi, pj;
  init_part(&pi, /*id=*/100, hi, mass_i, /*rho=*/1.0f, frac_i);
  init_part(&pj, /*id=*/200, hj, mass_j, /*rho=*/1.0f, frac_j);
  set_geometry(&pi, geom);
  set_geometry(&pj, geom);

  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    pi.chemistry_data.gradients.rhoZ[m][0] = 0.05;
    pj.chemistry_data.gradients.rhoZ[m][0] = -0.03;
  }

  if (fct_enabled) {
    for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
      pi.chemistry_data.fct_theta[m] = 0.3;
      pj.chemistry_data.fct_theta[m] = 0.3;
    }
  }

#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  if (hyp_flux_nonzero) {
    for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
      pi.chemistry_data.diffusion_flux[m][0] = 0.05;
      pi.chemistry_data.diffusion_flux[m][1] = -0.02;
      pi.chemistry_data.diffusion_flux[m][2] = 0.01;
      pj.chemistry_data.diffusion_flux[m][0] = -0.03;
      pj.chemistry_data.diffusion_flux[m][1] = 0.04;
      pj.chemistry_data.diffusion_flux[m][2] = -0.01;
      pi.chemistry_data.gradients.flux[m][0][0] = 0.02;
      pi.chemistry_data.gradients.flux[m][1][1] = -0.01;
      pi.chemistry_data.gradients.flux[m][2][2] = 0.005;
      pj.chemistry_data.gradients.flux[m][0][0] = -0.015;
      pj.chemistry_data.gradients.flux[m][1][1] = 0.02;
      pj.chemistry_data.gradients.flux[m][2][2] = -0.01;
    }
  }
#else
  (void)hyp_flux_nonzero;
#endif

  if (mixed_activity) {
    pi.chemistry_data.flux.dt = 1.0f;
    pj.chemistry_data.flux.dt = -1.0f;
  } else {
    pi.chemistry_data.flux.dt = 0.7f;
    pj.chemistry_data.flux.dt = 1.3f;
  }

  struct chemistry_global_data chem_data;
  build_chem_data(&chem_data, noise_gate);

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);

  const float dx[3] = {0.3f * hi, 0.0f, 0.0f};

  if (gate_out)
    *gate_out =
        measure_gate_residual(&pi, &pj, dx, &cosmo, with_cosmology, &chem_data);
  if (mindt_out)
    *mindt_out = measure_mindt_residual(&pi, &pj, dx, &cosmo, &chem_data);

  if (gate_out)
    message("  %-40s Fi=% .8e Fj=% .8e gate_residual=% .3e", label,
            gate_out->Fi, gate_out->Fj, gate_out->signed_residual);
  if (mindt_out) message("  %-40s mindt_residual=% .3e", label, *mindt_out);
}

/* ------------------------------------------------------------------------ */
/* Local mixed-band tie-break (arm 3) and cross-rank (arm 4) checks          */
/* ------------------------------------------------------------------------ */

/**
 * @brief Build a simple, well-behaved both-active pair and drive it through
 * the dispatcher's non-GATE arms.
 *
 * The position used for the tie-break (flux_exchange_precedes()) is
 * deliberately decoupled from the physical separation used for the flux
 * geometry (always a fixed dx = {0.3, 0, 0} in the caller): a real face's flux
 * can vanish for a particular dx sign under a pure-gradient (parabolic) driver
 * with a fixed gradient pair, which would make a naive "flip pi_x, keep dx =
 * pi.x - pj.x" setup accidentally test a degenerate face instead of the
 * tie-break logic. Using an arbitrary, dx-independent x-offset keeps every
 * configuration's face physically identical and non-degenerate; only the
 * tie-break outcome changes.
 *
 * @param order_sign +1 => flux_exchange_precedes(pi->x, pj->x) is FALSE (pj
 * "less"); -1 => TRUE (pi "less"); 0 => coincident-position tie.
 */
static void build_tiebreak_pair(struct part *pi, struct part *pj,
                                float order_sign) {
  init_part(pi, /*id=*/1, /*h=*/1.0f, /*mass=*/10.0f, /*rho=*/1.0f,
            /*metal_mass_frac=*/0.3);
  init_part(pj, /*id=*/2, /*h=*/1.0f, /*mass=*/10.0f, /*rho=*/1.0f,
            /*metal_mass_frac=*/0.1);
  pj->x[0] = 100.0;
  pi->x[0] = 100.0 + order_sign * 0.1;
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    pi->chemistry_data.gradients.rhoZ[m][0] = 0.05;
    pj->chemistry_data.gradients.rhoZ[m][0] = -0.03;
  }
  pi->chemistry_data.flux.dt = 0.7f;
  pj->chemistry_data.flux.dt = 1.3f;
}

/**
 * @brief Build a pair with hi != hj and the position order deliberately set
 * OPPOSITE of what h dictates, isolating h-primary designation (A1) from
 * the position fallback exercised by build_tiebreak_pair()/hi==hj above.
 *
 * @param hi_larger 1 => pi has the larger h; 0 => pj does. Either way,
 * position alone would pick the OTHER side, so a pass here requires h to
 * actually override position.
 */
static void build_h_primary_pair(struct part *pi, struct part *pj,
                                 int hi_larger) {
  const float h_big = 2.0f, h_small = 1.0f;
  init_part(pi, /*id=*/1, hi_larger ? h_big : h_small, /*mass=*/10.0f,
            /*rho=*/1.0f, /*metal_mass_frac=*/0.3);
  init_part(pj, /*id=*/2, hi_larger ? h_small : h_big, /*mass=*/10.0f,
            /*rho=*/1.0f, /*metal_mass_frac=*/0.1);
  pj->x[0] = 100.0;
  pi->x[0] = hi_larger ? 100.1 : 99.9;
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    pi->chemistry_data.gradients.rhoZ[m][0] = 0.05;
    pj->chemistry_data.gradients.rhoZ[m][0] = -0.03;
  }
  pi->chemistry_data.flux.dt = 0.7f;
  pj->chemistry_data.flux.dt = 1.3f;
}

/**
 * @brief Check that the local mixed-band arm's ownership is h-primary (A1):
 * with position deliberately reversed from what h dictates, the larger-h
 * side must still be the one that acts.
 *
 * @return 1 on pass, 0 on failure.
 */
static int check_h_primary_designation(
    const struct cosmology *cosmo,
    const struct chemistry_global_data *chem_data) {
  const float time_base = 1e-5;
  const integertime_t ti_current = 1;
  int all_ok = 1;

  for (int hi_larger = 0; hi_larger < 2; hi_larger++) {
    struct part pi, pj;
    build_h_primary_pair(&pi, &pj, hi_larger);
    const float dx[3] = {0.3f, 0.f, 0.f};
    const float r2 = dx[0] * dx[0];

    reset_flux_accumulators(&pi);
    reset_flux_accumulators(&pj);
    runner_iact_chemistry_flux_exchange(
        r2, dx, pi.h, pj.h, &pi, &pj, /*both_updatable_here=*/0,
        /*local_first=*/1, /*local_second=*/1, /*a=*/1.f, /*H=*/0.f, time_base,
        ti_current, cosmo, /*with_cosmology=*/0, chem_data);
    const int pi_acted = pi.chemistry_data.flux.metal_mass[TEST_METAL] != 0.0;
    const int pj_acted = pj.chemistry_data.flux.metal_mass[TEST_METAL] != 0.0;
    /* pi is entitled at this call site; it must act iff it has the larger
       h -- position was deliberately set to say the opposite. */
    const int ok = (pi_acted == hi_larger) && (pi_acted == pj_acted);
    message(
        "  h-primary designation hi_larger=%d (position reversed): "
        "pi_acted=%d pj_acted=%d -- %s",
        hi_larger, pi_acted, pj_acted, ok ? "PASS" : "FAIL");
    all_ok = all_ok && ok;
  }
  return all_ok;
}

/**
 * @brief Check the local mixed-band tie-break arm (both_updatable_here=0,
 * local_first=local_second=1): exactly one of the two argument orders must
 * act (apply a nonzero, exactly-conservative flux), the other must skip
 * (leave both particles' flux untouched); a tied position must skip both
 * ways.
 *
 * @return 1 on pass, 0 on failure (message()s the details either way).
 */
static int check_local_tiebreak(const struct cosmology *cosmo,
                                const struct chemistry_global_data *chem_data) {
  const float time_base = 1e-5;
  const integertime_t ti_current = 1;
  int all_ok = 1;

  const float test_order_signs[3] = {1.f, -1.f, 0.f};
  for (int t = 0; t < 3; t++) {
    struct part pi, pj;
    build_tiebreak_pair(&pi, &pj, test_order_signs[t]);
    const float dx[3] = {0.3f, 0.f, 0.f};
    const float r2 = dx[0] * dx[0];

    reset_flux_accumulators(&pi);
    reset_flux_accumulators(&pj);
    runner_iact_chemistry_flux_exchange(
        r2, dx, pi.h, pj.h, &pi, &pj, /*both_updatable_here=*/0,
        /*local_first=*/1, /*local_second=*/1, /*a=*/1.f, /*H=*/0.f, time_base,
        ti_current, cosmo, /*with_cosmology=*/0, chem_data);
    const int fwd_acted =
        (pi.chemistry_data.flux.metal_mass[TEST_METAL] != 0.0) ||
        (pj.chemistry_data.flux.metal_mass[TEST_METAL] != 0.0);

    struct part pi2, pj2;
    build_tiebreak_pair(&pi2, &pj2, test_order_signs[t]);
    const float mdx[3] = {-dx[0], 0.f, 0.f};
    reset_flux_accumulators(&pi2);
    reset_flux_accumulators(&pj2);
    runner_iact_chemistry_flux_exchange(
        r2, mdx, pj2.h, pi2.h, &pj2, &pi2, /*both_updatable_here=*/0,
        /*local_first=*/1, /*local_second=*/1, /*a=*/1.f, /*H=*/0.f, time_base,
        ti_current, cosmo, /*with_cosmology=*/0, chem_data);
    const int bwd_acted =
        (pi2.chemistry_data.flux.metal_mass[TEST_METAL] != 0.0) ||
        (pj2.chemistry_data.flux.metal_mass[TEST_METAL] != 0.0);

    int ok;
    const char *expect;
    if (test_order_signs[t] == 0.0f) {
      ok = !fwd_acted && !bwd_acted;
      expect = "tie -> both skip";
    } else {
      ok = (fwd_acted != bwd_acted); /* exactly one of the two acts */
      expect = "exactly one side acts";
    }
    if (ok && fwd_acted) {
      const double Fi = pi.chemistry_data.flux.metal_mass[TEST_METAL];
      const double Fj = pj.chemistry_data.flux.metal_mass[TEST_METAL];
      if (fabs(Fi + Fj) > 1e-12 * max(fabs(Fi), fabs(Fj))) ok = 0;
    }
    message(
        "  local tie-break order_sign=% .0f: fwd_acted=%d bwd_acted=%d (%s) "
        "-- %s",
        test_order_signs[t], fwd_acted, bwd_acted, expect,
        ok ? "PASS" : "FAIL");
    all_ok = all_ok && ok;
  }
  return all_ok;
}

/**
 * @brief Check the cross-rank arm (both_updatable_here=0, local_first=1,
 * local_second=0): the canonical-order computation must always apply only
 * to the local (first-argument) slot, whichever side wins the position
 * comparison, and a tied position must skip.
 *
 * @return 1 on pass, 0 on failure.
 */
static int check_cross_rank(const struct cosmology *cosmo,
                            const struct chemistry_global_data *chem_data) {
  const float time_base = 1e-5;
  const integertime_t ti_current = 1;
  int all_ok = 1;

  /* order_sign > 0 => pi is NOT position-canonical-left (pj is); order_sign
     < 0 => pi IS canonical-left. Both must update only pi (the local,
     first-argument, slot) and never touch pj (the foreign copy). */
  const float test_order_signs[3] = {1.f, -1.f, 0.f};
  for (int t = 0; t < 3; t++) {
    struct part pi, pj;
    build_tiebreak_pair(&pi, &pj, test_order_signs[t]);
    const float dx[3] = {0.3f, 0.f, 0.f};
    const float r2 = dx[0] * dx[0];

    reset_flux_accumulators(&pi);
    reset_flux_accumulators(&pj);
    runner_iact_chemistry_flux_exchange(
        r2, dx, pi.h, pj.h, &pi, &pj, /*both_updatable_here=*/0,
        /*local_first=*/1, /*local_second=*/0, /*a=*/1.f, /*H=*/0.f, time_base,
        ti_current, cosmo, /*with_cosmology=*/0, chem_data);

    const double Fi = pi.chemistry_data.flux.metal_mass[TEST_METAL];
    const double Fj_foreign = pj.chemistry_data.flux.metal_mass[TEST_METAL];

    int ok;
    if (test_order_signs[t] == 0.0f) {
      ok = (Fi == 0.0) && (Fj_foreign == 0.0);
    } else {
      /* The foreign slot must NEVER be written; the local slot must get a
         non-zero, correctly-signed share. */
      ok = (Fj_foreign == 0.0) && (Fi != 0.0);
    }
    message(
        "  cross-rank order_sign=% .0f: Fi(local)=% .3e Fj(foreign)=% .3e -- "
        "%s",
        test_order_signs[t], Fi, Fj_foreign, ok ? "PASS" : "FAIL");
    all_ok = all_ok && ok;
  }
  return all_ok;
}

/* ------------------------------------------------------------------------ */
/* Main sweep                                                               */
/* ------------------------------------------------------------------------ */

int main(int argc, char *argv[]) {

  message("=== Chemistry flux-exchange dispatcher gate test ===");
#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  message("Build: HYPERBOLIC diffusion");
  const int have_hyp_flux_axis = 1;
#else
  message("Build: PARABOLIC diffusion");
  const int have_hyp_flux_axis = 0;
#endif

  const enum test_geometry geoms[2] = {GEOM_WELL_BEHAVED, GEOM_ILL_CONDITIONED};
  const float h_ratios[2] = {1.0f, 8.0f};
  const enum test_limiter_regime limiters[5] = {
      LIMITER_INERT, LIMITER_BITING, LIMITER_NEAR_NOISE_GATE,
      LIMITER_NEAR_CAPACITY_POS, LIMITER_NEAR_CAPACITY_NEG};
  const int hyp_flux_options[2] = {0, 1};
  const int n_hyp_flux_options = have_hyp_flux_axis ? 2 : 1;
  const int fct_options[2] = {0, 1};

  double worst_gate_residual = 0.0;
  char worst_gate_label[256] = "";

  for (int fct_i = 0; fct_i < 2; fct_i++) {
    const int fct_enabled = fct_options[fct_i];
    message("--- GATE sweep (both_updatable_here=1), FCT %s ---",
            fct_enabled ? "ENABLED (theta=0.3)" : "disabled (theta=1.0)");

    for (int g = 0; g < 2; g++) {
      for (int hr = 0; hr < 2; hr++) {
        for (int l = 0; l < 5; l++) {
          for (int hf = 0; hf < n_hyp_flux_options; hf++) {
            char label[256];
            snprintf(label, sizeof(label),
                     "geom=%s h_ratio=%.0f limiter=%s hyp_flux=%d",
                     geom_name(geoms[g]), h_ratios[hr],
                     limiter_name(limiters[l]), hyp_flux_options[hf]);
            struct residual_result r;
            run_config(label, geoms[g], h_ratios[hr], limiters[l],
                       hyp_flux_options[hf], fct_enabled,
                       /*mixed_activity=*/0, /*with_cosmology=*/0, &r, NULL);
            if (fabs(r.signed_residual) > fabs(worst_gate_residual)) {
              worst_gate_residual = r.signed_residual;
              snprintf(worst_gate_label, sizeof(worst_gate_label), "%s [%s]",
                       label, fct_enabled ? "FCT on" : "FCT off");
            }
          }
        }
      }
    }
  }

  message("=== mindt order-independence check (mixed activity) ===");
  double worst_mindt_residual = 0.0;
  int mindt_pass = 1;
  for (int g = 0; g < 2; g++) {
    for (int hr = 0; hr < 2; hr++) {
      char label[256];
      snprintf(label, sizeof(label), "geom=%s h_ratio=%.0f",
               geom_name(geoms[g]), h_ratios[hr]);
      double mindt_residual;
      run_config(label, geoms[g], h_ratios[hr], LIMITER_INERT,
                 /*hyp_flux_nonzero=*/0, /*fct_enabled=*/0,
                 /*mixed_activity=*/1, /*with_cosmology=*/0, NULL,
                 &mindt_residual);
      if (fabs(mindt_residual) > fabs(worst_mindt_residual))
        worst_mindt_residual = mindt_residual;
      if (fabs(mindt_residual) > 1e-12) mindt_pass = 0;
    }
  }
  message("Worst mindt residual: %.3e (%s)", worst_mindt_residual,
          mindt_pass ? "PASS (order-independent)" : "FAIL (order-dependent)");

  message("=== local mixed-band tie-break arm ===");
  struct chemistry_global_data chem_data_tb;
  build_chem_data(&chem_data_tb, 1e-15);
  struct cosmology cosmo_tb;
  cosmology_init_no_cosmo(&cosmo_tb);
  const int tiebreak_pass = check_local_tiebreak(&cosmo_tb, &chem_data_tb);
  message("Local tie-break arm: %s", tiebreak_pass ? "PASS" : "FAIL");

  message("=== h-primary designation (A1) ===");
  const int h_primary_pass =
      check_h_primary_designation(&cosmo_tb, &chem_data_tb);
  message("h-primary designation: %s", h_primary_pass ? "PASS" : "FAIL");

  message("=== cross-rank arm ===");
  const int cross_rank_pass = check_cross_rank(&cosmo_tb, &chem_data_tb);
  message("Cross-rank arm: %s", cross_rank_pass ? "PASS" : "FAIL");

  message("=== Summary ===");
  message("Worst GATE residual (both_updatable_here=1): %.3e  [%s]",
          worst_gate_residual, worst_gate_label);
  message("Worst mindt residual (mixed activity): %.3e", worst_mindt_residual);
  message("Local tie-break arm: %s", tiebreak_pass ? "PASS" : "FAIL");
  message("h-primary designation: %s", h_primary_pass ? "PASS" : "FAIL");
  message("Cross-rank arm: %s", cross_rank_pass ? "PASS" : "FAIL");

  const int overall_pass =
      mindt_pass && tiebreak_pass && h_primary_pass && cross_rank_pass;
  if (!overall_pass) error("testChemistryFluxAntisymmetry: FAILED");

  return 0;
}

#else

int main(int argc, char **argv) { return 0; }

#endif
