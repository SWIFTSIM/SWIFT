/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026.
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

/* Some standard headers. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/* #radiation_dissipation_floor_relaxation_gate and its two neighbours,
 * #radiation_update_dissipation_alpha_band and
 * #radiation_dissipation_alpha_floor_band, on the degenerate inputs their own
 * doxygen documents a specific return value for: an exactly-zero residual
 * denominator, a disabled gate (eps_R <= 0), no relaxation timescale
 * (c_hyp <= 0 or kappa = H = 0), zero energy, and the alpha_prev boundaries
 * {0, alpha_max}. All three are `static INLINE` in radiation_isrf.c (not
 * declared in radiation_isrf.h), so they are exercised the same way
 * testRadiationISRFFluxUnderflow.c does: through
 * #radiation_end_gradient_propagation, reading back
 * #feedback_isrf_band_data.dissipation_alpha_trigger/dissipation_alpha_floor.
 */
#if defined(FEEDBACK_GEAR)

#include "feedback/GEAR/radiation_isrf.h"

/**
 * @brief Fail unless a float is finite, by its bit pattern: `isnan`/`isinf`
 * are unavailable under -ffast-math.
 *
 * @param name Name of the quantity, for the failure message.
 * @param x The value.
 */
static void check_finite(const char *name, float x) {

  uint32_t bits;
  memcpy(&bits, &x, sizeof(bits));
  if (((bits >> 23) & 0xFFu) == 0xFFu)
    error("%s is not finite (bits 0x%08x)", name, bits);
}

/**
 * @brief Fail unless two values agree to a relative tolerance, with an
 * absolute floor on the scale. Always checks finiteness first, so a NaN
 * never slips past by comparing false against the tolerance.
 *
 * @param name Name of the quantity, for the failure message.
 * @param expected The reference value.
 * @param value The value under test.
 * @param tol Relative tolerance.
 * @param floor Smallest scale the tolerance is applied to.
 */
static void check_close(const char *name, float expected, float value,
                        float tol, float floor) {

  check_finite(name, value);
  const float scale = fmaxf(fabsf(expected), floor);
  if (!(fabsf(value - expected) <= tol * scale))
    error("%s: expected %.9e, got %.9e", name, expected, value);
}

/**
 * @brief Assemble the minimal engine #radiation_end_gradient_propagation
 * reads, with every dissipation parameter exposed (the underflow test's
 * make_engine() fixes alpha_max/eps_1/eps_lambda; the gate's degenerate
 * branches need H and eps_R varied too).
 *
 * @param e (return) The engine.
 * @param cosmo The cosmology it points at.
 * @param fp The feedback properties it points at.
 * @param pc The physical constants it points at.
 * @param H #cosmology.H.
 * @param alpha_floor #feedback_props.ISRF_dissipation_alpha_floor.
 * @param eps_R #feedback_props.ISRF_dissipation_floor_relaxation_residual.
 */
static void make_engine(struct engine *e, struct cosmology *cosmo,
                        struct feedback_props *fp, struct phys_const *pc,
                        double H, float alpha_floor, float eps_R) {

  bzero(cosmo, sizeof(struct cosmology));
  bzero(fp, sizeof(struct feedback_props));
  bzero(pc, sizeof(struct phys_const));
  bzero(e, sizeof(struct engine));

  cosmo->a = 1.;
  cosmo->H = H;
  /* The gate's own `w = kappa + H/c` divides by the TRUE speed of light, not
   * c_hyp (radiation_isrf.c's radiation_dissipation_floor_relaxation_gate):
   * every H = 0 case below is unaffected by the value chosen here (0/c = 0
   * for any finite c > 0), but it must be nonzero to avoid a 0/0 when H is
   * also 0. */
  pc->const_speed_light_c = 1.e4;

  fp->ISRF_propagation = 1;
  fp->ISRF_dissipation_alpha_max = 1.f;
  fp->ISRF_dissipation_negativity_threshold = 0.1f;
  fp->ISRF_dissipation_alpha_floor = alpha_floor;
  fp->ISRF_dissipation_floor_h_over_lambda = 0.5f;
  fp->ISRF_dissipation_floor_relaxation_residual = eps_R;

  e->cosmology = cosmo;
  e->feedback_props = fp;
  e->physical_constants = pc;
}

/**
 * @brief Set one band's full pre-step state.
 *
 * @param p (in/out) The particle (already zeroed by the caller).
 * @param b Band index.
 * @param u This band's #feedback_isrf_band_data.u (and u_prev).
 * @param ngb_mean_abs_u_V This band's kernel-mean scratch.
 * @param kappa This band's absorption rate.
 * @param F This band's specific_flux, BEFORE this step's update.
 * @param grad_u This band's grad_u accumulator.
 * @param alpha_trigger This band's incoming dissipation_alpha_trigger
 * (the trigger's own memory of its previous output).
 */
static void set_band(struct part *p, int b, float u, float ngb_mean_abs_u_V,
                     float kappa, const float F[3], const float grad_u[3],
                     float alpha_trigger) {

  struct feedback_isrf_band_data *band = &p->feedback_data.isrf_band[b];
  band->u = u;
  band->u_prev = u;
  band->ngb_mean_abs_u_V = ngb_mean_abs_u_V;
  band->kappa = kappa;
  band->dissipation_alpha_trigger = alpha_trigger;
  for (int k = 0; k < 3; k++) {
    band->specific_flux[k] = F[k];
    band->grad_u[k] = grad_u[k];
  }
}

/**
 * @brief Zero a particle and set its particle-level (not per-band) state.
 *
 * @param p (return) The particle.
 * @param h Smoothing length.
 * @param c_hyp #feedback_part_data.c_hyp.
 * @param dt #feedback_part_data.dt_prev.
 */
static void init_part(struct part *p, float h, float c_hyp, float dt) {

  bzero(p, sizeof(struct part));
  p->h = h;
  p->feedback_data.c_hyp = c_hyp;
  p->feedback_data.dt_prev = dt;
  p->feedback_data.rho_prev = 1.f;
}

/* Common band-level defaults shared by most of the gate cases below: a
 * nonzero, non-quiescent energy (so "zero energy" is never conflated with
 * "zero flux"), and a trigger memory that does not itself perturb the
 * floor's own value (the trigger and floor components are independent
 * fields, see feedback_isrf_band_data's own doxygen). */
static const float u_default = 1.f;
static const float alpha_trigger_default = 0.f;

/**
 * @brief `eps_R <= 0` disables the gate outright (`s = 1` unconditionally),
 * which is checked BEFORE even the quiescent `F = grad_u = 0` case: a
 * disabled gate must not fall through to the "trivially at the fixed point"
 * branch and return 0 instead. Both bands get a different kappa, so a
 * baseline that had one accidentally take the same code path as the other
 * would not go unnoticed.
 */
static void test_eps_R_zero_disables_gate(void) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  const float alpha_floor = 0.5f;
  make_engine(&e, &cosmo, &fp, &pc, /*H=*/0., alpha_floor, /*eps_R=*/0.f);

  const float zero[3] = {0.f, 0.f, 0.f};
  const float kappa[ISRF_BAND_COUNT] = {1.f, 0.f};
  /* h = 1, eps_lambda = 0.5: x = kappa/0.5, floor_band = alpha_floor/(1+x^4).
   */
  const float expected[ISRF_BAND_COUNT] = {alpha_floor / 17.f, alpha_floor};

  struct part p;
  init_part(&p, /*h=*/1.f, /*c_hyp=*/2.f, /*dt=*/0.5f);
  for (int b = 0; b < ISRF_BAND_COUNT; b++)
    set_band(&p, b, u_default, /*ngb_mean_abs_u_V=*/1.f, kappa[b], zero, zero,
             alpha_trigger_default);

  radiation_end_gradient_propagation(&p, &e);

  for (int b = 0; b < ISRF_BAND_COUNT; b++)
    check_close("eps_R=0, quiescent F/grad_u", expected[b],
                p.feedback_data.isrf_band[b].dissipation_alpha_floor, 1e-5f,
                1e-3f);

  message("eps_R = 0: gate returns 1 even for a fully quiescent particle");
}

/**
 * @brief `c_hyp <= 0`: no relaxation timescale to settle F against, so the
 * gate returns 1 regardless of F/grad_u. c_hyp = 0 is not reachable in a
 * real run (radiation_snapshot_part_propagation derives it as
 * ISRF_c_hyp_margin*h/dt with ISRF_c_hyp_margin > 0 enforced at parse time,
 * h > 0 by the SPH density iteration's own invariant, dt floored at
 * FLT_MIN), but the gate's own doxygen documents this branch explicitly, so
 * it is checked here as a defensive contract, not a reachable-state claim.
 * Also confirms the surrounding update stays finite at c_hyp = 0: the M1
 * limiter's `c_M*u` term collapses to 0, so specific_flux is driven to
 * exactly zero (a limiter to "no propagation speed allows no flux", not a
 * division hazard).
 */
static void test_c_hyp_zero(void) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  const float alpha_floor = 0.5f;
  make_engine(&e, &cosmo, &fp, &pc, /*H=*/0., alpha_floor, /*eps_R=*/0.1f);

  const float F[3] = {3.f, 0.f, 0.f};
  const float grad_u[3] = {-5.f, 2.f, 1.f};

  struct part p;
  init_part(&p, /*h=*/1.f, /*c_hyp=*/0.f, /*dt=*/0.5f);
  for (int b = 0; b < ISRF_BAND_COUNT; b++)
    set_band(&p, b, u_default, /*ngb_mean_abs_u_V=*/1.f, /*kappa=*/1.f, F,
             grad_u, alpha_trigger_default);

  radiation_end_gradient_propagation(&p, &e);

  /* h = 1, kappa = 1, eps_lambda = 0.5: floor_band = alpha_floor/(1+2^4). */
  const float expected_floor = alpha_floor / 17.f;
  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    const struct feedback_isrf_band_data *band = &p.feedback_data.isrf_band[b];
    check_close("c_hyp=0 floor", expected_floor, band->dissipation_alpha_floor,
                1e-5f, 1e-3f);
    for (int k = 0; k < 3; k++) {
      check_finite("c_hyp=0 flux", band->specific_flux[k]);
      if (band->specific_flux[k] != 0.f)
        error("c_hyp=0 flux did not collapse to 0: got %.9e",
              band->specific_flux[k]);
    }
  }

  message("c_hyp = 0: gate returns 1, flux collapses to 0, both finite");
}

/**
 * @brief `w = kappa + H/c_hyp <= 0` (kappa = 0 and H = 0, the near-primordial
 * free-streaming limit Design B's own governing-equation derivation targets,
 * see this worktree's CLAUDE.md): the gate returns 1 unconditionally, before
 * ever forming F/grad_u's residual, so this is checked with an F/grad_u pair
 * that is NOT itself at R = 1 (unlike test_exact_R_equals_one below) to
 * confirm the w <= 0 short-circuit governs, not an accidental R = 1.
 */
static void test_w_zero_kappa_and_H_zero(void) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  const float alpha_floor = 0.5f;
  make_engine(&e, &cosmo, &fp, &pc, /*H=*/0., alpha_floor, /*eps_R=*/0.1f);

  const float F[3] = {3.f, -2.f, 0.5f};
  const float grad_u[3] = {1.f, 4.f, -1.f};

  struct part p;
  init_part(&p, /*h=*/1.f, /*c_hyp=*/2.f, /*dt=*/0.5f);
  for (int b = 0; b < ISRF_BAND_COUNT; b++)
    set_band(&p, b, u_default, /*ngb_mean_abs_u_V=*/1.f, /*kappa=*/0.f, F,
             grad_u, alpha_trigger_default);

  radiation_end_gradient_propagation(&p, &e);

  /* kappa = 0: x = 0, floor_band = alpha_floor exactly. */
  for (int b = 0; b < ISRF_BAND_COUNT; b++)
    check_close("kappa=H=0 floor", alpha_floor,
                p.feedback_data.isrf_band[b].dissipation_alpha_floor, 1e-5f,
                1e-3f);

  message("kappa = H = 0: gate returns 1 (no relaxation timescale)");
}

/**
 * @brief `kappa = 0` but `H > 0`: the Hubble term alone still supplies a
 * relaxation timescale (`w = H/c > 0`), so the gate falls through to the
 * real R computation rather than the `w <= 0` short-circuit above -- the
 * cosmological path that matters for a production (comoving) run at
 * near-primordial kappa. `H` is chosen equal to `c` (both otherwise
 * arbitrary here; the gate's own doxygen has why physically `c >> H`
 * always) purely so that `w = H/c = 1`, giving the same clean arithmetic as
 * before this function was fixed to divide by `c` instead of `c_hyp`: w =
 * H/c = 1e4/1e4 = 1. F = (1,0,0), grad_u = (-1,0,0): wx = 1*1 + 2*(-1) = -1,
 * num = 1; F_norm = G_norm = 1, den = 1*1 + 2*1 = 3; R = 1/3. eps_R = 0.5:
 * ratio2 = (R/eps_R)^2 = (2/3)^2 = 4/9.
 */
static void test_w_positive_via_hubble_term(void) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  const float alpha_floor = 0.5f;
  make_engine(&e, &cosmo, &fp, &pc, /*H=*/1.e4, alpha_floor, /*eps_R=*/0.5f);

  const float F[3] = {1.f, 0.f, 0.f};
  const float grad_u[3] = {-1.f, 0.f, 0.f};

  struct part p;
  init_part(&p, /*h=*/1.f, /*c_hyp=*/2.f, /*dt=*/0.5f);
  for (int b = 0; b < ISRF_BAND_COUNT; b++)
    set_band(&p, b, u_default, /*ngb_mean_abs_u_V=*/1.f, /*kappa=*/0.f, F,
             grad_u, alpha_trigger_default);

  radiation_end_gradient_propagation(&p, &e);

  const float expected = alpha_floor * (4.f / 9.f);
  for (int b = 0; b < ISRF_BAND_COUNT; b++)
    check_close("kappa=0, H>0 floor", expected,
                p.feedback_data.isrf_band[b].dissipation_alpha_floor, 1e-4f,
                1e-3f);

  message(
      "kappa = 0, H > 0: gate uses the Hubble term as the relaxation "
      "timescale (s = 4/9)");
}

/**
 * @brief The quiescent particle (`F = grad_u = 0` exactly, `w > 0`): the
 * only case the gate's own doxygen singles out as returning 0 rather than
 * falling through to the `R` formula -- avoiding an unfiltered 0/0 the
 * squared-ratio folding under -ffast-math would otherwise produce. Distinct
 * from test_eps_R_zero_disables_gate's identical F/grad_u: there eps_R <= 0
 * pre-empts this branch and returns 1 instead, so the two tests together
 * pin down the documented precedence between the two special cases.
 */
static void test_den_zero_quiescent(void) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  const float alpha_floor = 0.5f;
  make_engine(&e, &cosmo, &fp, &pc, /*H=*/0., alpha_floor, /*eps_R=*/0.1f);

  const float zero[3] = {0.f, 0.f, 0.f};

  struct part p;
  init_part(&p, /*h=*/1.f, /*c_hyp=*/2.f, /*dt=*/0.5f);
  for (int b = 0; b < ISRF_BAND_COUNT; b++)
    set_band(&p, b, u_default, /*ngb_mean_abs_u_V=*/1.f, /*kappa=*/1.f, zero,
             zero, alpha_trigger_default);

  radiation_end_gradient_propagation(&p, &e);

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    const struct feedback_isrf_band_data *band = &p.feedback_data.isrf_band[b];
    check_finite("quiescent floor", band->dissipation_alpha_floor);
    if (band->dissipation_alpha_floor != 0.f)
      error("quiescent gate did not return 0: floor %.9e",
            band->dissipation_alpha_floor);
  }

  message(
      "F = grad_u = 0, w > 0: gate returns 0 (exactly at the fixed "
      "point), not 0/0");
}

/**
 * @brief `R = 1` exactly whenever exactly one of `F`, `grad_u` is zero (the
 * gate's own doxygen guarantee): tested at both the maximum allowed eps_R
 * (1.0, GEARFeedback:ISRF_dissipation_floor_relaxation_residual's own parsed
 * upper bound) and a small one, since `s = min(1, (R/eps_R)^2)` saturates to
 * 1 either way once R = eps_R. w = kappa = 1 (H = 0, c_hyp = 2): zero flux,
 * grad_u = (3,0,0): wx = 2*3 = 6, num = 6; F_norm = 0, G_norm = 3, den =
 * 2*3 = 6; R = 1. Zero grad_u, F = (3,0,0): wx = 1*3 = 3, num = 3; F_norm =
 * 3, G_norm = 0, den = 1*3 = 3; R = 1 again.
 */
static void test_exact_R_equals_one(void) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  const float alpha_floor = 0.5f;
  const float eps_R_values[2] = {1.0f, 0.05f};
  const float zero[3] = {0.f, 0.f, 0.f};
  const float V[3] = {3.f, 0.f, 0.f};
  /* h = 1, kappa = 1, eps_lambda = 0.5: floor_band = alpha_floor/(1+2^4). */
  const float expected = alpha_floor / 17.f;

  for (int i = 0; i < 2; i++) {
    make_engine(&e, &cosmo, &fp, &pc, /*H=*/0., alpha_floor, eps_R_values[i]);

    struct part p_zero_flux, p_zero_grad;
    init_part(&p_zero_flux, /*h=*/1.f, /*c_hyp=*/2.f, /*dt=*/0.5f);
    init_part(&p_zero_grad, /*h=*/1.f, /*c_hyp=*/2.f, /*dt=*/0.5f);
    for (int b = 0; b < ISRF_BAND_COUNT; b++) {
      set_band(&p_zero_flux, b, u_default, /*ngb_mean_abs_u_V=*/1.f,
               /*kappa=*/1.f, zero, V, alpha_trigger_default);
      set_band(&p_zero_grad, b, u_default, /*ngb_mean_abs_u_V=*/1.f,
               /*kappa=*/1.f, V, zero, alpha_trigger_default);
    }

    radiation_end_gradient_propagation(&p_zero_flux, &e);
    radiation_end_gradient_propagation(&p_zero_grad, &e);

    for (int b = 0; b < ISRF_BAND_COUNT; b++) {
      check_close(
          "R=1, zero flux", expected,
          p_zero_flux.feedback_data.isrf_band[b].dissipation_alpha_floor, 1e-5f,
          1e-3f);
      check_close(
          "R=1, zero grad_u", expected,
          p_zero_grad.feedback_data.isrf_band[b].dissipation_alpha_floor, 1e-5f,
          1e-3f);
    }
  }

  message(
      "R = 1 exactly whenever exactly one of F, grad_u is zero, at "
      "eps_R in {1.0, 0.05}");
}

/**
 * @brief Exact Fickian cancellation (`R = 0`, `w*F + c_hyp*grad_u = 0`
 * componentwise): the discrete steady state the floor's cost formula
 * assumes, distinct from the quiescent `F = grad_u = 0` case above. By the
 * triangle inequality `num = |w*F + c_hyp*grad_u| <= w|F| + c_hyp|grad_u| =
 * den` always, so `R <= 1` and `s <= 1` unconditionally: `min(1, ratio2)` is
 * a defensive clamp against float roundoff pushing R fractionally above 1,
 * never a bound this formula needs to enforce on paper. w = kappa = 1
 * (H = 0, c_hyp = 2): F = (2,0,0), grad_u = (-1,0,0): wx = 1*2 + 2*(-1) = 0
 * exactly (power-of-two operands, exact under any FMA contraction), so
 * num = 0 and R = 0 regardless of den.
 */
static void test_exact_R_equals_zero(void) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  const float alpha_floor = 0.5f;
  make_engine(&e, &cosmo, &fp, &pc, /*H=*/0., alpha_floor, /*eps_R=*/0.1f);

  const float F[3] = {2.f, 0.f, 0.f};
  const float grad_u[3] = {-1.f, 0.f, 0.f};

  struct part p;
  init_part(&p, /*h=*/1.f, /*c_hyp=*/2.f, /*dt=*/0.5f);
  for (int b = 0; b < ISRF_BAND_COUNT; b++)
    set_band(&p, b, u_default, /*ngb_mean_abs_u_V=*/1.f, /*kappa=*/1.f, F,
             grad_u, alpha_trigger_default);

  radiation_end_gradient_propagation(&p, &e);

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    const struct feedback_isrf_band_data *band = &p.feedback_data.isrf_band[b];
    check_finite("R=0 floor", band->dissipation_alpha_floor);
    if (band->dissipation_alpha_floor != 0.f)
      error("R=0 gate did not return 0: floor %.9e",
            band->dissipation_alpha_floor);
  }

  message(
      "F, grad_u exactly Fickian-balanced: gate returns 0 (R = 0, not "
      "the quiescent branch)");
}

/**
 * @brief #radiation_update_dissipation_alpha_band's trigger at the
 * documented `alpha_prev` boundaries {0, alpha_max}, crossed with "zero
 * energy" (`u_V = 0`, so `eps = 0`, `alpha_aim = 0`) and "deep negativity"
 * (`u_V << 0` with `ngb_mean_abs_u_V = 0`, so `eps = 1 >= eps_1`,
 * `alpha_aim = alpha_max`). All four combinations resolve through the
 * `alpha_aim >= alpha_prev` branch except zero-energy/alpha_prev=alpha_max,
 * which decays: kappa = 0 so a_kappa = 0, decay = exp(-c_hyp*dt/
 * (RADIATION_ISRF_DISSIPATION_DECAY_LENGTH*h)) = exp(-2*0.5/(5*1)) =
 * exp(-0.2).
 */
static void test_alpha_trigger_boundaries(void) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  make_engine(&e, &cosmo, &fp, &pc, /*H=*/0., /*alpha_floor=*/0.f,
              /*eps_R=*/0.f);
  const float alpha_max = fp.ISRF_dissipation_alpha_max;
  const float zero[3] = {0.f, 0.f, 0.f};

  struct {
    float u_V;
    float ngb_mean_abs_u_V;
    float alpha_prev;
    float expected;
    const char *label;
  } cases[4] = {
      {0.f, 1.f, 0.f, 0.f, "zero energy, alpha_prev=0"},
      {0.f, 1.f, alpha_max, alpha_max * expf(-0.2f),
       "zero energy, alpha_prev=alpha_max"},
      {-1.f, 0.f, 0.f, alpha_max, "deep negativity, alpha_prev=0"},
      {-1.f, 0.f, alpha_max, alpha_max,
       "deep negativity, alpha_prev=alpha_max"},
  };

  for (int i = 0; i < 4; i++) {
    struct part p;
    init_part(&p, /*h=*/1.f, /*c_hyp=*/2.f, /*dt=*/0.5f);
    for (int b = 0; b < ISRF_BAND_COUNT; b++)
      set_band(&p, b, cases[i].u_V, cases[i].ngb_mean_abs_u_V, /*kappa=*/0.f,
               zero, zero, cases[i].alpha_prev);

    radiation_end_gradient_propagation(&p, &e);

    for (int b = 0; b < ISRF_BAND_COUNT; b++)
      check_close(cases[i].label, cases[i].expected,
                  p.feedback_data.isrf_band[b].dissipation_alpha_trigger, 1e-5f,
                  1e-3f);
  }

  message(
      "alpha_prev at {0, alpha_max}, crossed with zero energy and deep "
      "negativity: all four finite and match the closed form");
}

int main(int argc, char *argv[]) {

  test_eps_R_zero_disables_gate();
  test_c_hyp_zero();
  test_w_zero_kappa_and_H_zero();
  test_w_positive_via_hubble_term();
  test_den_zero_quiescent();
  test_exact_R_equals_one();
  test_exact_R_equals_zero();
  test_alpha_trigger_boundaries();

  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* FEEDBACK_GEAR */
