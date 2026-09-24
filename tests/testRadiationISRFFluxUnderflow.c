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

/* The M1 closure, flux limiter and floor relaxation-residual gate at flux
 * magnitudes whose float32 square underflows (|F| < sqrt(FLT_MIN) ~ 1.1e-19).
 * The closure and the limiter are homogeneous of degree zero in (u, F), and
 * the gate in (F, grad_u), so each tiny-scale result must match the same
 * state scaled up by a large factor, where no underflow is possible. */
#if defined(FEEDBACK_GEAR)

#include "feedback/GEAR/radiation_isrf.h"
#include "feedback/GEAR/radiation_propagation_iact.h"

/* Unit flux direction, deliberately not along an axis. */
static const float F_dir[3] = {0.8571429f, -0.4285714f, 0.2857143f};

/* Scale that lifts every tiny state above float32 underflow of its square. */
static const float lift = 1e20f;

static const float c_hyp = 2.f;

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
 * absolute floor on the scale.
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
 * @brief Assemble the minimal engine #radiation_end_gradient_propagation reads.
 *
 * @param e (return) The engine.
 * @param cosmo The cosmology it points at.
 * @param fp The feedback properties it points at.
 * @param pc The physical constants it points at.
 * @param alpha_floor #feedback_props.ISRF_dissipation_alpha_floor.
 * @param eps_R #feedback_props.ISRF_dissipation_floor_relaxation_residual.
 */
static void make_engine(struct engine *e, struct cosmology *cosmo,
                        struct feedback_props *fp, struct phys_const *pc,
                        float alpha_floor, float eps_R) {

  bzero(cosmo, sizeof(struct cosmology));
  bzero(fp, sizeof(struct feedback_props));
  bzero(pc, sizeof(struct phys_const));
  bzero(e, sizeof(struct engine));

  cosmo->a = 1.;
  cosmo->H = 0.;

  fp->ISRF_propagation = 1;
  fp->ISRF_extinction_path_in_kernel_radii = 2.0f;
  fp->ISRF_dissipation_alpha_max = 1.f;
  fp->ISRF_dissipation_negativity_threshold = 0.1f;
  fp->ISRF_dissipation_alpha_floor = alpha_floor;
  fp->ISRF_dissipation_floor_h_over_lambda = 0.5f;
  fp->ISRF_dissipation_floor_relaxation_residual = eps_R;

  pc->const_speed_light_c = 1.e4;

  e->cosmology = cosmo;
  e->feedback_props = fp;
  e->physical_constants = pc;
}

/**
 * @brief Set both bands of a particle to the same `(u, F, grad_u, kappa)`.
 *
 * @param p (return) The particle.
 * @param u The specific field.
 * @param F The specific flux.
 * @param grad_u The gradient accumulator.
 * @param kappa The absorption rate.
 */
static void set_part(struct part *p, float u, const float F[3],
                     const float grad_u[3], float kappa) {

  bzero(p, sizeof(struct part));
  p->h = 1.f;

  struct feedback_part_data *fd = &p->feedback_data;
  fd->dt_prev = 0.5f;
  fd->c_hyp = c_hyp;
  fd->rho_prev = 1.f;
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
    struct feedback_isrf_operator_data *op =
        &fd->isrf_operator[radiation_isrf_moment_to_operator[m]];
    op->kappa = kappa;
    moment->u = u;
    moment->u_prev = u;
    op->ngb_mean_abs_u_V = 1.f;
    for (int k = 0; k < 3; k++) {
      moment->specific_flux[k] = F[k];
      moment->grad_u[k] = grad_u[k];
    }
  }
}

/**
 * @brief Closure tensor at flux magnitude `F_mag` and reduced flux 0.9,
 * against the same state lifted above underflow.
 *
 * @param F_mag The flux magnitude, internal units.
 */
static void test_closure(float F_mag) {

  const float f_target = 0.9f;
  const float u = F_mag / (c_hyp * f_target);
  const float F[3] = {F_mag * F_dir[0], F_mag * F_dir[1], F_mag * F_dir[2]};
  const float F_ref[3] = {F[0] * lift, F[1] * lift, F[2] * lift};

  float D[3][3], D_ref[3][3];
  radiation_get_m1_closure_tensor_band(u, F, c_hyp, D);
  radiation_get_m1_closure_tensor_band(u * lift, F_ref, c_hyp, D_ref);

  /* At f = 0.9 the beam tensor is far from isotropic: D_xx - 1/3 is O(0.1). */
  if (fabsf(D_ref[0][0] - 1.f / 3.f) < 0.05f)
    error("reference closure is unexpectedly close to isotropic");

  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      check_close("closure D", D_ref[a][b], D[a][b], 1e-5f, 1e-3f);

  message("|F| = %.2e: closure matches the lifted state (D_xx = %.6f)", F_mag,
          D[0][0]);
}

/**
 * @brief M1 flux limiter at flux magnitude `F_mag` with `c_hyp*u = |F|/3.5`,
 * driven through #radiation_end_gradient_propagation with no decay and no
 * gradient, so the flux reaches the limiter unchanged.
 *
 * @param F_mag The flux magnitude, internal units.
 */
static void test_limiter(float F_mag) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  make_engine(&e, &cosmo, &fp, &pc, /*alpha_floor=*/0.f, /*eps_R=*/0.f);

  const float ratio = 1.f / 3.5f;
  const float u = ratio * F_mag / c_hyp;
  const float F[3] = {F_mag * F_dir[0], F_mag * F_dir[1], F_mag * F_dir[2]};
  const float zero[3] = {0.f, 0.f, 0.f};

  struct part p;
  set_part(&p, u, F, zero, /*kappa=*/0.f);
  radiation_end_gradient_propagation(&p, &e);

  for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
    for (int k = 0; k < 3; k++)
      check_close("limited flux", ratio * F[k],
                  p.feedback_data.isrf_moment[m].specific_flux[k], 1e-5f,
                  1e-3f * F_mag);

  message("|F| = %.2e: limiter clamps |F| to c_hyp*u", F_mag);
}

/**
 * @brief Floor relaxation-residual gate with a tiny flux and a zero gradient
 * (`R = 1`, so the floor must survive at full strength), and the same with
 * the roles swapped.
 *
 * @param F_mag The flux (or gradient) magnitude, internal units.
 */
static void test_gate(float F_mag) {

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  const float alpha_floor = 0.5f;
  make_engine(&e, &cosmo, &fp, &pc, alpha_floor, /*eps_R=*/0.1f);

  /* h = 1, kappa = 1, eps_lambda = 0.5: floor = alpha_floor/(1 + 2^4). */
  const float kappa = 1.f;
  const float expected = alpha_floor / 17.f;
  const float u = 1.f;
  const float V[3] = {F_mag * F_dir[0], F_mag * F_dir[1], F_mag * F_dir[2]};
  const float zero[3] = {0.f, 0.f, 0.f};

  struct part p;
  set_part(&p, u, V, zero, kappa);
  radiation_end_gradient_propagation(&p, &e);
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
    check_close(
        "floor, tiny flux", expected,
        p.feedback_data.isrf_operator[radiation_isrf_moment_to_operator[m]]
            .dissipation_alpha_floor,
        1e-5f, 1e-3f);

  set_part(&p, u, zero, V, kappa);
  radiation_end_gradient_propagation(&p, &e);
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++)
    check_close(
        "floor, tiny gradient", expected,
        p.feedback_data.isrf_operator[radiation_isrf_moment_to_operator[m]]
            .dissipation_alpha_floor,
        1e-5f, 1e-3f);

  message("|F| = %.2e: gate keeps the full floor", F_mag);
}

/**
 * @brief Exactly zero flux: isotropic closure, limiter no-op, and a quiescent
 * particle's gate returning 0.
 */
static void test_zero_flux(void) {

  const float zero[3] = {0.f, 0.f, 0.f};

  float D[3][3];
  radiation_get_m1_closure_tensor_band(1e-20f, zero, c_hyp, D);
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      check_close("zero-flux closure", (a == b) ? 1.f / 3.f : 0.f, D[a][b],
                  1e-6f, 1.f);

  struct engine e;
  struct cosmology cosmo;
  struct feedback_props fp;
  struct phys_const pc;
  make_engine(&e, &cosmo, &fp, &pc, /*alpha_floor=*/0.5f, /*eps_R=*/0.1f);

  struct part p;
  set_part(&p, /*u=*/1e-20f, zero, zero, /*kappa=*/1.f);
  radiation_end_gradient_propagation(&p, &e);
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    const struct feedback_isrf_moment_data *moment =
        &p.feedback_data.isrf_moment[m];
    const struct feedback_isrf_operator_data *op =
        &p.feedback_data.isrf_operator[radiation_isrf_moment_to_operator[m]];
    for (int k = 0; k < 3; k++) {
      check_finite("zero flux", moment->specific_flux[k]);
      if (moment->specific_flux[k] != 0.f)
        error("zero flux did not stay zero: %.9e", moment->specific_flux[k]);
    }
    check_finite("quiescent floor", op->dissipation_alpha_floor);
    if (op->dissipation_alpha_floor != 0.f)
      error("quiescent gate did not return 0: floor %.9e",
            op->dissipation_alpha_floor);
  }

  message("zero flux: isotropic closure, limiter no-op, quiescent gate 0");
}

int main(int argc, char *argv[]) {

  /* 1e-18 is above the float32 underflow of |F|^2; the others are below. */
  const float magnitudes[4] = {1e-18f, 1e-19f, 3e-20f, 1e-20f};

  for (int i = 0; i < 4; i++) {
    test_closure(magnitudes[i]);
    test_limiter(magnitudes[i]);
    test_gate(magnitudes[i]);
  }
  test_zero_flux();

  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* FEEDBACK_GEAR */
