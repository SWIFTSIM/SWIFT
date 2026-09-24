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
#include <stdio.h>
#include <stdlib.h>

/* Local headers. */
#include "swift.h"

#if defined(FEEDBACK_GEAR)

#include "feedback/GEAR/radiation_propagation_iact.h"

/* One physical two-particle state, expressed below at several scale factors.
 * Every input is deliberately asymmetric, so an accidentally symmetric
 * expression cannot hide a missing factor. */
static const float phys_dx[3] = {0.5f, -0.4f, 0.3f};
static const float phys_h_i = 0.9f;
static const float phys_h_j = 1.1f;
static const float phys_rho_i = 0.7f;
static const float phys_rho_j = 1.3f;
static const float mass_i = 1.0f;
static const float mass_j = 2.5f;
static const float u_i = 3.0f;
static const float u_j = 1.4f;
static const float F_i[3] = {0.4f, 0.1f, -0.2f};
static const float F_j[3] = {-0.3f, 0.25f, 0.05f};
static const float c_hyp = 2.0f;
/* Passed as the trigger coefficients, with the floor held at 0: one
 * component is enough, this test measuring scale-factor scaling only. */
static const float alpha_i = 0.3f;
static const float alpha_j = 0.45f;
static const float dt = 0.05f;

static const float scale_factors[4] = {1.f, 0.5f, 0.25f, 0.1f};
#define NUM_SCALE_FACTORS 4

/* Everything one evaluation of the three operators produces, at one scale
 * factor. */
struct operator_outputs {
  float div_F_i, div_F_j;
  float grad_u_i[3], grad_u_j[3];
  float dissipation_u_i, dissipation_u_j;
  float reduced_flux_i;
};

/**
 * @brief Evaluate the three spatial operators on the physical state above,
 * expressed in the comoving frame of a given scale factor.
 *
 * @param a The scale factor.
 * @param a_factor The conversion factor to hand the operators. Normally
 * `1/a`; passing `1.f` instead reproduces the pre-fix arithmetic exactly,
 * which is what the discrimination check below needs.
 * @param out (return) The operator outputs.
 */
static void evaluate_operators(float a, float a_factor,
                               struct operator_outputs *out) {

  /* Comoving expression of the same physical state: lengths divide by `a`,
   * densities multiply by `a^dim`. */
  const float dx[3] = {phys_dx[0] / a, phys_dx[1] / a, phys_dx[2] / a};
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = 1.f / r;
  const float hi = phys_h_i / a;
  const float hj = phys_h_j / a;
  const float rho_i = phys_rho_i * pow_dimension(a);
  const float rho_j = phys_rho_j * pow_dimension(a);

  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;
  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);

  out->div_F_i = 0.f;
  out->div_F_j = 0.f;
  radiation_divergence_accumulate_band(dx, r_inv, wi_dr, wj_dr, mass_i, mass_j,
                                       rho_i, rho_j, F_i, F_j, c_hyp, c_hyp,
                                       a_factor, &out->div_F_i, &out->div_F_j);

  float D_i[3][3], D_j[3][3];
  radiation_get_m1_closure_tensor_band(u_i, F_i, c_hyp, D_i);
  radiation_get_m1_closure_tensor_band(u_j, F_j, c_hyp, D_j);

  for (int k = 0; k < 3; k++) {
    out->grad_u_i[k] = 0.f;
    out->grad_u_j[k] = 0.f;
  }
  radiation_gradient_accumulate_band(dx, r_inv, wi_dr, wj_dr, mass_i, mass_j,
                                     rho_i, rho_j, u_i, u_j, D_i, D_j, a_factor,
                                     out->grad_u_i, out->grad_u_j);

  out->dissipation_u_i = 0.f;
  out->dissipation_u_j = 0.f;
  radiation_dissipation_force_accumulate_band(
      wi_dr, wj_dr, mass_i, mass_j, rho_i, rho_j, c_hyp, c_hyp, alpha_i,
      alpha_j, /*alpha_floor_i=*/0.f, /*alpha_floor_j=*/0.f, u_i, u_j, a_factor,
      &out->dissipation_u_i, &out->dissipation_u_j);

  /* One flux relaxation step at zero opacity (decay = phi = 1), then the M1
   * reduced flux the closure branches on. This is the quantity the bug
   * actually biased: under-reading it by a factor `a` selects the diffusive
   * branch too early. */
  float F_new[3];
  for (int k = 0; k < 3; k++)
    F_new[k] = F_i[k] - c_hyp * c_hyp * dt * out->grad_u_i[k];
  const float F_new_norm =
      sqrtf(F_new[0] * F_new[0] + F_new[1] * F_new[1] + F_new[2] * F_new[2]);
  out->reduced_flux_i = F_new_norm / (c_hyp * u_i);
}

/**
 * @brief Fail unless two values agree to a relative tolerance.
 *
 * @param name Name of the quantity, for the failure message.
 * @param a The scale factor the second value was produced at.
 * @param reference The value at a = 1.
 * @param value The value at scale factor a.
 */
static void check_same(const char *name, float a, float reference,
                       float value) {

  const float scale = fmaxf(fabsf(reference), 1e-30f);
  const float rel = fabsf(value - reference) / scale;
  if (rel > 1e-5f)
    error("%s is not scale-factor independent: %.9e at a = 1 vs %.9e at a = %g",
          name, reference, value, a);
}

int main(int argc, char *argv[]) {

  struct operator_outputs out[NUM_SCALE_FACTORS];
  for (int i = 0; i < NUM_SCALE_FACTORS; i++) {
    const float a = scale_factors[i];
    evaluate_operators(a, 1.f / a, &out[i]);
  }

  for (int i = 1; i < NUM_SCALE_FACTORS; i++) {
    const float a = scale_factors[i];
    check_same("div(F)_i", a, out[0].div_F_i, out[i].div_F_i);
    check_same("div(F)_j", a, out[0].div_F_j, out[i].div_F_j);
    for (int k = 0; k < 3; k++) {
      check_same("grad(u)_i", a, out[0].grad_u_i[k], out[i].grad_u_i[k]);
      check_same("grad(u)_j", a, out[0].grad_u_j[k], out[i].grad_u_j[k]);
    }
    check_same("dissipation_u_i", a, out[0].dissipation_u_i,
               out[i].dissipation_u_i);
    check_same("dissipation_u_j", a, out[0].dissipation_u_j,
               out[i].dissipation_u_j);
    check_same("reduced flux f", a, out[0].reduced_flux_i,
               out[i].reduced_flux_i);
  }

  message("All operators scale-factor independent. f = %.8e",
          out[0].reduced_flux_i);

  /* Discrimination check, so a vacuous test cannot pass: the pre-fix
   * arithmetic (conversion factor held at 1) must show the `a`-bias this
   * whole exercise exists to remove. */
  struct operator_outputs unfixed;
  evaluate_operators(0.5f, 1.f, &unfixed);
  const float expected = out[0].div_F_i * 0.5f;
  if (fabsf(unfixed.div_F_i - expected) > 1e-5f * fabsf(expected))
    error(
        "Discrimination check broken: unconverted div(F)_i at a = 0.5 is "
        "%.9e, expected a times the physical value, %.9e",
        unfixed.div_F_i, expected);
  if (fabsf(unfixed.div_F_i - out[0].div_F_i) < 1e-3f * fabsf(out[0].div_F_i))
    error("Discrimination check is vacuous: the unconverted operator agrees");

  message("Discrimination check passed: unconverted div(F)_i = %.8e vs %.8e",
          unfixed.div_F_i, out[0].div_F_i);

  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* FEEDBACK_GEAR */
