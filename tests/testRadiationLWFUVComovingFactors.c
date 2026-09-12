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

/* Stage 2 is off in a default build, and its midpoint reconstruction carries
 * the only `a`-factor of this file that is not the shared
 * `a_factor_comoving_to_physical`. Turned on here so that factor is covered
 * too. Safe in THIS translation unit only because every function it drives is
 * a header-side inline taking plain scalars and arrays: no `struct part` is
 * ever handed to the library, so the guarded field the macro adds to
 * `struct feedback_part_data` cannot cause a layout mismatch. Keep it that
 * way; the expansion-term test next door is the one that exercises real
 * particles. */
#ifndef RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION
#define RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION
#endif

/* Some standard headers. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
static const float phys_grad_prev_i[3] = {0.2f, -0.1f, 0.05f};
static const float phys_grad_prev_j[3] = {-0.15f, 0.3f, 0.1f};
static const float c_hyp = 2.0f;
/* Passed as the ungated trigger coefficients, with the floor held at 0:
 * this test measures the operators' scale-factor scaling, which the pair
 * contrast gate must not enter. */
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
                                       rho_i, rho_j, F_i, F_j, a_factor,
                                       &out->div_F_i, &out->div_F_j);

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
      dx, r, hi, hj, wi_dr, wj_dr, mass_i, mass_j, rho_i, rho_j, c_hyp, c_hyp,
      alpha_i, alpha_j, /*alpha_floor_i=*/0.f, /*alpha_floor_j=*/0.f,
      /*ngb_mean_abs_u_V_i=*/0.f, /*ngb_mean_abs_u_V_j=*/0.f, u_i, u_j,
      phys_grad_prev_i, phys_grad_prev_j, a_factor, a, &out->dissipation_u_i,
      &out->dissipation_u_j);

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

/**
 * @brief Read a float's raw bits.
 *
 * `isnan()`/`isfinite()` are useless here: the tests build with the same
 * `-ffast-math` as the library, under which the compiler may fold them to a
 * constant. The bit pattern is the only reliable check.
 *
 * @param x The value to inspect.
 * @return Its IEEE-754 binary32 representation.
 */
static uint32_t float_bits(float x) {
  uint32_t u;
  memcpy(&u, &x, sizeof(u));
  return u;
}

/**
 * @brief Fail unless a pair carrying no field at all contributes exactly
 * nothing, for every gate knee across the calibrated range and beyond it.
 *
 * This is the whole box's state on step 0, and it is the one input for
 * which the contrast ratio's denominator reduces to
 * #RADIATION_LW_FUV_DISSIPATION_U_V_ABSOLUTE_FLOOR alone. Under
 * `-ffast-math` the gate's own division by `q0` is reassociated into that
 * denominator, and flush-to-zero then turns a denormal product into an
 * exact zero, so a guard chosen too close to the underflow threshold makes
 * the ratio `0/0` and poisons every particle's field with NaN from the
 * first step. The check is on the exact bit pattern, both because a NaN
 * must be detected without `isnan()` and because the physically correct
 * answer here is a bitwise zero, not a small number.
 */
static void check_degenerate_zero_field_pair(void) {

  const float knees[] = {1e-4f, 0.01f, 0.2f, 0.3f, 0.5f, 1.f};
  const float no_gradient[3] = {0.f, 0.f, 0.f};
  const float saved_q0 = radiation_lw_fuv_dissipation_pair_gate_q0;

  for (int i = 0; i < (int)(sizeof(knees) / sizeof(knees[0])); i++) {
    radiation_lw_fuv_dissipation_pair_gate_q0 = knees[i];

    float acc_i = 0.f, acc_j = 0.f;
    radiation_dissipation_force_accumulate_band(
        phys_dx, 1.f, phys_h_i, phys_h_j, /*wi_dr=*/-1.f, /*wj_dr=*/-1.f,
        mass_i, mass_j, phys_rho_i, phys_rho_j, c_hyp, c_hyp,
        /*alpha_trigger_i=*/0.f, /*alpha_trigger_j=*/0.f,
        /*alpha_floor_i=*/0.5f, /*alpha_floor_j=*/0.5f,
        /*ngb_mean_abs_u_V_i=*/0.f, /*ngb_mean_abs_u_V_j=*/0.f, /*u_i=*/0.f,
        /*u_j=*/0.f, no_gradient, no_gradient, /*a_factor=*/1.f, /*a=*/1.f,
        &acc_i, &acc_j);

    if (float_bits(acc_i) != 0u || float_bits(acc_j) != 0u)
      error(
          "A pair with no field contributed at gate knee q0 = %g: "
          "dissipation_u_i bits 0x%08x, dissipation_u_j bits 0x%08x",
          knees[i], float_bits(acc_i), float_bits(acc_j));
  }

  radiation_lw_fuv_dissipation_pair_gate_q0 = saved_q0;

  /* Checked separately from the call above, because the `min(q_ij_raw, 1)`
     bound in the gate can hide a guard with no headroom by forcing the
     ratio to be materialised before the division by `q0`. This asserts the
     constant's own margin, whatever the surrounding expression folds to.
     This static assertion, not the dynamic loop above, is the confirmed
     regression guard against the FLT_MIN underflow bug: do not delete it
     as "redundant" with the loop. */
  const float smallest_supported_knee = knees[0];
  const float guarded_denominator =
      RADIATION_LW_FUV_DISSIPATION_U_V_ABSOLUTE_FLOOR * smallest_supported_knee;
  if ((float_bits(guarded_denominator) & 0x7f800000u) == 0u)
    error(
        "RADIATION_LW_FUV_DISSIPATION_U_V_ABSOLUTE_FLOOR (%g) leaves no "
        "underflow headroom: times a gate knee of %g it is denormal (bits "
        "0x%08x) and flush-to-zero makes the contrast ratio 0/0",
        (double)RADIATION_LW_FUV_DISSIPATION_U_V_ABSOLUTE_FLOOR,
        (double)smallest_supported_knee, float_bits(guarded_denominator));

  message("Degenerate zero-field pair contributes exactly zero at every knee");
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

  check_degenerate_zero_field_pair();

  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* FEEDBACK_GEAR */
