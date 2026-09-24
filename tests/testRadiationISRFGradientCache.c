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

/* The gradient loop's real dispatch, `runner_iact_[nonsym_]isrf_gradient`,
 * reads its M1 closure tensor from #feedback_isrf_moment_data's per-particle
 * cache (#radiation_cache_m1_closure_part), not from the from-scratch
 * #radiation_get_m1_closure_tensor_band every other ISRF test exercises. A
 * test that only ever calls the from-scratch path would keep passing even if
 * the cache silently wrote zeros. This drives the real dispatch on two
 * particles and checks it against the from-scratch tensor built on the same
 * (u, F, c_hyp) state. */
#if defined(FEEDBACK_GEAR)

#include "feedback/GEAR/radiation_propagation_iact.h"

#define NUM_SCENARIOS 3

struct scenario {
  const char *name;
  float u_i, u_j;
  float F_i[3], F_j[3];
  float c_hyp;
};

static const struct scenario scenarios[NUM_SCENARIOS] = {
    {"generic nonzero flux",
     3.0f,
     1.4f,
     {0.4f, 0.1f, -0.2f},
     {-0.3f, 0.25f, 0.05f},
     2.0f},
    /* i: |F|/(c_hyp*u) = 3.0/2.0, clamped at the limiter (f = 1).
     * j: |F|/(c_hyp*u) = 2.39/2.4 = 0.9958, just under it. */
    {"near/at the flux limiter",
     1.0f,
     1.2f,
     {3.0f, 0.f, 0.f},
     {0.f, 2.39f, 0.f},
     2.0f},
    {"zero flux (isotropic guard)",
     2.0f,
     1.5f,
     {0.f, 0.f, 0.f},
     {0.f, 0.f, 0.f},
     2.0f},
};

/* Comoving pair geometry, deliberately asymmetric in every input so an
 * accidentally swapped i/j term cannot hide. */
static const float dx[3] = {0.5f, -0.4f, 0.3f};
static const float hi = 0.9f, hj = 1.1f;
static const float rho_i = 0.7f, rho_j = 1.3f;
static const float mass_i = 1.0f, mass_j = 2.5f;

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
 */
static void check_close(const char *name, float expected, float value) {

  check_finite(name, value);
  const float scale = fmaxf(fabsf(expected), 1e-6f);
  if (!(fabsf(value - expected) <= 1e-5f * scale))
    error("%s: expected %.9e, got %.9e", name, expected, value);
}

/**
 * @brief Fail unless two floats are bitwise identical, on their raw bit
 * patterns (not `==`, which -ffast-math is free to fold differently). Both
 * are checked finite first.
 *
 * @param name Name of the quantity, for the failure message.
 * @param a First value.
 * @param b Second value.
 */
static void check_bits_equal(const char *name, float a, float b) {

  check_finite(name, a);
  check_finite(name, b);
  uint32_t bits_a, bits_b;
  memcpy(&bits_a, &a, sizeof(bits_a));
  memcpy(&bits_b, &b, sizeof(bits_b));
  if (bits_a != bits_b)
    error("%s: not bitwise identical (0x%08x vs 0x%08x, %.9e vs %.9e)", name,
          bits_a, bits_b, (double)a, (double)b);
}

/**
 * @brief Set both ISRF bands of a particle to the same `(u, F, c_hyp)`.
 *
 * @param p (return) The particle.
 * @param u The specific field.
 * @param F The specific flux.
 * @param rho_prev The comoving density snapshot.
 * @param c_hyp The hyperbolic propagation speed.
 * @param mass The particle mass.
 */
static void set_part(struct part *p, float u, const float F[3], float rho_prev,
                     float c_hyp, float mass) {

  bzero(p, sizeof(struct part));
  p->mass = mass;
  p->id = 1;
  p->time_bin = 1;
#ifdef SWIFT_DEBUG_CHECKS
  p->ti_drift = 8;
  p->ti_kick = 8;
#endif

  struct feedback_part_data *fd = &p->feedback_data;
  fd->rho_prev = rho_prev;
  fd->c_hyp = c_hyp;
  for (int b = 0; b < ISRF_MOMENT_COUNT; b++) {
    struct feedback_isrf_moment_data *band = &fd->isrf_moment[b];
    band->u = u;
    for (int k = 0; k < 3; k++) band->specific_flux[k] = F[k];
  }
}

/**
 * @brief The kernel-derivative terms `runner_iact_isrf_gradient` computes
 * internally from `dx`, `hi`, `hj`, replicated so the from-scratch reference
 * below uses the identical values.
 */
static void kernel_terms(float *r_inv, float *wi_dr, float *wj_dr) {

  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);
  *r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;
  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  *wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  *wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  (void)wi;
  (void)wj;
}

/**
 * @brief `grad(u)` from the from-scratch closure tensor
 * (#radiation_get_m1_closure_tensor_band), the reference the cached dispatch
 * must reproduce.
 *
 * @param sc The scenario.
 * @param grad_u_i (return) Particle i's gradient.
 * @param grad_u_j (return) Particle j's gradient.
 */
static void reference_grad(const struct scenario *sc, float grad_u_i[3],
                           float grad_u_j[3]) {

  float r_inv, wi_dr, wj_dr;
  kernel_terms(&r_inv, &wi_dr, &wj_dr);

  float D_i[3][3], D_j[3][3];
  radiation_get_m1_closure_tensor_band(sc->u_i, sc->F_i, sc->c_hyp, D_i);
  radiation_get_m1_closure_tensor_band(sc->u_j, sc->F_j, sc->c_hyp, D_j);

  for (int k = 0; k < 3; k++) {
    grad_u_i[k] = 0.f;
    grad_u_j[k] = 0.f;
  }
  radiation_gradient_accumulate_band(dx, r_inv, wi_dr, wj_dr, mass_i, mass_j,
                                     rho_i, rho_j, sc->u_i, sc->u_j, D_i, D_j,
                                     1.f, grad_u_i, grad_u_j);
}

int main(int argc, char *argv[]) {

  const float a = 1.f, H = 0.f;
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  for (int s = 0; s < NUM_SCENARIOS; s++) {
    const struct scenario *sc = &scenarios[s];

    struct part p_i0, p_j0;
    set_part(&p_i0, sc->u_i, sc->F_i, rho_i, sc->c_hyp, mass_i);
    set_part(&p_j0, sc->u_j, sc->F_j, rho_j, sc->c_hyp, mass_j);
    radiation_cache_m1_closure_part(&p_i0);
    radiation_cache_m1_closure_part(&p_j0);

    float ref_grad_u_i[3], ref_grad_u_j[3];
    reference_grad(sc, ref_grad_u_i, ref_grad_u_j);

    /* Symmetric dispatch, from the cached closure tensor on both sides. */
    struct part p_i = p_i0, p_j = p_j0;
    runner_iact_isrf_gradient(r2, dx, hi, hj, &p_i, &p_j, a, H);
    for (int b = 0; b < ISRF_MOMENT_COUNT; b++) {
      for (int k = 0; k < 3; k++) {
        check_close("symmetric grad_u_i", ref_grad_u_i[k],
                    p_i.feedback_data.isrf_moment[b].grad_u[k]);
        check_close("symmetric grad_u_j", ref_grad_u_j[k],
                    p_j.feedback_data.isrf_moment[b].grad_u[k]);
      }
    }
    /* grad_u has no snapshot field in any build (Gate 1a's within-run
     * identity leg cannot observe it there), so this C unit test is the
     * only place ISRF_MOMENT_LW and ISRF_MOMENT_LW_PHOTON's grad_u are
     * checked bitwise identical against EACH OTHER, not merely against the
     * same tolerance-bound reference: set_part() above gives every moment
     * the identical (u, F), matching Stage 1's direct-assignment
     * injection guarantee. */
    for (int k = 0; k < 3; k++) {
      check_bits_equal(
          "symmetric grad_u_i, LW vs LW_PHOTON",
          p_i.feedback_data.isrf_moment[ISRF_MOMENT_LW].grad_u[k],
          p_i.feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].grad_u[k]);
      check_bits_equal(
          "symmetric grad_u_j, LW vs LW_PHOTON",
          p_j.feedback_data.isrf_moment[ISRF_MOMENT_LW].grad_u[k],
          p_j.feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].grad_u[k]);
    }

    /* Non-symmetric dispatch: only i's accumulator is written. */
    struct part p_i2 = p_i0, p_j2 = p_j0;
    runner_iact_nonsym_isrf_gradient(r2, dx, hi, hj, &p_i2, &p_j2, a, H);
    for (int b = 0; b < ISRF_MOMENT_COUNT; b++) {
      for (int k = 0; k < 3; k++) {
        check_close("non-symmetric grad_u_i", ref_grad_u_i[k],
                    p_i2.feedback_data.isrf_moment[b].grad_u[k]);
        /* j is not the writable side of this dispatch: its accumulator must
         * stay at its zero seed. */
        check_close("non-symmetric grad_u_j (must stay zero)", 0.f,
                    p_j2.feedback_data.isrf_moment[b].grad_u[k]);
      }
    }
    for (int k = 0; k < 3; k++)
      check_bits_equal(
          "non-symmetric grad_u_i, LW vs LW_PHOTON",
          p_i2.feedback_data.isrf_moment[ISRF_MOMENT_LW].grad_u[k],
          p_i2.feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].grad_u[k]);

    message("cached gradient IACT matches the from-scratch closure: %s",
            sc->name);
  }

  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif /* FEEDBACK_GEAR */
