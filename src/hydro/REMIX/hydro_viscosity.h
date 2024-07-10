/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_PLANETARY_HYDRO_VISCOSITY_H
#define SWIFT_PLANETARY_HYDRO_VISCOSITY_H

/**
 * @file Planetary/hydro_viscosity.h
 * @brief Utilities for hydro viscosity calculations under various options.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Prepares extra viscosity parameters for a particle for the density
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_init_part_extra_viscosity(struct part *restrict p) {}

/**
 * @brief Extra density interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_density_extra_viscosity(struct part *restrict pi,
                                          struct part *restrict pj,
                                          const float dx[3], const float wi,
                                          const float wj, const float wi_dx,
                                          const float wj_dx) {}

/**
 * @brief Extra density interaction between two particles (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_density_extra_viscosity(struct part *restrict pi,
                                                 const struct part *restrict pj,
                                                 const float dx[3],
                                                 const float wi,
                                                 const float wi_dx) {}

/**
 * @brief Finishes extra viscosity parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_density_extra_viscosity(struct part *restrict p) {}

/**
 * @brief Prepares extra viscosity parameters for a particle for the gradient
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_extra_viscosity(struct part *restrict p) {

  int i, j;
  for (i = 0; i < 3; ++i) {
    p->du_norm_kernel[i] = 0.f;
    p->drho_norm_kernel[i] = 0.f;
    p->dh_norm_kernel[i] = 0.f;
    for (j = 0; j < 3; ++j) {
      p->dv_norm_kernel[i][j] = 0.f;
    }
  }
}

/**
 * @brief Extra gradient interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_gradient_extra_viscosity(struct part *restrict pi,
                                           struct part *restrict pj,
                                           const float dx[3], const float wi,
                                           const float wj, const float wi_dx,
                                           const float wj_dx) {

  int i, j;
  float volume_i = pi->mass / pi->rho_evolved;
  float volume_j = pj->mass / pj->rho_evolved;

  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float hi = pi->h;
  const float hi_inv = 1.0f / hi;                        /* 1/h */
  const float hi_inv_dim = pow_dimension(hi_inv);        /* 1/h^d */
  const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */

  const float hj = pj->h;
  const float hj_inv = 1.0f / hj;                        /* 1/h */
  const float hj_inv_dim = pow_dimension(hj_inv);        /* 1/h^d */
  const float hj_inv_dim_plus_one = hj_inv_dim * hj_inv; /* 1/h^(d+1) */

  float wi_dx_term[3], wj_dx_term[3];
  for (i = 0; i < 3; ++i) {
    wi_dx_term[i] = dx[i] * r_inv * wi_dx * hi_inv_dim_plus_one;
    wj_dx_term[i] = -dx[i] * r_inv * wj_dx * hj_inv_dim_plus_one;

    wi_dx_term[i] *= (1.f / pi->m0_no_mean_kernel);
    wj_dx_term[i] *= (1.f / pj->m0_no_mean_kernel);

    wi_dx_term[i] += -wi * hi_inv_dim * pi->grad_m0_no_mean_kernel[i] /
                     (pi->m0_no_mean_kernel * pi->m0_no_mean_kernel);
    wj_dx_term[i] += -wj * hj_inv_dim * pj->grad_m0_no_mean_kernel[i] /
                     (pj->m0_no_mean_kernel * pj->m0_no_mean_kernel);
  }

  /* Set velocity derivative elements */
  for (i = 0; i < 3; ++i) {

    pi->dh_norm_kernel[i] += (pj->h - pi->h) * wi_dx_term[i] * volume_j;
    pj->dh_norm_kernel[i] += (pi->h - pj->h) * wj_dx_term[i] * volume_i;

    if (pi->mat_id == pj->mat_id) {

      pi->du_norm_kernel[i] += (pj->u - pi->u) * wi_dx_term[i] * volume_j;
      pj->du_norm_kernel[i] += (pi->u - pj->u) * wj_dx_term[i] * volume_i;

      pi->drho_norm_kernel[i] +=
          (pj->rho_evolved - pi->rho_evolved) * wi_dx_term[i] * volume_j;
      pj->drho_norm_kernel[i] +=
          (pi->rho_evolved - pj->rho_evolved) * wj_dx_term[i] * volume_i;
    }
    for (j = 0; j < 3; ++j) {
      /* Gradients from eq 18 in Rosswog 2020 (without C multiplied) */

      pi->dv_norm_kernel[i][j] +=
          (pj->v[j] - pi->v[j]) * wi_dx_term[i] * volume_j;
      pj->dv_norm_kernel[i][j] +=
          (pi->v[j] - pj->v[j]) * wj_dx_term[i] * volume_i;
    }
  }
}

/**
 * @brief Extra gradient interaction between two particles (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_gradient_extra_viscosity(
    struct part *restrict pi, const struct part *restrict pj, const float dx[3],
    const float wi, const float wi_dx) {

  int i, j;

  float volume_j = pj->mass / pj->rho_evolved;

  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float hi = pi->h;
  const float hi_inv = 1.0f / hi;                        /* 1/h */
  const float hi_inv_dim = pow_dimension(hi_inv);        /* 1/h^d */
  const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */

  float wi_dx_term[3];
  for (i = 0; i < 3; ++i) {
    wi_dx_term[i] = dx[i] * r_inv * wi_dx * hi_inv_dim_plus_one;

    wi_dx_term[i] *= (1.f / pi->m0_no_mean_kernel);

    wi_dx_term[i] += -wi * hi_inv_dim * pi->grad_m0_no_mean_kernel[i] /
                     (pi->m0_no_mean_kernel * pi->m0_no_mean_kernel);
  }

  /* Set velocity derivative elements */
  for (i = 0; i < 3; ++i) {

    pi->dh_norm_kernel[i] += (pj->h - pi->h) * wi_dx_term[i] * volume_j;

    if (pi->mat_id == pj->mat_id) {
      pi->du_norm_kernel[i] += (pj->u - pi->u) * wi_dx_term[i] * volume_j;
      pi->drho_norm_kernel[i] +=
          (pj->rho_evolved - pi->rho_evolved) * wi_dx_term[i] * volume_j;
    }
    for (j = 0; j < 3; ++j) {
      pi->dv_norm_kernel[i][j] +=
          (pj->v[j] - pi->v[j]) * wi_dx_term[i] * volume_j;
    }
  }
}

/**
 * @brief Finishes extra viscosity parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_extra_viscosity(struct part *restrict p) {}

/**
 * @brief Returns particle viscous pressures
 */
__attribute__((always_inline)) INLINE static void hydro_set_Qi_Qj(
    float *Qi, float *Qj, float *visc_signal_velocity,
    float *difn_signal_velocity, const struct part *restrict pi,
    const struct part *restrict pj, const float dx[3], const float a,
    const float H) {

  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);

  const float hi_inv = 1.0f / pi->h;
  const float hj_inv = 1.0f / pj->h;

  /* Quadratically reconstructed velocities at the halfway point between
   * particles */
  float vtilde_i[3], vtilde_j[3];

  /* Some parameters for artificial visc. Taken from Rosswog 2020 */
  const float alpha = planetary_reconst_visc_alpha;
  const float beta = planetary_reconst_visc_beta;
  const float epsilon = planetary_reconst_visc_epsilon;

  /* Square of eta (eq 16 in Rosswog 2020) */
  float eta_i_2 = r2 * hi_inv * hi_inv;
  float eta_j_2 = r2 * hj_inv * hj_inv;

  /* If h=h_max don't do anything fancy. Things like using m/rho to calculate
   * the volume stops working */
  if (!pi->is_h_max && !pj->is_h_max) {

    /* eq 23 in Rosswog 2020 */
    float eta_ab = min(r * hi_inv, r * hj_inv);  // min(r * hi_inv, r * hj_inv);

    /* A numerators and denominators (eq 22 in Rosswog 2020) */
    float A_i_v = 0.f;
    float A_j_v = 0.f;

    /* Terms in square brackets in Rosswog 2020 eq 17 */
    float v_reconst_i[3] = {0};
    float v_reconst_j[3] = {0};

    /* eq 23 in Rosswog 2020 set to constant */
    const float eta_crit = max(pi->force.eta_crit, pj->force.eta_crit);

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        /* Get the A numerators and denominators (eq 22 in Rosswog 2020). dv
         * is from eq 18 */
        A_i_v += pi->dv_norm_kernel[i][j] * dx[i] * dx[j];
        A_j_v += pj->dv_norm_kernel[i][j] * dx[i] * dx[j];

        /* Terms in square brackets in Rosswog 2020 eq 17. Add in FIRST
         * derivative terms */
        v_reconst_i[j] -= 0.5 * pi->dv_norm_kernel[i][j] * dx[i];
        v_reconst_j[j] += 0.5 * pj->dv_norm_kernel[i][j] * dx[i];
      }
    }

    float phi_i_v, phi_j_v;

    if (A_i_v == 0.f && A_j_v == 0.f) {
      phi_i_v = 1.f;
      phi_j_v = 1.f;

    } else if ((A_i_v == 0.f && A_j_v != 0.f) ||
               (A_j_v == 0.f && A_i_v != 0.f) || (A_i_v == -A_j_v)) {
      phi_i_v = 0.f;
      phi_j_v = 0.f;
    } else {
      /* Slope limiter (eq 21 in Rosswog 2020) */
      phi_i_v = min(1.f, 4.f * A_i_v / A_j_v / (1.f + A_i_v / A_j_v) /
                             (1.f + A_i_v / A_j_v));
      phi_i_v = max(0.f, phi_i_v);

      phi_j_v = min(1.f, 4.f * A_j_v / A_i_v / (1.f + A_j_v / A_i_v) /
                             (1.f + A_j_v / A_i_v));
      phi_j_v = max(0.f, phi_j_v);
    }

    if (eta_ab < eta_crit) {
      phi_i_v *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25.f);
      phi_j_v *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25.f);
    }

    for (int i = 0; i < 3; ++i) {
      /* Assemble the reconstructed velocity (eq 17 in Rosswog 2020) */
      vtilde_i[i] =
          pi->v[i] + (1.f - pi->force.balsara) * phi_i_v * v_reconst_i[i];
      vtilde_j[i] =
          pj->v[i] + (1.f - pj->force.balsara) * phi_j_v * v_reconst_j[i];
    }
  } else {
    for (int i = 0; i < 3; ++i) {
      /* If h=h_max don't reconstruct velocity */
      vtilde_i[i] = 0.f;
      vtilde_j[i] = 0.f;
    }
  }

  /* Finally assemble eq 15 in Rosswog 2020 */
  float mu_i = min(0.f, ((vtilde_i[0] - vtilde_j[0]) * dx[0] +
                         (vtilde_i[1] - vtilde_j[1]) * dx[1] +
                         (vtilde_i[2] - vtilde_j[2]) * dx[2]) *
                            hi_inv / (eta_i_2 + epsilon * epsilon));
  float mu_j = min(0.f, ((vtilde_i[0] - vtilde_j[0]) * dx[0] +
                         (vtilde_i[1] - vtilde_j[1]) * dx[1] +
                         (vtilde_i[2] - vtilde_j[2]) * dx[2]) *
                            hj_inv / (eta_j_2 + epsilon * epsilon));

  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  float balsara_i = pi->force.balsara / 3.f + 2.f / 3.f;
  float balsara_j = pj->force.balsara / 3.f + 2.f / 3.f;

  /* Get viscous pressure terms (eq 14 in Rosswog 2020) */
  *Qi = balsara_i * 0.5f * pi->rho * (-alpha * ci * mu_i + beta * mu_i * mu_i);
  *Qj = balsara_j * 0.5f * pj->rho * (-alpha * cj * mu_j + beta * mu_j * mu_j);

  float different_form_beta = 2.f * beta / alpha;
  *visc_signal_velocity = ci + cj - different_form_beta * min(mu_i, mu_j);

  float mean_balsara = 0.5f * (pi->force.balsara + pj->force.balsara);

  *difn_signal_velocity =
      (0.95f * mean_balsara + 0.05f) *
      sqrtf((vtilde_i[0] - vtilde_j[0]) * (vtilde_i[0] - vtilde_j[0]) +
            (vtilde_i[1] - vtilde_j[1]) * (vtilde_i[1] - vtilde_j[1]) +
            (vtilde_i[2] - vtilde_j[2]) * (vtilde_i[2] - vtilde_j[2]));
}

__attribute__((always_inline)) INLINE static void hydro_set_u_rho_difn(
    float *utilde_i, float *utilde_j, float *rhotilde_i, float *rhotilde_j,
    const struct part *restrict pi, const struct part *restrict pj,
    const float dx[3], const float a, const float H) {

  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);

  const float hi_inv = 1.0f / pi->h;
  const float hj_inv = 1.0f / pj->h;

  /* If h=h_max don't do anything fancy. Things like using m/rho to calculate
   * the volume stops working */
  if (!pi->is_h_max && !pj->is_h_max) {

    /* eq 23 in Rosswog 2020 */
    float eta_ab = min(r * hi_inv, r * hj_inv);  // min(r * hi_inv, r * hj_inv);

    /* A numerators and denominators (eq 22 in Rosswog 2020) */
    float A_i_u = 0.f;
    float A_j_u = 0.f;
    float A_i_rho = 0.f;
    float A_j_rho = 0.f;

    /* Terms in square brackets in Rosswog 2020 eq 17 */
    float u_reconst_i = 0.f;
    float u_reconst_j = 0.f;
    float rho_reconst_i = 0.f;
    float rho_reconst_j = 0.f;

    /* eq 23 in Rosswog 2020 set to constant */
    const float eta_crit = max(pi->force.eta_crit, pj->force.eta_crit);

    for (int i = 0; i < 3; ++i) {
      /* Get the A numerators and denominators (eq 22 in Rosswog 2020). dv
       * is from eq 18 */
      A_i_u += pi->du_norm_kernel[i] * dx[i];
      A_j_u += pj->du_norm_kernel[i] * dx[i];
      A_i_rho += pi->drho_norm_kernel[i] * dx[i];
      A_j_rho += pj->drho_norm_kernel[i] * dx[i];

      /* Terms in square brackets in Rosswog 2020 eq 17. Add in FIRST
       * derivative terms */
      u_reconst_i -= 0.5 * pi->du_norm_kernel[i] * dx[i];
      u_reconst_j += 0.5 * pj->du_norm_kernel[i] * dx[i];
      rho_reconst_i -= 0.5 * pi->drho_norm_kernel[i] * dx[i];
      rho_reconst_j += 0.5 * pj->drho_norm_kernel[i] * dx[i];
    }

    float phi_i_u, phi_j_u, phi_i_rho, phi_j_rho;

    if (A_i_u == 0.f && A_j_u == 0.f) {
      phi_i_u = 1.f;
      phi_j_u = 1.f;

    } else if ((A_i_u == 0.f && A_j_u != 0.f) ||
               (A_j_u == 0.f && A_i_u != 0.f) || (A_i_u == -A_j_u)) {
      phi_i_u = 0.f;
      phi_j_u = 0.f;
    } else {
      /* Slope limiter (eq 21 in Rosswog 2020) */
      phi_i_u = min(1.f, 4.f * A_i_u / A_j_u / (1.f + A_i_u / A_j_u) /
                             (1.f + A_i_u / A_j_u));
      phi_i_u = max(0.f, phi_i_u);

      phi_j_u = min(1.f, 4.f * A_j_u / A_i_u / (1.f + A_j_u / A_i_u) /
                             (1.f + A_j_u / A_i_u));
      phi_j_u = max(0.f, phi_j_u);
    }

    if (A_i_rho == 0.f && A_j_rho == 0.f) {
      phi_i_rho = 1.f;
      phi_j_rho = 1.f;

    } else if ((A_i_rho == 0.f && A_j_rho != 0.f) ||
               (A_j_rho == 0.f && A_i_rho != 0.f) || (A_i_rho == -A_j_rho)) {
      phi_i_rho = 0.f;
      phi_j_rho = 0.f;
    } else {
      /* Slope limiter (eq 21 in Rosswog 2020) */
      phi_i_rho = min(1.f, 4.f * A_i_rho / A_j_rho / (1.f + A_i_rho / A_j_rho) /
                               (1.f + A_i_rho / A_j_rho));
      phi_i_rho = max(0.f, phi_i_rho);

      phi_j_rho = min(1.f, 4.f * A_j_rho / A_i_rho / (1.f + A_j_rho / A_i_rho) /
                               (1.f + A_j_rho / A_i_rho));
      phi_j_rho = max(0.f, phi_j_rho);
    }

    if (eta_ab < eta_crit) {
      phi_i_u *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25.f);
      phi_j_u *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25.f);
      phi_i_rho *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25.f);
      phi_j_rho *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25.f);
    }

    /* Assemble the reconstructed velocity (eq 17 in Rosswog 2020) */
    *utilde_i = pi->u + phi_i_u * u_reconst_i;
    *utilde_j = pj->u + phi_j_u * u_reconst_j;
    *rhotilde_i = pi->rho_evolved + phi_i_rho * rho_reconst_i;
    *rhotilde_j = pj->rho_evolved + phi_j_rho * rho_reconst_j;

  } else {
    for (int i = 0; i < 3; ++i) {
      /* If h=h_max don't reconstruct velocity */
      *utilde_i = 0.f;
      *utilde_j = 0.f;
      *rhotilde_i = 0.f;
      *rhotilde_j = 0.f;
    }
  }
}

#endif /* SWIFT_PLANETARY_HYDRO_VISCOSITY_H */
