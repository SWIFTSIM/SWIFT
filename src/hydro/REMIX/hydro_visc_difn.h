/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_REMIX_HYDRO_VISC_DIFN_H
#define SWIFT_REMIX_HYDRO_VISC_DIFN_H

/**
 * @file REMIX/hydro_visc_difn.h
 * @brief Utilities for REMIX artificial viscosity and diffusion calculations.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Prepares extra artificial viscosity and artificial diffusion
 * parameters for a particle for the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_extra_visc_difn(struct part *restrict p) {

  memset(p->du_norm_kernel, 0.f, 3 * sizeof(float));
  memset(p->drho_norm_kernel, 0.f, 3 * sizeof(float));
  memset(p->dh_norm_kernel, 0.f, 3 * sizeof(float));
  memset(p->dv_norm_kernel, 0.f, 3 * 3 * sizeof(float));
}

/**
 * @brief Extra artificial viscosity and artificial diffusion gradient
 * interaction between two particles
 *
 * @param pi First particle.
 * @param pj Second particle.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param wi The value of the unmodified kernel function W(r, hi) * hi^d.
 * @param wj The value of the unmodified kernel function W(r, hj) * hj^d.
 * @param wi_dx The norm of the gradient of wi: dW(r, hi)/dr * hi^(d+1).
 * @param wj_dx The norm of the gradient of wj: dW(r, hj)/dr * hj^(d+1).
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_gradient_extra_visc_difn(struct part *restrict pi,
                                           struct part *restrict pj,
                                           const float dx[3], const float wi,
                                           const float wj, const float wi_dx,
                                           const float wj_dx) {

  /* Get r and 1/r. */
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

  const float volume_i = pi->mass / pi->rho_evol;
  const float volume_j = pj->mass / pj->rho_evol;

  /* Gradients of normalised kernel (Sandnes+2025 Eqn. 30) */
  float wi_dx_term[3], wj_dx_term[3];
  for (int i = 0; i < 3; i++) {
    wi_dx_term[i] = dx[i] * r_inv * wi_dx * hi_inv_dim_plus_one;
    wj_dx_term[i] = -dx[i] * r_inv * wj_dx * hj_inv_dim_plus_one;

    wi_dx_term[i] *= (1.f / pi->m0);
    wj_dx_term[i] *= (1.f / pj->m0);

    wi_dx_term[i] += -wi * hi_inv_dim * pi->grad_m0[i] / (pi->m0 * pi->m0);
    wj_dx_term[i] += -wj * hj_inv_dim * pj->grad_m0[i] / (pj->m0 * pj->m0);
  }

  /* Gradient estimates of h (Sandnes+2025 Eqn. 31), v (Eqn. 35), u (Eqn. 46),
   * rho (Eqn. 47) using normalised kernel */
  for (int i = 0; i < 3; i++) {
    pi->dh_norm_kernel[i] += (pj->h - pi->h) * wi_dx_term[i] * volume_j;
    pj->dh_norm_kernel[i] += (pi->h - pj->h) * wj_dx_term[i] * volume_i;

    /* Contributions only from same-material particles for u and rho diffusion
     */
    if (pi->mat_id == pj->mat_id) {
      pi->du_norm_kernel[i] += (pj->u - pi->u) * wi_dx_term[i] * volume_j;
      pj->du_norm_kernel[i] += (pi->u - pj->u) * wj_dx_term[i] * volume_i;

      pi->drho_norm_kernel[i] +=
          (pj->rho_evol - pi->rho_evol) * wi_dx_term[i] * volume_j;
      pj->drho_norm_kernel[i] +=
          (pi->rho_evol - pj->rho_evol) * wj_dx_term[i] * volume_i;
    }

    for (int j = 0; j < 3; j++) {
      pi->dv_norm_kernel[i][j] +=
          (pj->v[j] - pi->v[j]) * wi_dx_term[i] * volume_j;
      pj->dv_norm_kernel[i][j] +=
          (pi->v[j] - pj->v[j]) * wj_dx_term[i] * volume_i;
    }
  }
}

/**
 * @brief Extra artificial viscosity and artificial diffusion gradient
 * interaction between two particles (non-symmetric)
 *
 * @param pi First particle.
 * @param pj Second particle.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param wi The value of the unmodified kernel function W(r, hi) * hi^d.
 * @param wi_dx The norm of the gradient of wi: dW(r, hi)/dr * hi^(d+1).
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_gradient_extra_visc_difn(
    struct part *restrict pi, const struct part *restrict pj, const float dx[3],
    const float wi, const float wi_dx) {

  /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float hi = pi->h;
  const float hi_inv = 1.0f / hi;                        /* 1/h */
  const float hi_inv_dim = pow_dimension(hi_inv);        /* 1/h^d */
  const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */

  const float volume_j = pj->mass / pj->rho_evol;

  /* Gradients of normalised kernel (Sandnes+2025 Eqn. 30) */
  float wi_dx_term[3];
  for (int i = 0; i < 3; i++) {
    wi_dx_term[i] = dx[i] * r_inv * wi_dx * hi_inv_dim_plus_one;
    wi_dx_term[i] *= (1.f / pi->m0);
    wi_dx_term[i] += -wi * hi_inv_dim * pi->grad_m0[i] / (pi->m0 * pi->m0);
  }

  /* Gradient estimates of h (Sandnes+2025 Eqn. 31), v (Eqn. 35), u (Eqn. 46),
   * rho (Eqn. 47) using normalised kernel */
  for (int i = 0; i < 3; i++) {
    pi->dh_norm_kernel[i] += (pj->h - pi->h) * wi_dx_term[i] * volume_j;

    /* Contributions only from same-material particles for u and rho diffusion
     */
    if (pi->mat_id == pj->mat_id) {
      pi->du_norm_kernel[i] += (pj->u - pi->u) * wi_dx_term[i] * volume_j;
      pi->drho_norm_kernel[i] +=
          (pj->rho_evol - pi->rho_evol) * wi_dx_term[i] * volume_j;
    }

    for (int j = 0; j < 3; j++) {
      pi->dv_norm_kernel[i][j] +=
          (pj->v[j] - pi->v[j]) * wi_dx_term[i] * volume_j;
    }
  }
}

/**
 * @brief Finishes extra artificial viscosity and artificial diffusion parts of
 * the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_extra_visc_difn(struct part *restrict p) {}

/**
 * @brief Returns particle viscous pressures
 *
 * @param Qi (return) Viscous pressure for first particle.
 * @param Qj (return) Viscous pressure for second particle.
 * @param visc_signal_velocity (return) Signal velocity of artificial viscosity.
 * @param difn_signal_velocity (return) Signal velocity of artificial diffusion.
 * @param pi First particle.
 * @param pj Second particle.
 * @param dx Comoving vector separating both particles (pi - pj).
 */
__attribute__((always_inline)) INLINE static void hydro_set_Qi_Qj(
    float *Qi, float *Qj, float *visc_signal_velocity,
    float *difn_signal_velocity, const struct part *restrict pi,
    const struct part *restrict pj, const float dx[3]) {

  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);

  const float hi_inv = 1.0f / pi->h;
  const float hj_inv = 1.0f / pj->h;

  /* Reconstructed velocities at the halfway point between particles
   * (Sandnes+2025 Eqn. 36) */
  float vtilde_i[3], vtilde_j[3];

  /* Viscosity parameters. These are set in hydro_parameters.h  */
  const float alpha = const_remix_visc_alpha;
  const float beta = const_remix_visc_beta;
  const float epsilon = const_remix_visc_epsilon;
  const float a_visc = const_remix_visc_a;
  const float b_visc = const_remix_visc_b;
  const float eta_crit = 0.5f * (pi->force.eta_crit + pj->force.eta_crit);
  const float slope_limiter_exp_denom = const_remix_slope_limiter_exp_denom;

  if ((pi->is_h_max) || (pj->is_h_max)) {
    /* Don't reconstruct velocity if either particle has h=h_max */
    vtilde_i[0] = 0.f;
    vtilde_i[1] = 0.f;
    vtilde_i[2] = 0.f;

    vtilde_j[0] = 0.f;
    vtilde_j[1] = 0.f;
    vtilde_j[2] = 0.f;

  } else {
    /* A numerators and denominators (Sandnes+2025 Eqn. 38) */
    float A_i_v = 0.f;
    float A_j_v = 0.f;

    /* 1/2 * (r_j - r_i) * dv/dr in second term of Sandnes+2025 Eqn. 38 */
    float v_reconst_i[3] = {0};
    float v_reconst_j[3] = {0};

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        /* Get the A numerators and denominators (Sandnes+2025 Eqn. 38).
         * dv_norm_kernel is from Eqn. 35 */
        A_i_v += pi->dv_norm_kernel[i][j] * dx[i] * dx[j];
        A_j_v += pj->dv_norm_kernel[i][j] * dx[i] * dx[j];

        v_reconst_i[j] -= 0.5 * pi->dv_norm_kernel[i][j] * dx[i];
        v_reconst_j[j] += 0.5 * pj->dv_norm_kernel[i][j] * dx[i];
      }
    }

    /* Slope limiter (Sandnes+2025 Eqn. 37) special cases */
    float phi_i_v, phi_j_v;
    if ((A_i_v == 0.f) && (A_j_v == 0.f)) {
      phi_i_v = 1.f;
      phi_j_v = 1.f;

    } else if ((A_i_v == 0.f && A_j_v != 0.f) ||
               (A_j_v == 0.f && A_i_v != 0.f) || (A_i_v == -A_j_v)) {
      phi_i_v = 0.f;
      phi_j_v = 0.f;
    } else {
      /* Slope limiter (Sandnes+2025 Eqn. 37) */
      phi_i_v = min(1.f, 4.f * A_i_v / A_j_v / (1.f + A_i_v / A_j_v) /
                             (1.f + A_i_v / A_j_v));
      phi_i_v = max(0.f, phi_i_v);

      phi_j_v = min(1.f, 4.f * A_j_v / A_i_v / (1.f + A_j_v / A_i_v) /
                             (1.f + A_j_v / A_i_v));
      phi_j_v = max(0.f, phi_j_v);
    }

    /* exp in slope limiter (middle case in Sandnes+2025 Eqn. 37) */
    const float eta_ab = min(r * hi_inv, r * hj_inv);
    if (eta_ab < eta_crit) {
      phi_i_v *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) /
                      slope_limiter_exp_denom);
      phi_j_v *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) /
                      slope_limiter_exp_denom);
    }

    for (int i = 0; i < 3; i++) {
      /* Assemble the reconstructed velocity (Sandnes+2025 Eqn. 36) */
      vtilde_i[i] =
          pi->v[i] + (1.f - pi->force.balsara) * phi_i_v * v_reconst_i[i];
      vtilde_j[i] =
          pj->v[i] + (1.f - pj->force.balsara) * phi_j_v * v_reconst_j[i];
    }
  }

  /* Assemble Sandnes+2025 Eqn. 40 */
  const float mu_i =
      min(0.f, ((vtilde_i[0] - vtilde_j[0]) * dx[0] +
                (vtilde_i[1] - vtilde_j[1]) * dx[1] +
                (vtilde_i[2] - vtilde_j[2]) * dx[2]) *
                   hi_inv / (r2 * hi_inv * hi_inv + epsilon * epsilon));
  const float mu_j =
      min(0.f, ((vtilde_i[0] - vtilde_j[0]) * dx[0] +
                (vtilde_i[1] - vtilde_j[1]) * dx[1] +
                (vtilde_i[2] - vtilde_j[2]) * dx[2]) *
                   hj_inv / (r2 * hj_inv * hj_inv + epsilon * epsilon));

  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Finally assemble viscous pressure terms (Sandnes+2025 41) */
  *Qi = (a_visc + b_visc * pi->force.balsara) * 0.5f * pi->rho *
        (-alpha * ci * mu_i + beta * mu_i * mu_i);
  *Qj = (a_visc + b_visc * pj->force.balsara) * 0.5f * pj->rho *
        (-alpha * cj * mu_j + beta * mu_j * mu_j);

  /* Account for alpha being outside brackets in timestep code */
  const float viscosity_parameter_factor = (alpha == 0.f) ? 0.f : beta / alpha;
  *visc_signal_velocity =
      ci + cj - 2.f * viscosity_parameter_factor * min(mu_i, mu_j);

  /* Signal velocity used for the artificial diffusion (Sandnes+2025 Eqns. 42
   * and 43) */
  *difn_signal_velocity =
      sqrtf((vtilde_i[0] - vtilde_j[0]) * (vtilde_i[0] - vtilde_j[0]) +
            (vtilde_i[1] - vtilde_j[1]) * (vtilde_i[1] - vtilde_j[1]) +
            (vtilde_i[2] - vtilde_j[2]) * (vtilde_i[2] - vtilde_j[2]));
}

/**
 * @brief Returns midpoint reconstructions of internal energies and densities
 *
 * @param utilde_i (return) u reconstructed to midpoint from first particle.
 * @param utilde_j (return) u reconstructed to midpoint from second particle.
 * @param rhotilde_i (return) rho reconstructed to midpoint from first particle.
 * @param rhotilde_j (return) rho reconstructed to midpoint from second
 * particle.
 * @param pi First particle.
 * @param pj Second particle.
 * @param dx Comoving vector separating both particles (pi - pj).
 */
__attribute__((always_inline)) INLINE static void hydro_set_u_rho_difn(
    float *utilde_i, float *utilde_j, float *rhotilde_i, float *rhotilde_j,
    const struct part *restrict pi, const struct part *restrict pj,
    const float dx[3]) {

  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);

  const float hi_inv = 1.0f / pi->h;
  const float hj_inv = 1.0f / pj->h;

  const float eta_crit = 0.5f * (pi->force.eta_crit + pj->force.eta_crit);
  const float slope_limiter_exp_denom = const_remix_slope_limiter_exp_denom;

  if ((pi->is_h_max) || (pj->is_h_max)) {
    /* Don't reconstruct internal energy of density if either particle has
     * h=h_max */
    *utilde_i = 0.f;
    *utilde_j = 0.f;
    *rhotilde_i = 0.f;
    *rhotilde_j = 0.f;

    return;
  }

  /* A numerators and denominators (Sandnes+2025 Eqns. 48 and 49) */
  float A_i_u = 0.f;
  float A_j_u = 0.f;
  float A_i_rho = 0.f;
  float A_j_rho = 0.f;

  /* 1/2 * (r_j - r_i) * du/dr in second term of Sandnes+2025 Eqn. 44 */
  float u_reconst_i = 0.f;
  float u_reconst_j = 0.f;

  /* 1/2 * (r_j - r_i) * drho/dr in second term of Sandnes+2025 Eqn. 45 */
  float rho_reconst_i = 0.f;
  float rho_reconst_j = 0.f;

  for (int i = 0; i < 3; i++) {
    /* Get the A numerators and denominators (Sandnes+2025 Eqns. 48 and 49).
     * du_norm_kernel is from Eqn. 46 and drho_norm_kernel is from Eqn. 47 */
    A_i_u += pi->du_norm_kernel[i] * dx[i];
    A_j_u += pj->du_norm_kernel[i] * dx[i];
    A_i_rho += pi->drho_norm_kernel[i] * dx[i];
    A_j_rho += pj->drho_norm_kernel[i] * dx[i];

    u_reconst_i -= 0.5 * pi->du_norm_kernel[i] * dx[i];
    u_reconst_j += 0.5 * pj->du_norm_kernel[i] * dx[i];
    rho_reconst_i -= 0.5 * pi->drho_norm_kernel[i] * dx[i];
    rho_reconst_j += 0.5 * pj->drho_norm_kernel[i] * dx[i];
  }

  float phi_i_u, phi_j_u, phi_i_rho, phi_j_rho;
  /* Slope limiter (Sandnes+2025 Eqn. 37) special cases */
  if ((A_i_u == 0.f) && (A_j_u == 0.f)) {
    phi_i_u = 1.f;
    phi_j_u = 1.f;

  } else if ((A_i_u == 0.f && A_j_u != 0.f) || (A_j_u == 0.f && A_i_u != 0.f) ||
             (A_i_u == -A_j_u)) {
    phi_i_u = 0.f;
    phi_j_u = 0.f;
  } else {
    /* Slope limiter (Sandnes+2025 Eqn. 37) */
    phi_i_u = min(1.f, 4.f * A_i_u / A_j_u / (1.f + A_i_u / A_j_u) /
                           (1.f + A_i_u / A_j_u));
    phi_i_u = max(0.f, phi_i_u);

    phi_j_u = min(1.f, 4.f * A_j_u / A_i_u / (1.f + A_j_u / A_i_u) /
                           (1.f + A_j_u / A_i_u));
    phi_j_u = max(0.f, phi_j_u);
  }

  /* Slope limiter (Sandnes+2025 Eqn. 37) special cases */
  if ((A_i_rho == 0.f) && (A_j_rho == 0.f)) {
    phi_i_rho = 1.f;
    phi_j_rho = 1.f;

  } else if ((A_i_rho == 0.f && A_j_rho != 0.f) ||
             (A_j_rho == 0.f && A_i_rho != 0.f) || (A_i_rho == -A_j_rho)) {
    phi_i_rho = 0.f;
    phi_j_rho = 0.f;
  } else {
    /* Slope limiter (Sandnes+2025 Eqn. 37) */
    phi_i_rho = min(1.f, 4.f * A_i_rho / A_j_rho / (1.f + A_i_rho / A_j_rho) /
                             (1.f + A_i_rho / A_j_rho));
    phi_i_rho = max(0.f, phi_i_rho);

    phi_j_rho = min(1.f, 4.f * A_j_rho / A_i_rho / (1.f + A_j_rho / A_i_rho) /
                             (1.f + A_j_rho / A_i_rho));
    phi_j_rho = max(0.f, phi_j_rho);
  }

  /* exp in slope limiter (middle case in Sandnes+2025 Eqn. 37) */
  const float eta_ab = min(r * hi_inv, r * hj_inv);
  if (eta_ab < eta_crit) {
    phi_i_u *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) /
                    slope_limiter_exp_denom);
    phi_j_u *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) /
                    slope_limiter_exp_denom);
    phi_i_rho *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) /
                      slope_limiter_exp_denom);
    phi_j_rho *= expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) /
                      slope_limiter_exp_denom);
  }

  /* Assemble the reconstructed internal energy (Sandnes+2025 Eqn. 44) and
   * density (Sandnes+2025 Eqn. 45) */
  *utilde_i = pi->u + phi_i_u * u_reconst_i;
  *utilde_j = pj->u + phi_j_u * u_reconst_j;
  *rhotilde_i = pi->rho_evol + phi_i_rho * rho_reconst_i;
  *rhotilde_j = pj->rho_evol + phi_j_rho * rho_reconst_j;
}

#endif /* SWIFT_REMIX_HYDRO_VISC_DIFN_H */
