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
#ifndef SWIFT_PLANETARY_HYDRO_DENSITY_ESTIMATE_H
#define SWIFT_PLANETARY_HYDRO_DENSITY_ESTIMATE_H

/**
 * @file Planetary/hydro_density_estimate.h
 * @brief Utilities for hydro density estimate under various options.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Prepares extra density estimate parameters for a particle for the
 * density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_init_part_extra_density_estimate(struct part *restrict p) {

#ifdef PLANETARY_IMBALANCE
  // ### Don't think we need to set these two to 0?
  p->P = 0.f;
  p->T = 0.f;

  p->sum_rij[0] = 0.f;
  p->sum_rij[1] = 0.f;
  p->sum_rij[2] = 0.f;
  p->I = 0.f;
  p->sum_wij = 0.f;
#elif PLANETARY_SMOOTHING_CORRECTION
  p->drho_dh = 0.f;
  p->grad_rho[0] = 0.f;
  p->grad_rho[1] = 0.f;
  p->grad_rho[2] = 0.f;
#endif /* PLANETARY_IMBALANCE */
}

/**
 * @brief Extra density interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_density_extra_density_estimate(
    struct part *restrict pi, struct part *restrict pj, const float dx[3],
    const float wi, const float wj, const float wi_dx, const float wj_dx) {

#ifdef PLANETARY_IMBALANCE
  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Add contribution to kernel averages */
  pi->sum_wij += wi * mj;
  pj->sum_wij += wj * mi;

  /* Add contribution r_ij * mj * sqrt(Wij) */
  if (pi->mat_id == pj->mat_id) {
    pi->sum_rij[0] += -dx[0] * wi * mj;
    pi->sum_rij[1] += -dx[1] * wi * mj;
    pi->sum_rij[2] += -dx[2] * wi * mj;

    pj->sum_rij[0] += dx[0] * wj * mi;
    pj->sum_rij[1] += dx[1] * wj * mi;
    pj->sum_rij[2] += dx[2] * wj * mi;
  }

  if (pi->mat_id != pj->mat_id) {
    pi->sum_rij[0] += dx[0] * wi * mj;
    pi->sum_rij[1] += dx[1] * wi * mj;
    pi->sum_rij[2] += dx[2] * wi * mj;

    pj->sum_rij[0] += -dx[0] * wj * mi;
    pj->sum_rij[1] += -dx[1] * wj * mi;
    pj->sum_rij[2] += -dx[2] * wj * mi;
  }
#elif PLANETARY_SMOOTHING_CORRECTION
  /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  const float ui = r / pi->h;
  const float uj = r / pj->h;

  pi->drho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pj->drho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);

  pi->grad_rho[0] += dx[0] * wi_dx * r_inv * mj;
  pi->grad_rho[1] += dx[1] * wi_dx * r_inv * mj;
  pi->grad_rho[2] += dx[2] * wi_dx * r_inv * mj;

  pj->grad_rho[0] += -dx[0] * wj_dx * r_inv * mi;
  pj->grad_rho[1] += -dx[1] * wj_dx * r_inv * mi;
  pj->grad_rho[2] += -dx[2] * wj_dx * r_inv * mi;
#endif /* PLANETARY_IMBALANCE */
}

/**
 * @brief Extra density interaction between two particles (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_density_extra_density_estimate(
    struct part *restrict pi, const struct part *restrict pj, const float dx[3],
    const float wi, const float wi_dx) {

#ifdef PLANETARY_IMBALANCE
  const float mj = pj->mass;

  /* Add contribution to kernel averages */
  pi->sum_wij += wi * mj;

  /* Add contribution r_ij*mj*sqrt(Wij) */
  if (pi->mat_id == pj->mat_id) {
    pi->sum_rij[0] += -dx[0] * wi * mj;
    pi->sum_rij[1] += -dx[1] * wi * mj;
    pi->sum_rij[2] += -dx[2] * wi * mj;
  }

  if (pi->mat_id != pj->mat_id) {
    pi->sum_rij[0] += dx[0] * wi * mj;
    pi->sum_rij[1] += dx[1] * wi * mj;
    pi->sum_rij[2] += dx[2] * wi * mj;
  }
#elif PLANETARY_SMOOTHING_CORRECTION
  /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float mj = pj->mass;

  const float ui = r / pi->h;

  pi->drho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->grad_rho[0] += dx[0] * wi_dx * r_inv * mj;
  pi->grad_rho[1] += dx[1] * wi_dx * r_inv * mj;
  pi->grad_rho[2] += dx[2] * wi_dx * r_inv * mj;
#endif /* PLANETARY_IMBALANCE */
}

/**
 * @brief Finishes extra density estimate parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_density_extra_density_estimate(struct part *restrict p) {

#ifdef PLANETARY_IMBALANCE
  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h; /* 1/h */

  /* Final operation on sum_wij (add self-contribution) */
  // p->sum_wij += sqrtf(kernel_root)*p->mass; // sqrt variation
  p->sum_wij += kernel_root * p->mass;  // nosqrt variation

  /* Compute norm sum_rij */
  float sum_rij_norm = 0.f;
  sum_rij_norm +=
      p->sum_rij[0] * p->sum_rij[0] * h_inv * h_inv / p->sum_wij / p->sum_wij;
  sum_rij_norm +=
      p->sum_rij[1] * p->sum_rij[1] * h_inv * h_inv / p->sum_wij / p->sum_wij;
  sum_rij_norm +=
      p->sum_rij[2] * p->sum_rij[2] * h_inv * h_inv / p->sum_wij / p->sum_wij;
  p->I = sqrtf(sum_rij_norm) * planetary_imbalance_alpha;

#elif PLANETARY_SMOOTHING_CORRECTION
  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  p->drho_dh -= p->mass * hydro_dimension * kernel_root;
  p->drho_dh *= h_inv_dim_plus_one;

  p->grad_rho[0] *= h_inv_dim_plus_one;
  p->grad_rho[1] *= h_inv_dim_plus_one;
  p->grad_rho[2] *= h_inv_dim_plus_one;
#endif /* PLANETARY_IMBALANCE */
}

/**
 * @brief Prepares extra density estimate parameters for a particle for the
 * gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_extra_density_estimate(struct part *restrict p) {

#ifdef PLANETARY_IMBALANCE
  /* Initialize kernel averages to 0 */
  p->sum_wij_exp_T = 0.f;
  p->sum_wij_exp_P = 0.f;
  p->sum_wij_exp = 0.f;

  // Compute the pressure
  const float pressure =
      gas_pressure_from_internal_energy(p->rho, p->u, p->mat_id);

  // Compute the temperature
  const float temperature =
      gas_temperature_from_internal_energy(p->rho, p->u, p->mat_id);

  p->P = pressure;
  p->T = temperature;

  /* Turn Imbalance to 0 if h == h_max */
  if (p->is_h_max) {
    p->I = 0.f;
  }

  // p->imbalance_flag = 0;
  // if (p->h < 0.999f * hydro_props->h_max){
  //      p->imbalance_flag = 1;
  //}
#elif PLANETARY_SMOOTHING_CORRECTION
  // Compute the pressure
  const float pressure =
      gas_pressure_from_internal_energy(p->rho, p->u, p->mat_id);

  p->P = pressure;

  p->P_tilde_numerator = 0.f;
  p->P_tilde_denominator = 0.f;

  p->max_ngb_sph_rho = p->rho;
  p->min_ngb_sph_rho = p->rho;

  p->grad_drho_dh[0] = 0.f;
  p->grad_drho_dh[1] = 0.f;
  p->grad_drho_dh[2] = 0.f;
#endif /* PLANETARY_IMBALANCE */
}

/**
 * @brief Extra gradient interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_gradient_extra_density_estimate(
    struct part *restrict pi, struct part *restrict pj, const float dx[3],
    const float wi, const float wj, const float wi_dx, const float wj_dx) {

#ifdef PLANETARY_IMBALANCE
  /* Compute kernel averages */
  pi->sum_wij_exp += wi * expf(-pj->I * pj->I);
  pi->sum_wij_exp_P += pj->P * wi * expf(-pj->I * pj->I);
  pi->sum_wij_exp_T += pj->T * wi * expf(-pj->I * pj->I);

  pj->sum_wij_exp += wj * expf(-pi->I * pi->I);
  pj->sum_wij_exp_P += pi->P * wj * expf(-pi->I * pi->I);
  pj->sum_wij_exp_T += pi->T * wj * expf(-pi->I * pi->I);
#elif PLANETARY_SMOOTHING_CORRECTION
  /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  float si = (pi->h / pi->rho) * sqrtf(pi->grad_rho[0] * pi->grad_rho[0] +
                                       pi->grad_rho[1] * pi->grad_rho[1] +
                                       pi->grad_rho[2] * pi->grad_rho[2]);
  float sj = (pj->h / pj->rho) * sqrtf(pj->grad_rho[0] * pj->grad_rho[0] +
                                       pj->grad_rho[1] * pj->grad_rho[1] +
                                       pj->grad_rho[2] * pj->grad_rho[2]);

  float f_gi = 1.f / (si + 0.001f);
  float f_gj = 1.f / (sj + 0.001f);

  pi->P_tilde_numerator += pj->P * f_gj * sqrtf(wi);
  pj->P_tilde_numerator += pi->P * f_gi * sqrtf(wj);
  pi->P_tilde_denominator += f_gj * sqrtf(wi);
  pj->P_tilde_denominator += f_gi * sqrtf(wj);

  pi->max_ngb_sph_rho = max(pi->max_ngb_sph_rho, pj->rho);
  pi->min_ngb_sph_rho = min(pi->min_ngb_sph_rho, pj->rho);
  pj->max_ngb_sph_rho = max(pj->max_ngb_sph_rho, pi->rho);
  pj->min_ngb_sph_rho = min(pj->min_ngb_sph_rho, pi->rho);

  pi->grad_drho_dh[0] += (pj->drho_dh - pi->drho_dh) * (dx[0] * wi_dx * r_inv) *
                         (pj->mass / pj->rho);
  pi->grad_drho_dh[1] += (pj->drho_dh - pi->drho_dh) * (dx[1] * wi_dx * r_inv) *
                         (pj->mass / pj->rho);
  pi->grad_drho_dh[2] += (pj->drho_dh - pi->drho_dh) * (dx[2] * wi_dx * r_inv) *
                         (pj->mass / pj->rho);

  pj->grad_drho_dh[0] += (pi->drho_dh - pj->drho_dh) *
                         (-dx[0] * wj_dx * r_inv) * (pi->mass / pi->rho);
  pj->grad_drho_dh[1] += (pi->drho_dh - pj->drho_dh) *
                         (-dx[1] * wj_dx * r_inv) * (pi->mass / pi->rho);
  pj->grad_drho_dh[2] += (pi->drho_dh - pj->drho_dh) *
                         (-dx[2] * wj_dx * r_inv) * (pi->mass / pi->rho);
#endif /* PLANETARY_IMBALANCE */
}

/**
 * @brief Extra gradient interaction between two particles (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_gradient_extra_density_estimate(
    struct part *restrict pi, const struct part *restrict pj, const float dx[3],
    const float wi, const float wi_dx) {

#ifdef PLANETARY_IMBALANCE
  /* Compute kernel averages */
  pi->sum_wij_exp += wi * expf(-pj->I * pj->I);
  pi->sum_wij_exp_P += pj->P * wi * expf(-pj->I * pj->I);
  pi->sum_wij_exp_T += pj->T * wi * expf(-pj->I * pj->I);
#elif PLANETARY_SMOOTHING_CORRECTION
  /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  float sj = (pj->h / pj->rho) * sqrtf(pj->grad_rho[0] * pj->grad_rho[0] +
                                       pj->grad_rho[1] * pj->grad_rho[1] +
                                       pj->grad_rho[2] * pj->grad_rho[2]);

  float f_gj = 1.f / (sj + 0.001f);

  pi->P_tilde_numerator += pj->P * f_gj * sqrtf(wi);
  pi->P_tilde_denominator += f_gj * sqrtf(wi);

  pi->max_ngb_sph_rho = max(pi->max_ngb_sph_rho, pj->rho);
  pi->min_ngb_sph_rho = min(pi->min_ngb_sph_rho, pj->rho);

  pi->grad_drho_dh[0] += (pj->drho_dh - pi->drho_dh) * (dx[0] * wi_dx * r_inv) *
                         (pj->mass / pj->rho);
  pi->grad_drho_dh[1] += (pj->drho_dh - pi->drho_dh) * (dx[1] * wi_dx * r_inv) *
                         (pj->mass / pj->rho);
  pi->grad_drho_dh[2] += (pj->drho_dh - pi->drho_dh) * (dx[2] * wi_dx * r_inv) *
                         (pj->mass / pj->rho);
#endif /* PLANETARY_IMBALANCE */
}

/**
 * @brief Finishes extra density estimate parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_extra_density_estimate(struct part *restrict p) {

#ifdef PLANETARY_IMBALANCE
  /* Add self contribution to kernel averages*/
  float I2 = p->I * p->I;
  p->sum_wij_exp += kernel_root * expf(-I2);
  p->sum_wij_exp_P += p->P * kernel_root * expf(-I2);
  p->sum_wij_exp_T += p->T * kernel_root * expf(-I2);

  /* compute minimum SPH quantities */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  const float rho_min = p->mass * kernel_root * h_inv_dim;

  if (p->sum_wij_exp > 0.f && p->sum_wij_exp_P > 0.f && p->I > 0.f) {
    /* End computation */
    p->sum_wij_exp_P /= p->sum_wij_exp;
    p->sum_wij_exp_T /= p->sum_wij_exp;

    /* Compute new P */
    float P_new = expf(-I2) * p->P + (1.f - expf(-I2)) * p->sum_wij_exp_P;

    /* Compute new T */
    float T_new = expf(-I2) * p->T + (1.f - expf(-I2)) * p->sum_wij_exp_T;

    /* Compute new density */
    float rho_new =
        gas_density_from_pressure_and_temperature(P_new, T_new, p->mat_id);

    if (rho_new > rho_min) {
      p->rho = rho_new;
    } else {
      p->rho = rho_min;
    }
  }

  // finish computations
  const float P = gas_pressure_from_internal_energy(p->rho, p->u, p->mat_id);
  const float T = gas_temperature_from_internal_energy(p->rho, p->u, p->mat_id);
  p->P = P;
  p->T = T;
#elif PLANETARY_SMOOTHING_CORRECTION
  /* compute minimum SPH quantities */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */
  const float rho_min = p->mass * kernel_root * h_inv_dim;

  p->grad_drho_dh[0] *= h_inv_dim_plus_one;
  p->grad_drho_dh[1] *= h_inv_dim_plus_one;
  p->grad_drho_dh[2] *= h_inv_dim_plus_one;

  float s = (p->h / p->rho) * sqrtf(p->grad_rho[0] * p->grad_rho[0] +
                                    p->grad_rho[1] * p->grad_rho[1] +
                                    p->grad_rho[2] * p->grad_rho[2]);

  float f_g = 1.f / (s + 0.001f);

  p->P_tilde_numerator += p->P * f_g * sqrtf(kernel_root);
  p->P_tilde_denominator += f_g * sqrtf(kernel_root);
  if (p->P_tilde_numerator == 0.f || p->P_tilde_denominator == 0.f) {
    p->P_tilde_numerator = p->P;
    p->P_tilde_denominator = 1.f;
  }
  float P_tilde = p->P_tilde_numerator / p->P_tilde_denominator;

  float grad_drho_dh = sqrtf(p->grad_drho_dh[0] * p->grad_drho_dh[0] +
                             p->grad_drho_dh[1] * p->grad_drho_dh[1] +
                             p->grad_drho_dh[2] * p->grad_drho_dh[2]);
  float S = (p->h / p->rho) * (fabs(p->drho_dh) + p->h * grad_drho_dh);
  float f_S = 0.5f * (1.f + tanhf(3.f - 3.f * S / (0.1f)));

  /* Turn S to 0 if h == h_max */
  if (f_S < 0.99f) {
    if (p->is_h_max) {
      S = 0.f;  // This is only for output files
      f_S = 1.f;
    } else {
      // Compute new P
      float P_new = f_S * p->P + (1.f - f_S) * P_tilde;

      float rho_ref;
      if (p->last_corrected_rho) {
        rho_ref = p->last_corrected_rho;
      } else {
        rho_ref = p->rho;
      }

      // Compute rho from u, P_new
      float rho_new_from_u = gas_density_from_pressure_and_internal_energy(
          P_new, p->u, rho_ref, p->rho, p->mat_id);

      if (rho_new_from_u > p->max_ngb_sph_rho) {
        rho_new_from_u = p->max_ngb_sph_rho;
      }
      if (rho_new_from_u < p->min_ngb_sph_rho) {
        rho_new_from_u = p->min_ngb_sph_rho;
      }

      // Ensure new density is not lower than minimum SPH density
      if (rho_new_from_u > rho_min) {
        p->rho = rho_new_from_u;
      } else {
        p->rho = rho_min;
      }
    }
  }

  p->smoothing_error = S;

  p->last_corrected_rho = p->rho;
  p->last_f_S = f_S;
#endif /* PLANETARY_IMBALANCE */
}

#endif /* SWIFT_PLANETARY_HYDRO_DENSITY_ESTIMATE_H */
