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
#ifndef SWIFT_PLANETARY_HYDRO_KERNELS_H
#define SWIFT_PLANETARY_HYDRO_KERNELS_H

/**
 * @file Planetary/hydro_kernels.h
 * @brief Utilities for REMIX hydro kernels.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Prepares extra kernel parameters for a particle for the density
 * calculation.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part_extra_kernel(
    struct part *restrict p) {

  p->m0 = 0.f;

  p->grad_m0[0] = 0.f;
  p->grad_m0[1] = 0.f;
  p->grad_m0[2] = 0.f;
}

/**
 * @brief Extra kernel density interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_density_extra_kernel(struct part *restrict pi,
                                       struct part *restrict pj,
                                       const float dx[3], const float wi,
                                       const float wj, const float wi_dx,
                                       const float wj_dx) {

  /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Geometric moments and gradients (unmodified kernel, for normalisation) */
  pi->m0 += pj->mass * wi / pj->rho_evol;
  pj->m0 += pi->mass * wj / pi->rho_evol;
  for (int i = 0; i < 3; i++) {
    pi->grad_m0[i] += (pj->mass / pj->rho_evol) * dx[i] * wi_dx * r_inv;
    pj->grad_m0[i] += -(pi->mass / pi->rho_evol) * dx[i] * wj_dx * r_inv;
  }
}

/**
 * @brief Extra kernel density interaction between two particles (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_density_extra_kernel(struct part *restrict pi,
                                              const struct part *restrict pj,
                                              const float dx[3], const float wi,
                                              const float wi_dx) {

  /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Geometric moments and gradients (unmodified kernel, for normalisation) */
  pi->m0 += pj->mass * wi / pj->rho_evol;
  for (int i = 0; i < 3; i++) {
    pi->grad_m0[i] += (pj->mass / pj->rho_evol) * dx[i] * wi_dx * r_inv;
  }
}

/**
 * @brief Finishes extra kernel parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_density_extra_kernel(struct part *restrict p) {

  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Geometric moments and gradients (unmodified kernel, for normalisation) */
  p->m0 += p->mass * kernel_root / p->rho_evol;
  p->m0 *= h_inv_dim;
  for (int i = 0; i < 3; i++) {
    p->grad_m0[i] *= h_inv_dim_plus_one;
  }
}

/**
 * @brief Prepares extra kernel parameters for a particle for the gradient
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_extra_kernel(struct part *restrict p) {

  /* Initialise geometric moment matrices (ij-mean, for linear-order repr.
   * kernel) */
  zero_sym_matrix(&p->gradient.m2_bar);
  zero_sym_matrix(&p->gradient.grad_m2_bar[0]);
  zero_sym_matrix(&p->gradient.grad_m2_bar[1]);
  zero_sym_matrix(&p->gradient.grad_m2_bar[2]);
  zero_sym_matrix(&p->gradient.grad_m2_bar_gradhterm);

  /* Geometric moments and gradients (ij-mean, for linear-order repr. kernel) */
  p->gradient.m0_bar = 0.f;
  p->gradient.grad_m0_bar_gradhterm = 0.f;
  for (int i = 0; i < 3; i++) {
    p->gradient.m1_bar[i] = 0.f;
    p->gradient.grad_m0_bar[i] = 0.f;
    p->gradient.grad_m1_bar_gradhterm[i] = 0.f;
    for (int j = 0; j < 3; j++) {
      p->gradient.grad_m1_bar[i][j] = 0.f;
    }
  }
}

/**
 * @brief Extra kernel gradient interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_gradient_extra_kernel(struct part *restrict pi,
                                        struct part *restrict pj,
                                        const float dx[3], const float wi,
                                        const float wj, const float wi_dx,
                                        const float wj_dx) {

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

  // Volume elements
  float volume_i = pi->mass / pi->rho_evol;
  float volume_j = pj->mass / pj->rho_evol;

  // Mean ij kernels and gradients
  float wi_term = 0.5f * (wi * hi_inv_dim + wj * hj_inv_dim);
  float wj_term = wi_term;
  float wi_dx_term[3], wj_dx_term[3];
  for (int i = 0; i < 3; i++) {
    wi_dx_term[i] = dx[i] * r_inv * 0.5f *
                    (wi_dx * hi_inv_dim_plus_one + wj_dx * hj_inv_dim_plus_one);
    wj_dx_term[i] = -dx[i] * r_inv * 0.5f *
                    (wi_dx * hi_inv_dim_plus_one + wj_dx * hj_inv_dim_plus_one);
  }

  // Grad-h term, dW/dh
  float wi_dx_gradhterm =
      -0.5f * (hydro_dimension * wi + (r / hi) * wi_dx) * hi_inv_dim_plus_one;
  float wj_dx_gradhterm =
      -0.5f * (hydro_dimension * wj + (r / hj) * wj_dx) * hj_inv_dim_plus_one;

  // Geometric moments m_0, m_1, and m_2, their gradients and grad-h terms
  pi->gradient.m0_bar += wi_term * volume_j;
  pj->gradient.m0_bar += wj_term * volume_i;

  pi->gradient.grad_m0_bar_gradhterm += wi_dx_gradhterm * volume_j;
  pj->gradient.grad_m0_bar_gradhterm += wj_dx_gradhterm * volume_i;
  for (int i = 0; i < 3; i++) {
    pi->gradient.m1_bar[i] += dx[i] * wi_term * volume_j;
    pj->gradient.m1_bar[i] += -dx[i] * wj_term * volume_i;

    pi->gradient.grad_m0_bar[i] += wi_dx_term[i] * volume_j;
    pj->gradient.grad_m0_bar[i] += wj_dx_term[i] * volume_i;

    pi->gradient.grad_m1_bar_gradhterm[i] += dx[i] * wi_dx_gradhterm * volume_j;
    pj->gradient.grad_m1_bar_gradhterm[i] +=
        -dx[i] * wj_dx_gradhterm * volume_i;

    for (int j = 0; j < 3; j++) {

      pi->gradient.grad_m1_bar[i][j] += dx[j] * wi_dx_term[i] * volume_j;
      pj->gradient.grad_m1_bar[i][j] += -dx[j] * wj_dx_term[i] * volume_i;
    }
  }

  for (int j = 0; j < 3; j++) {
    for (int k = j; k < 3; k++) {
      int i = (j == k) ? j : 2 + j + k;
      pi->gradient.m2_bar.elements[i] += dx[j] * dx[k] * wi_term * volume_j;
      pj->gradient.m2_bar.elements[i] += dx[j] * dx[k] * wj_term * volume_i;

      pi->gradient.grad_m2_bar_gradhterm.elements[i] +=
          dx[j] * dx[k] * wi_dx_gradhterm * volume_j;
      pj->gradient.grad_m2_bar_gradhterm.elements[i] +=
          dx[j] * dx[k] * wj_dx_gradhterm * volume_i;

      for (int l = 0; l < 3; l++) {
        pi->gradient.grad_m2_bar[l].elements[i] +=
            dx[j] * dx[k] * wi_dx_term[l] * volume_j;
        pj->gradient.grad_m2_bar[l].elements[i] +=
            dx[j] * dx[k] * wj_dx_term[l] * volume_i;
      }
    }
  }
}

/**
 * @brief Extra kernel gradient interaction between two particles
 * (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_gradient_extra_kernel(
    struct part *restrict pi, struct part *restrict pj, const float dx[3],
    const float wi, const float wj, const float wi_dx, const float wj_dx) {

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

  // Volume elements
  float volume_j = pj->mass / pj->rho_evol;

  // Mean ij kernel and gradients
  float wi_term = 0.5f * (wi * hi_inv_dim + wj * hj_inv_dim);
  float wi_dx_term[3];
  for (int i = 0; i < 3; i++) {
    wi_dx_term[i] = dx[i] * r_inv * 0.5f *
                    (wi_dx * hi_inv_dim_plus_one + wj_dx * hj_inv_dim_plus_one);
  }

  // Grad-h term, dW/dh
  float wi_dx_gradhterm =
      -0.5f * (hydro_dimension * wi + (r / hi) * wi_dx) * hi_inv_dim_plus_one;

  // Geometric moments m_0, m_1, and m_2, their gradients and grad-h terms
  pi->gradient.m0_bar += wi_term * volume_j;

  pi->gradient.grad_m0_bar_gradhterm += wi_dx_gradhterm * volume_j;
  for (int i = 0; i < 3; i++) {
    pi->gradient.m1_bar[i] += dx[i] * wi_term * volume_j;

    pi->gradient.grad_m0_bar[i] += wi_dx_term[i] * volume_j;

    pi->gradient.grad_m1_bar_gradhterm[i] += dx[i] * wi_dx_gradhterm * volume_j;

    for (int j = 0; j < 3; j++) {

      pi->gradient.grad_m1_bar[i][j] += dx[j] * wi_dx_term[i] * volume_j;
    }
  }

  for (int j = 0; j < 3; j++) {
    for (int k = j; k < 3; k++) {
      int i = (j == k) ? j : 2 + j + k;
      pi->gradient.m2_bar.elements[i] += dx[j] * dx[k] * wi_term * volume_j;

      pi->gradient.grad_m2_bar_gradhterm.elements[i] +=
          dx[j] * dx[k] * wi_dx_gradhterm * volume_j;

      for (int l = 0; l < 3; l++) {
        pi->gradient.grad_m2_bar[l].elements[i] +=
            dx[j] * dx[k] * wi_dx_term[l] * volume_j;
      }
    }
  }
}

/**
 * @brief Finishes extra kernel parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_extra_kernel(struct part *restrict p) {

  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  // Volume elements
  float volume = p->mass / p->rho_evol;

  // Self contribution to geometric moments and gradients
  p->gradient.m0_bar += volume * kernel_root * h_inv_dim;
  p->gradient.grad_m0_bar_gradhterm -=
      0.5f * volume * hydro_dimension * kernel_root * h_inv_dim_plus_one;

  // Multiply dh/dr into grad-h terms
  for (int i = 0; i < 3; i++) {
    p->gradient.grad_m0_bar[i] +=
        p->gradient.grad_m0_bar_gradhterm * p->dh_norm_kernel[i];
    for (int j = 0; j < 3; j++) {
      p->gradient.grad_m1_bar[i][j] +=
          p->gradient.grad_m1_bar_gradhterm[j] * p->dh_norm_kernel[i];
    }
    for (int k = 0; k < 6; k++) {
      p->gradient.grad_m2_bar[i].elements[k] +=
          p->gradient.grad_m2_bar_gradhterm.elements[k] * p->dh_norm_kernel[i];
    }
  }
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_force_extra_kernel(struct part *restrict p) {

  if (!p->is_h_max) {

    // Add second terms to complete the geometric moment gradients
    for (int i = 0; i < 3; i++) {
      p->gradient.grad_m1_bar[i][i] += p->gradient.m0_bar;
    }

    p->gradient.grad_m2_bar[0].xx += 2.f * p->gradient.m1_bar[0];
    p->gradient.grad_m2_bar[0].xy += p->gradient.m1_bar[1];
    p->gradient.grad_m2_bar[0].xz += p->gradient.m1_bar[2];

    p->gradient.grad_m2_bar[1].yy += 2.f * p->gradient.m1_bar[1];
    p->gradient.grad_m2_bar[1].xy += p->gradient.m1_bar[0];
    p->gradient.grad_m2_bar[1].yz += p->gradient.m1_bar[2];

    p->gradient.grad_m2_bar[2].zz += 2.f * p->gradient.m1_bar[2];
    p->gradient.grad_m2_bar[2].xz += p->gradient.m1_bar[0];
    p->gradient.grad_m2_bar[2].yz += p->gradient.m1_bar[1];

    // Inverse of symmetric geometric moment m_2 (bar)
    struct sym_matrix m2_bar_inv;
    // Make m2_bar dimensionless for calculation of inverse
    struct sym_matrix m2_bar_over_h2;
    for (int i = 0; i < 6; i++) m2_bar_over_h2.elements[i] =
                                p->gradient.m2_bar.elements[i] / (p->h * p->h);
    sym_matrix_invert(&m2_bar_inv, &m2_bar_over_h2);
    for (int i = 0; i < 6; i++) m2_bar_inv.elements[i] /= (p->h * p->h);

    // Components for constructing linear-order kernel's A and B, and gradients
    float m2_bar_inv_mult_m1_bar[3];
    float m2_bar_inv_mult_grad_m1_bar[3][3];
    struct sym_matrix m2_bar_inv_mult_grad_m2_bar_mult_m2_bar_inv[3];
    float ABA_mult_m1_bar[3][3];

    sym_matrix_multiply_by_vector(m2_bar_inv_mult_m1_bar, &m2_bar_inv,
                                  p->gradient.m1_bar);
    for (int i = 0; i < 3; i++) {
      sym_matrix_multiply_by_vector(m2_bar_inv_mult_grad_m1_bar[i], &m2_bar_inv,
                                    p->gradient.grad_m1_bar[i]);
      sym_matrix_multiplication_ABA(
          &m2_bar_inv_mult_grad_m2_bar_mult_m2_bar_inv[i], &m2_bar_inv,
          &p->gradient.grad_m2_bar[i]);
      sym_matrix_multiply_by_vector(
          ABA_mult_m1_bar[i], &m2_bar_inv_mult_grad_m2_bar_mult_m2_bar_inv[i],
          p->gradient.m1_bar);
    }

    // Linear-order reproducing kernel's A and B components and gradients
    float A, B[3], grad_A[3], grad_B[3][3];

    // Calculate A
    A = p->gradient.m0_bar;
    for (int i = 0; i < 3; i++) {
      A -= m2_bar_inv_mult_m1_bar[i] * p->gradient.m1_bar[i];
    }
    A = 1.f / A;

    // Calculate B
    for (int i = 0; i < 3; i++) {
      B[i] = -m2_bar_inv_mult_m1_bar[i];
    }

    // Calculate grad A
    for (int i = 0; i < 3; i++) {
      grad_A[i] = p->gradient.grad_m0_bar[i];

      for (int j = 0; j < 3; j++) {
        grad_A[i] +=
            -2 * m2_bar_inv_mult_m1_bar[j] * p->gradient.grad_m1_bar[i][j] +
            ABA_mult_m1_bar[i][j] * p->gradient.m1_bar[j];
      }

      grad_A[i] *= -A * A;
    }

    // Calculate grad B
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        grad_B[j][i] =
            -m2_bar_inv_mult_grad_m1_bar[j][i] + ABA_mult_m1_bar[j][i];
      }
    }

    // Store final values
    p->force.A = A;
    for (int i = 0; i < 3; i++) {
      p->force.B[i] = B[i];
      p->force.grad_A[i] = grad_A[i];
      for (int j = 0; j < 3; j++) {
        p->force.grad_B[i][j] = grad_B[i][j];
      }
    }

    // Vacuum-boundary proximity switch
    p->force.vac_switch = 1.f;
    float hB = p->h * sqrtf(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
    float offset = 0.8f;

    if (hB > offset) {
      float sigma = 0.2f;
      p->force.vac_switch =
          expf(-(hB - offset) * (hB - offset) / (2.f * sigma * sigma));
    }

  } else {
    // p->is_h_max
    p->force.A = 1.f;
    p->force.vac_switch = 1.f;
    for (int i = 0; i < 3; i++) {
      p->force.B[i] = 0.f;
      p->force.grad_A[i] = 0.f;
      for (int j = 0; j < 3; j++) {
        p->force.grad_B[i][j] = 0.f;
      }
    }
  }
}

/**
 * @brief Set gradient terms for linear-order reproducing kernel.
 *
 * Note `G` here corresponds to `d/dr tilde{mathcal{W}}` in Sandnes et al.
 * (2024)
 */
__attribute__((always_inline)) INLINE static void hydro_set_Gi_Gj_forceloop(
    float Gi[3], float Gj[3], const struct part *restrict pi,
    const struct part *restrict pj, const float dx[3], const float wi,
    const float wj, const float wi_dx, const float wj_dx) {

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

  const float wi_dr = hi_inv_dim_plus_one * wi_dx;
  const float wj_dr = hj_inv_dim_plus_one * wj_dx;

  if ((!pi->is_h_max) && (!pj->is_h_max)) {
    // Mean ij kernels and gradients with grad-h terms
    float wi_term = 0.5f * (wi * hi_inv_dim + wj * hj_inv_dim);
    float wj_term = wi_term;
    float wi_dx_term[3], wj_dx_term[3];
    float Ai = pi->force.A;
    float Aj = pj->force.A;
    float grad_Ai[3], grad_Aj[3], Bi[3], Bj[3], grad_Bi[3][3], grad_Bj[3][3];
    for (int i = 0; i < 3; i++) {
      grad_Ai[i] = pi->force.grad_A[i];
      grad_Aj[i] = pj->force.grad_A[i];
      Bi[i] = pi->force.B[i];
      Bj[i] = pj->force.B[i];
      for (int j = 0; j < 3; j++) {
        grad_Bi[i][j] = pi->force.grad_B[i][j];
        grad_Bj[i][j] = pj->force.grad_B[i][j];
      }
    }

    for (int i = 0; i < 3; i++) {
      wi_dx_term[i] =
          dx[i] * r_inv * 0.5f *
          (wi_dx * hi_inv_dim_plus_one + wj_dx * hj_inv_dim_plus_one);
      wj_dx_term[i] =
          -dx[i] * r_inv * 0.5f *
          (wi_dx * hi_inv_dim_plus_one + wj_dx * hj_inv_dim_plus_one);

      wi_dx_term[i] += -0.5f * (hydro_dimension * wi + (r / pi->h) * wi_dx) *
                       hi_inv_dim_plus_one * pi->dh_norm_kernel[i];
      wj_dx_term[i] += -0.5f * (hydro_dimension * wj + (r / pj->h) * wj_dx) *
                       hj_inv_dim_plus_one * pj->dh_norm_kernel[i];

      Gi[i] = Ai * wi_dx_term[i] + grad_Ai[i] * wi_term + Ai * Bi[i] * wi_term;
      Gj[i] = Aj * wj_dx_term[i] + grad_Aj[i] * wj_term + Aj * Bj[i] * wj_term;

      for (int j = 0; j < 3; j++) {
        Gi[i] += Ai * Bi[j] * dx[j] * wi_dx_term[i] +
                 grad_Ai[i] * Bi[j] * dx[j] * wi_term +
                 Ai * grad_Bi[i][j] * dx[j] * wi_term;

        Gj[i] += -(Aj * Bj[j] * dx[j] * wj_dx_term[i] +
                   grad_Aj[i] * Bj[j] * dx[j] * wj_term +
                   Aj * grad_Bj[i][j] * dx[j] * wj_term);
      }
    }

    // Standard-kernel gradients, to be used for vacuum boundary switch
    float wi_dx_term_vac[3], wj_dx_term_vac[3];
    for (int i = 0; i < 3; i++) {
      wi_dx_term_vac[i] = dx[i] * r_inv * wi_dx * hi_inv_dim_plus_one;
      wj_dx_term_vac[i] = -dx[i] * r_inv * wj_dx * hj_inv_dim_plus_one;

      wi_dx_term_vac[i] += -(hydro_dimension * wi + (r / pi->h) * wi_dx) *
                           hi_inv_dim_plus_one * pi->dh_norm_kernel[i];
      wj_dx_term_vac[i] += -(hydro_dimension * wj + (r / pj->h) * wj_dx) *
                           hj_inv_dim_plus_one * pj->dh_norm_kernel[i];
    }

    // Gradients, including vacuum boundary switch
    for (int i = 0; i < 3; i++) {
      Gi[i] = pi->force.vac_switch * Gi[i] +
              (1.f - pi->force.vac_switch) * wi_dx_term_vac[i];
      Gj[i] = pj->force.vac_switch * Gj[i] +
              (1.f - pj->force.vac_switch) * wj_dx_term_vac[i];
    }

  } else {
    // One or both particles at h_max
    for (int i = 0; i < 3; i++) {
      Gi[i] = wi_dr * dx[i] * r_inv;
      Gj[i] = -wj_dr * dx[i] * r_inv;
    }
  }
}

#endif /* SWIFT_PLANETARY_HYDRO_KERNELS_H */
