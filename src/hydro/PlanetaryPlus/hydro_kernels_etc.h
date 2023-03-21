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
#ifndef SWIFT_PLANETARY_HYDRO_KERNELS_ETC_H
#define SWIFT_PLANETARY_HYDRO_KERNELS_ETC_H

/**
 * @file Planetary/hydro_kernels_etc.h
 * @brief Utilities for hydro kernels and e.g. gradient estimates.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Prepares extra kernel parameters for a particle for the density
 * calculation.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part_extra_kernel(
    struct part *restrict p) {}

/**
 * @brief Extra kernel density interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_density_extra_kernel(struct part *restrict pi,
                                       struct part *restrict pj,
                                       const float dx[3], const float wi,
                                       const float wj, const float wi_dx,
                                       const float wj_dx) {}

/**
 * @brief Extra kernel density interaction between two particles (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_density_extra_kernel(struct part *restrict pi,
                                              const struct part *restrict pj,
                                              const float dx[3], const float wi,
                                              const float wi_dx) {}

/**
 * @brief Finishes extra kernel parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_density_extra_kernel(struct part *restrict p) {}

/**
 * @brief Prepares extra kernel parameters for a particle for the gradient
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_extra_kernel(struct part *restrict p) {

#if defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC
  int i, j;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      p->Cinv[i][j] = 0.f;
    }
  }
#endif
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

#if defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC

  float volume_i = pi->mass / pi->rho;
  float volume_j = pj->mass / pj->rho;

  planetary_smoothing_correction_tweak_volume(&volume_i, pi);
  planetary_smoothing_correction_tweak_volume(&volume_j, pj);

  int i, j;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      /* Inverse of C matrix (eq 6 in Rosswog 2020) */
      pi->Cinv[i][j] += dx[i] * dx[j] * wi * volume_j;
      pj->Cinv[i][j] += dx[i] * dx[j] * wj * volume_i;
    }
  }

#if defined(HYDRO_DIMENSION_2D)
  /* This is so we can do 3x3 matrix inverse even when 2D */
  pi->Cinv[2][2] = 1.f;
  pj->Cinv[2][2] = 1.f;
#endif

#endif /* defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC */
}

/**
 * @brief Extra kernel gradient interaction between two particles
 * (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_gradient_extra_kernel(struct part *restrict pi,
                                               const struct part *restrict pj,
                                               const float dx[3],
                                               const float wi,
                                               const float wi_dx) {

#if defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC

  float volume_j = pj->mass / pj->rho;

  planetary_smoothing_correction_tweak_volume(&volume_j, pj);

  int i, j;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      /* Inverse of C matrix (eq 6 in Rosswog 2020) */
      pi->Cinv[i][j] += dx[i] * dx[j] * wi * volume_j;
    }
  }

#if defined(HYDRO_DIMENSION_2D)
  /* This is so we can do 3x3 matrix inverse even when 2D */
  pi->Cinv[2][2] = 1.f;
#endif

#endif /* defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC */
}

/**
 * @brief Finishes extra kernel parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_extra_kernel(struct part *restrict p) {

#if defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC
  int i, j;

  /* Find the inverse of the Cinv matrix */

  /* If h=h_max don't do anything fancy. Things like using m/rho to calculate
   * the volume stops working */

  if (!p->is_h_max) {

    float determinant = 0.f;
    /* We normalise the Cinv matrix to the mean of its 9 elements to stop us
     * hitting float precision limits during matrix inversion process */
    float mean_Cinv = (p->Cinv[0][0] + p->Cinv[0][1] + p->Cinv[0][2] +
                       p->Cinv[1][0] + p->Cinv[1][1] + p->Cinv[1][2] +
                       p->Cinv[2][0] + p->Cinv[2][1] + p->Cinv[2][2]) /
                      9.f;

    /* Calculate det and normalise Cinv */
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        p->Cinv[i][j] = p->Cinv[i][j] / mean_Cinv;
      }
    }

    for (i = 0; i < 3; i++) {
      determinant +=
          (p->Cinv[0][i] * (p->Cinv[1][(i + 1) % 3] * p->Cinv[2][(i + 2) % 3] -
                            p->Cinv[1][(i + 2) % 3] * p->Cinv[2][(i + 1) % 3]));
    }

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        /* Find C from inverse of Cinv */
        p->C[i][j] = ((p->Cinv[(i + 1) % 3][(j + 1) % 3] *
                       p->Cinv[(i + 2) % 3][(j + 2) % 3]) -
                      (p->Cinv[(i + 1) % 3][(j + 2) % 3] *
                       p->Cinv[(i + 2) % 3][(j + 1) % 3])) /
                     (determinant * mean_Cinv);
        if (isnan(p->C[i][j]) || isinf(p->C[i][j])) {
          p->C[i][j] = 0.f;
          // printf("C error"); //##
          // exit(0);
        }
      }
    }

  } else {

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        p->C[i][j] = 0.f;
      }
    }
  }
#endif
}

/**
 * @brief Returns particle Gs, equivalent to kernel gradients
 */
__attribute__((always_inline)) INLINE static void hydro_set_Gi_Gj(
    float Gi[3], float Gj[3], const struct part *restrict pi,
    const struct part *restrict pj, const float dx[3], const float wi,
    const float wj, const float wi_dx, const float wj_dx) {

  /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float hi_inv = 1.0f / pi->h;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */

  const float hj_inv = 1.0f / pj->h;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */

  const float wi_dr = hid_inv * wi_dx;
  const float wj_dr = hjd_inv * wj_dx;

#ifdef PLANETARY_MATRIX_INVERSION
  if (!pi->is_h_max && !pj->is_h_max) {
    for (int i = 0; i < 3; ++i) {
      /* eq 4 and 5 in Rosswog 2020. These replace the gradient of the kernel */
      Gi[i] =
          -(pi->C[i][0] * dx[0] + pi->C[i][1] * dx[1] + pi->C[i][2] * dx[2]) *
          wi;
      Gj[i] =
          -(pj->C[i][0] * dx[0] + pj->C[i][1] * dx[1] + pj->C[i][2] * dx[2]) *
          wj;
    }
  } else {
    for (int i = 0; i < 3; ++i) {
      /* If h=h_max use the standard kernel gradients */
      Gi[i] = wi_dr * dx[i] * r_inv;
      Gj[i] = wj_dr * dx[i] * r_inv;
    }
  }
#elseif PLANETARY_GDF
  /* Standard GDF kernel gradients, Wadsley+2017 Eqn. 7, in Rosswog2020
   * framework */
  /* Include the dx and r_inv here instead of later */
  for (int i = 0; i < 3; i++) {
    Gi[i] = wi_dr * dx[i] * r_inv * pi->f_gdf;
    Gj[i] = wj_dr * dx[i] * r_inv * pj->f_gdf;
  }
#else  /* !PLANETARY_MATRIX_INVERSION, !PLANETARY_GDF */
  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / pj->mass;
  const float f_ji = 1.f - pj->force.f / pi->mass;

  for (int i = 0; i < 3; i++) {
    Gi[i] = wi_dr * dx[i] * r_inv * f_ij;
    Gj[i] = wj_dr * dx[i] * r_inv * f_ji;
  }
#endif /* PLANETARY_MATRIX_INVERSION */
}

/**
 * @brief Returns kernel gradient terms used in evolution equations
 */
__attribute__((always_inline)) INLINE static void
hydro_set_kernel_gradient_terms(float kernel_gradient_i[3],
                                float kernel_gradient_j[3],
                                float Q_kernel_gradient_i[3],
                                float Q_kernel_gradient_j[3], const float Gi[3],
                                const float Gj[3]) {

#ifdef PLANETARY_GDF
  /* In GDF we use average of Gi and Gj. */
  kernel_gradient_i[0] = 0.5f * (Gi[0] + Gj[0]);
  kernel_gradient_i[1] = 0.5f * (Gi[1] + Gj[1]);
  kernel_gradient_i[2] = 0.5f * (Gi[2] + Gj[2]);

  kernel_gradient_j[0] = 0.5f * (Gi[0] + Gj[0]);
  kernel_gradient_j[1] = 0.5f * (Gi[1] + Gj[1]);
  kernel_gradient_j[2] = 0.5f * (Gi[2] + Gj[2]);
#else  /* !PLANETARY_GDF */
  kernel_gradient_i[0] = Gi[0];
  kernel_gradient_i[1] = Gi[1];
  kernel_gradient_i[2] = Gi[2];

  kernel_gradient_j[0] = Gj[0];
  kernel_gradient_j[1] = Gj[1];
  kernel_gradient_j[2] = Gj[2];
#endif /* PLANETARY_GDF */

#ifdef PLANETARY_QUAD_VISC
  Q_kernel_gradient_i[0] = kernel_gradient_i[0];
  Q_kernel_gradient_i[1] = kernel_gradient_i[1];
  Q_kernel_gradient_i[2] = kernel_gradient_i[2];

  Q_kernel_gradient_j[0] = kernel_gradient_j[0];
  Q_kernel_gradient_j[1] = kernel_gradient_j[1];
  Q_kernel_gradient_j[2] = kernel_gradient_j[2];
#else  /* !PLANETARY_QUAD_VISC */
  Q_kernel_gradient_i[0] = 0.5f * (kernel_gradient_i[0] + kernel_gradient_j[0]);
  Q_kernel_gradient_i[1] = 0.5f * (kernel_gradient_i[1] + kernel_gradient_j[1]);
  Q_kernel_gradient_i[2] = 0.5f * (kernel_gradient_i[2] + kernel_gradient_j[2]);

  Q_kernel_gradient_j[0] = 0.5f * (kernel_gradient_i[0] + kernel_gradient_j[0]);
  Q_kernel_gradient_j[1] = 0.5f * (kernel_gradient_i[1] + kernel_gradient_j[1]);
  Q_kernel_gradient_j[2] = 0.5f * (kernel_gradient_i[2] + kernel_gradient_j[2]);
#endif /* PLANETARY_QUAD_VISC */
}

#endif /* SWIFT_PLANETARY_HYDRO_KERNELS_ETC_H */
