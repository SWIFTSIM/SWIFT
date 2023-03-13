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
#ifndef SWIFT_PLANETARY_HYDRO_DENSITY_EXTRA_H
#define SWIFT_PLANETARY_HYDRO_DENSITY_EXTRA_H

/**
 * @file Planetary/hydro_density_extra.h
 * @brief Modifications and additions for hydro density estimates.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Prepares extra parameters for a particle for the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_init_part_extra(
    struct part *restrict p) {

#ifdef PLANETARY_IMBALANCE
  p->P = 0.f;
  p->T = 0.f;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  p->N_density = 1; /* Self contribution */
  p->N_force = 0;
  p->N_density_exact = 0;
  p->N_force_exact = 0;
  p->rho_exact = 0.f;
  p->n_density = 0.f;
  p->n_density_exact = 0.f;
  p->n_force = 0.f;
  p->n_force_exact = 0.f;
  p->inhibited_exact = 0;
  p->limited_part = 0;
#endif

#ifdef PLANETARY_IMBALANCE
  p->sum_rij[0] = 0.f;
  p->sum_rij[1] = 0.f;
  p->sum_rij[2] = 0.f;
  p->I = 0.f;
  p->sum_wij = 0.f;
#endif

#ifdef PLANETARY_SMOOTHING_CORRECTION
  p->drho_dh = 0.f;
  p->grad_rho[0] = 0.f;
  p->grad_rho[1] = 0.f;
  p->grad_rho[2] = 0.f;
#endif

#ifdef PLANETARY_QUAD_VISC
  int i, j;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      p->Dinv[i][j] = 0.f;
      p->E[i][j] = 0.f;
    }
  }
#endif
}

/**
 * @brief Finishes extra parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_end_density_extra(
    struct part *restrict p) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

#ifdef PLANETARY_IMBALANCE
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
#endif

#ifdef PLANETARY_SMOOTHING_CORRECTION
  p->drho_dh -= p->mass * hydro_dimension * kernel_root;
  p->drho_dh *= h_inv_dim_plus_one;

  p->grad_rho[0] *= h_inv_dim_plus_one;
  p->grad_rho[1] *= h_inv_dim_plus_one;
  p->grad_rho[2] *= h_inv_dim_plus_one;
#endif

#ifdef PLANETARY_QUAD_VISC
  int i, j, k;

  /* In this section we carry out matrix inversion to find D and calculate
   * dv_aux */

  float determinant = 0.f;
  /* We normalise the Dinv matrix to the mean of its 9 elements to stop us
   * hitting float precision limits during matrix inversion process */
  float mean_Dinv = (p->Dinv[0][0] + p->Dinv[0][1] + p->Dinv[0][2] +
                     p->Dinv[1][0] + p->Dinv[1][1] + p->Dinv[1][2] +
                     p->Dinv[2][0] + p->Dinv[2][1] + p->Dinv[2][2]) /
                    9.f;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Normalise Dinv to mean of its values */
      p->Dinv[i][j] = p->Dinv[i][j] / mean_Dinv;

      /* Aux dv (eq 19 in Rosswog 2020) */
      p->dv_aux[i][j] = 0.f;
    }
  }

  for (i = 0; i < 3; i++) {
    /* Matrix Dinv det */
    determinant +=
        (p->Dinv[0][i] * (p->Dinv[1][(i + 1) % 3] * p->Dinv[2][(i + 2) % 3] -
                          p->Dinv[1][(i + 2) % 3] * p->Dinv[2][(i + 1) % 3]));
  }

  float D[3][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Find D from inverse of Dinv */
      D[i][j] = ((p->Dinv[(i + 1) % 3][(j + 1) % 3] *
                  p->Dinv[(i + 2) % 3][(j + 2) % 3]) -
                 (p->Dinv[(i + 1) % 3][(j + 2) % 3] *
                  p->Dinv[(i + 2) % 3][(j + 1) % 3])) /
                (determinant * mean_Dinv);
      if (isnan(D[i][j]) || isinf(D[i][j])) {
        D[i][j] = 0.f;
      }

      for (k = 0; k < 3; ++k) {
        /* Calculate dv_aux (eq 19 in Rosswog 2020) */
        p->dv_aux[i][k] += D[i][j] * p->E[k][j];
      }
    }
  }
#endif
}

#endif /* SWIFT_PLANETARY_HYDRO_DENSITY_EXTRA_H */
