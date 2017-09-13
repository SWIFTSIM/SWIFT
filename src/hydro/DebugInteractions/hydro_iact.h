/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DEBUG_INTERACTIONS_HYDRO_IACT_H
#define SWIFT_DEBUG_INTERACTIONS_HYDRO_IACT_H

/**
 * @file DebugInteractions/hydro_iact.h
 * @brief Empty SPH implementation used solely to test the SELF/PAIR routines.
 */

/**
 * @brief Density loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;
  float wj, wj_dx;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  /* Compute the kernel function for pi */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute the kernel function for pj */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Compute contribution to the number of neighbours */
  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  /* Update ngb counters */
  pi->ids_ngbs_density[pi->num_ngb_density] = pj->id;
  ++pi->num_ngb_density;
  pj->ids_ngbs_density[pj->num_ngb_density] = pi->id;
  ++pj->num_ngb_density;
}

/**
 * @brief Density loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Update ngb counters */
  pi->ids_ngbs_density[pi->num_ngb_density] = pj->id;
  ++pi->num_ngb_density;
}

/**
 * @brief Force loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* Update ngb counters */
  pi->ids_ngbs_force[pi->num_ngb_force] = pj->id;
  ++pi->num_ngb_force;
  pj->ids_ngbs_force[pj->num_ngb_force] = pi->id;
  ++pj->num_ngb_force;
}

/**
 * @brief Force loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* Update ngb counters */
  pi->ids_ngbs_force[pi->num_ngb_force] = pj->id;
  ++pi->num_ngb_force;
}

#endif /* SWIFT_DEBUG_INTERACTIONS_HYDRO_IACT_H */
