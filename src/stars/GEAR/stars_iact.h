/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_STARS_IACT_H
#define SWIFT_GEAR_STARS_IACT_H

/* Local includes */
#include "hydro.h"
#include "feedback.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_density(const float r2, const float *dx,
                                 const float hi, const float hj,
                                 struct spart *restrict si,
                                 const struct part *restrict pj, const float a,
                                 const float H) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  si->density.wcount += wi;
  si->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

#ifdef DEBUG_INTERACTIONS_STARS
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_STARS)
    si->ids_ngbs_density[si->num_ngb_density] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_density;
#endif

  /* Accumulate the sum in the numerator and denominator of f_plus and f_minus */
  double dx_ij_plus[3];
  double dx_ij_minus[3];
  double scalar_weight_j = feedback_compute_scalar_weight(r2, dx, hi, hj, si, pj,
							  dx_ij_plus, dx_ij_minus);

  si->feedback_data.f_plus_num[0] += scalar_weight_j*fabs(dx_ij_minus[0]);
  si->feedback_data.f_plus_num[1] += scalar_weight_j*fabs(dx_ij_minus[1]);
  si->feedback_data.f_plus_num[2] += scalar_weight_j*fabs(dx_ij_minus[2]);

  si->feedback_data.f_plus_denom[0] += scalar_weight_j*fabs(dx_ij_plus[0]);
  si->feedback_data.f_plus_denom[1] += scalar_weight_j*fabs(dx_ij_plus[1]);
  si->feedback_data.f_plus_denom[2] += scalar_weight_j*fabs(dx_ij_plus[2]);


  si->feedback_data.f_minus_num[0] += scalar_weight_j*fabs(dx_ij_plus[0]);
  si->feedback_data.f_minus_num[1] += scalar_weight_j*fabs(dx_ij_plus[1]);
  si->feedback_data.f_minus_num[2] += scalar_weight_j*fabs(dx_ij_plus[2]);

  si->feedback_data.f_minus_denom[0] += scalar_weight_j*fabs(dx_ij_minus[0]);
  si->feedback_data.f_minus_denom[1] += scalar_weight_j*fabs(dx_ij_minus[1]);
  si->feedback_data.f_minus_denom[2] += scalar_weight_j*fabs(dx_ij_minus[2]);
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_feedback(const float r2, const float *dx,
                                  const float hi, const float hj,
                                  struct spart *restrict si,
                                  struct part *restrict pj, const float a,
                                  const float H) {
#ifdef DEBUG_INTERACTIONS_STARS
  /* Update ngb counters */
  if (si->num_ngb_force < MAX_NUM_OF_NEIGHBOURS_STARS)
    si->ids_ngbs_force[si->num_ngb_force] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_force;
#endif
}

#endif /* SWIFT_GEAR_STARS_IACT_H */
