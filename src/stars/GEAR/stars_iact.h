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


/* For A, we need the d Kernel / d r_j two times (once with r_j, h_i and the
     second with r_j, h_j) */
  /* float r_i2 = si->x[0]*si->x[0] + si->x[1]*si->x[1] + si->x[1]*si->x[1]; */
  const float r_j2 = pj->x[0]*si->x[0] + pj->x[1]*si->x[1] + pj->x[1]*si->x[1];

  float dW_ij_dr_j;
  float dW_jj_dr_j;
  const float u_ij = sqrt(r_j2)/hi;
  const float u_jj = sqrt(r_j2)/hj;
  float dummy_W;

  kernel_deval(u_ij, &dummy_W, &dW_ij_dr_j);
  kernel_deval(u_jj, &dummy_W, &dW_jj_dr_j);

  dW_ij_dr_j *= pow_dimension(1.0/hi); /* 1/h_i^d */
  dW_jj_dr_j *= pow_dimension(1.0/hj); /* 1/h_j^d */

  /* Compute the projection vectors */
  const double dx_ij_plus[3] = {max(dx[0], 0.0)/r,
			  max(dx[1], 0.0)/r,
			  max(dx[2], 0.0)/r};

  const double dx_ij_minus[3] = {min(dx[0], 0.0)/r,
				 min(dx[1], 0.0)/r,
				 min(dx[2], 0.0)/r};

  const double dx_ij_hat[3] = {(dx_ij_plus[0] + dx_ij_minus[0])/r,
			       (dx_ij_plus[1] + dx_ij_minus[1])/r,
			       (dx_ij_plus[2] + dx_ij_minus[2])/r};

  /* I also need the x_ij_hat... which must be calculated before... */
  double n_bar_i_2_inv = 1.0/(si->density.wcount*si->density.wcount);
  double n_bar_j_2_inv = 1.0/(pj->density.wcount*pj->density.wcount);

  double A_j[3] = {(n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[0],
		   (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[1],
		   (n_bar_i_2_inv * dW_ij_dr_j + n_bar_j_2_inv * dW_jj_dr_j)*dx_ij_hat[2]};

  double number_1 = A_j[0]*dx_ij_hat[0] *  A_j[1]*dx_ij_hat[1] *  A_j[2]*dx_ij_hat[2];
  double number_2 = M_PI * r2;
  double denom = sqrt(1 + number_1/number_2);

  /* Store at least this value. I would not be against storing dx_ij and A to
     avoid computing them many times and making mistakes... */
  double scalar_weight_j = 0.5*(1.0 + 1.0/denom) ;

  si->feedback_data.f_plus_denom[0] += scalar_weight_j*fabs(dx_ij_minus[0]);
  si->feedback_data.f_plus_denom[1] += scalar_weight_j*fabs(dx_ij_minus[1]);
  si->feedback_data.f_plus_denom[2] += scalar_weight_j*fabs(dx_ij_minus[2]);

  si->feedback_data.f_minus_denom[0] += scalar_weight_j*fabs(dx_ij_plus[0]);
  si->feedback_data.f_minus_denom[1] += scalar_weight_j*fabs(dx_ij_plus[1]);
  si->feedback_data.f_minus_denom[2] += scalar_weight_j*fabs(dx_ij_plus[2]);
}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_prep1(const float r2, const float dx[3],
                               const float hi, const float hj, struct spart *si,
                               const struct part *pj, const float a,
                               const float H) {}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_prep2(const float r2, const float dx[3],
                               const float hi, const float hj, struct spart *si,
                               const struct part *pj, const float a,
                               const float H) {}

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
