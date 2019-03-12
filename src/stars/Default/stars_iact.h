/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_DEFAULT_STARS_IACT_H
#define SWIFT_DEFAULT_STARS_IACT_H

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
runner_iact_nonsym_stars_density(float r2, const float *dx, float hi, float hj,
                                 struct spart *restrict si,
                                 const struct part *restrict pj, float a,
                                 float H) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

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
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
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
runner_iact_nonsym_stars_feedback(float r2, const float *dx, float hi, float hj,
                                  struct spart *restrict si,
                                  struct part *restrict pj, float a, float H) {

  const float mj = hydro_get_mass(pj);
  const float rhoj = hydro_get_comoving_density(pj);
  const float r = sqrtf(r2);
  const float ri = 1.f / r;

  /* Get the kernel for hi. */
  float hi_inv = 1.0f / hi;
  float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  float wi_dr = hid_inv * wi_dx;

  /* Compute dv dot r */
  float dvdr = (si->v[0] - pj->v[0]) * dx[0] + (si->v[1] - pj->v[1]) * dx[1] +
               (si->v[2] - pj->v[2]) * dx[2];

  /* Get the time derivative for h. */
  si->feedback.h_dt -= mj * dvdr * ri / rhoj * wi_dr;

#ifdef DEBUG_INTERACTIONS_STARS
  /* Update ngb counters */
  if (si->num_ngb_force < MAX_NUM_OF_NEIGHBOURS_STARS)
    si->ids_ngbs_force[si->num_ngb_force] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_force;
#endif
}

#endif /* SWIFT_DEFAULT_STARS_IACT_H */
