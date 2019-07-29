/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2019 Fabien Jeanquartier (fabien.jeanquartier@epfl.ch)
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
#ifndef SWIFT_GEAR_STAR_FORMATION_IACT_H
#define SWIFT_GEAR_STAR_FORMATION_IACT_H

/**
 * @file GEAR/star_formation_iact.h
 * @brief Density computation
 */

/**
 * @brief do star_formation computation after the runner_iact_density (symmetric
 * version)
 *
 * Compute the velocity dispersion follow eq. 2 in Revaz & Jablonka 2018.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_star_formation(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float wi;
  float wj;
  /* Evaluation of the SPH kernel */
  kernel_eval(sqrt(r2) / hi, &wi);
  kernel_eval(sqrt(r2) / hj, &wj);

  /* Delta v */
  float dv[3] = {
    pi->v[0] - pj->v[0],
    pi->v[1] - pj->v[1],
    pi->v[2] - pj->v[2]
  };

  /* Norms at power 2 */
  const float norm_v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
  const float norm_x2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  /* Compute the velocity dispersion */
  const float sigma2 = norm_v2 + H * norm_x2;

  /* Compute the velocity dispersion */
  pi->sf_data.sigma2 += sigma2 * wi * hydro_get_mass(pj);
  pj->sf_data.sigma2 += sigma2 * wj * hydro_get_mass(pi);
}

/**
 * @brief do star_formation computation after the runner_iact_density (non
 * symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_star_formation(float r2, const float *dx, float hi, float hj,
                                  struct part *restrict pi,
                                  const struct part *restrict pj, float a,
                                  float H) {
  float wi;
  /* Evaluation of the SPH kernel */
  kernel_eval(sqrt(r2) / hi, &wi);

  /* Delta v */
  float dv[3] = {
    pi->v[0] - pj->v[0],
    pi->v[1] - pj->v[1],
    pi->v[2] - pj->v[2]
  };

  /* Norms at power 2 */
  const float norm_v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
  const float norm_x2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  /* Compute the velocity dispersion */
  const float sigma2 = norm_v2 + H * norm_x2;

  /* Compute the velocity dispersion */
  pi->sf_data.sigma2 += sigma2 * wi * hydro_get_mass(pj);
}

#endif /* SWIFT_GEAR_STAR_FORMATION_IACT_H */
