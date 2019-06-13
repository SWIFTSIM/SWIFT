/*******************************************************************************
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
    struct part *restrict pj, struct xpart *restrict xpi,
    struct xpart *restrict xpj, float a, float H) {
  /* The goal here is to estimate the local turbulence */
  /* Value used to evaluate the SPH kernel */
  float wi;
  float wj;
  /* Evaluation of the SPH kernel */
  kernel_eval(sqrt(r2) / hi, &wi);
  kernel_eval(sqrt(r2) / hj, &wj);
  /* Square of the velocity norm between particles i and j */
  float norm_v2 = pow(pi->v[0] - pj->v[0], 2) + pow(pi->v[1] - pj->v[1], 2) +
                  pow(pi->v[2] - pj->v[2], 2);

  /* Estimation of local turbulence for pi and pj, see Revaz & Jablonka, eq (2)
   */
  xpi->sf_data.sigma2 += pow(hi, -3) * norm_v2 * wi * hydro_get_mass(pj);
  xpj->sf_data.sigma2 += pow(hj, -3) * norm_v2 * wj * hydro_get_mass(pi);
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
                                  const struct part *restrict pj,
                                  struct xpart *restrict xpi,
                                  const struct xpart *restrict xpj, float a,
                                  float H)
/* The goal here is to estimate the local turbulence */
/* Value used to evaluate the SPH kernel */
{
  float wi;
  /* Evaluation of the SPH kernel */
  kernel_eval(sqrt(r2) / hi, &wi);
  /* Square of the velocity norm */
  float norm_v2 = pow(pi->v[0] - pj->v[0], 2) + pow(pi->v[1] - pj->v[1], 2) +
                  pow(pi->v[2] - pj->v[2], 2);
  /* Estimation of local turbulence for pi */
  xpi->sf_data.sigma2 += pow(hi, -3) * norm_v2 * wi * hydro_get_mass(pj);
}
#endif /* SWIFT_GEAR_STAR_FORMATION_IACT_H */
