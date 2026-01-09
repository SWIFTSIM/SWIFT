/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Katy Proctor (katy.proctor@fysik.su.se)
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
#ifndef SWIFT_BASIC_SIDM_IACT_H
#define SWIFT_BASIC_SIDM_IACT_H

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (sipi - sipj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param sipi First part*icle.
 * @param sipj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_sidm_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct sipart *restrict sipi, struct sipart *restrict sipj, const float a,
    const float H) {

  float wi, wj, wi_dx, wj_dx;

  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mi = sipi->mass;
  const float mj = sipj->mass;

  /* Compute density of sipi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);

  sipi->rho += mj * wi;
  sipi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  sipi->density.wcount += wi;
  sipi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  sipj->rho += mi * wj;
  sipj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
  sipj->density.wcount += wj;
  sipj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

#ifdef SWIFT_SIDM_DENSITY_CHECKS
  sipi->n += wi;
  sipj->n += wj;
  sipi->N_density++;
  sipj->N_density++;
#endif
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sidm_density(const float r2, const float dx[3],
                                const float hi, const float hj,
                                struct sipart *restrict sipi,
                                const struct sipart *restrict sipj,
                                const float a, const float H) {

  float wi, wi_dx;

  /* Get the masses. */
  const float mj = sipj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  sipi->rho += mj * wi;
  sipi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  sipi->density.wcount += wi;
  sipi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

#ifdef SWIFT_SIDM_DENSITY_CHECKS
  sipi->n += wi;
  sipi->N_density++;
#endif
}

#endif /* SWIFT_BASIC_SIDM_IACT_H */