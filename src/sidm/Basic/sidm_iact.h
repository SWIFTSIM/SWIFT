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

/* Local headers. */
#include "sidm_properties.h"
#include "timeline.h"

/**
 * @brief Kernel overlap assuming top-hat kernel. Equation 2 of Fisher et al
 * 2026.
 *
 * @param r Comoving distance between the two particles.
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 */
INLINE static float sidm_kernel_overlap_tophat(const float r, const float hi,
                                               const float hj) {

  const float Vconst = (4.f / 3.f) * M_PI;
  const float Vi = Vconst * hi * hi * hi;
  const float Vj = Vconst * hj * hj * hj;

  float V_overlap;

  if (r >= hi + hj) {
    return 0.f;
  }

  // sphere within sphere case
  if (r <= fabsf(hi - hj)) {
    const float hmin = hi < hj ? hi : hj;
    V_overlap = Vconst * hmin * hmin * hmin;
  } else {
    // usual case
    const float term = hi + hj - r;
    V_overlap = M_PI * term * term *
                (r * r + 2 * r * (hi + hj) - 3 * (hi - hj) * (hi - hj)) /
                (12.f * r);
  }

  return V_overlap / (Vi * Vj);
}

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
    const float H, const int with_cosmology, const struct cosmology *cosmo,
    const struct sidm_props *sidm_props, const integertime_t ti_current,
    const double time_base) {

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
runner_iact_nonsym_sidm_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct sipart *restrict sipi, const struct sipart *restrict sipj,
    const float a, const float H, const int with_cosmology,
    const struct cosmology *cosmo, const struct sidm_props *sidm_props,
    const integertime_t ti_current, const double time_base) {

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

/**
 * @brief Force interaction between two particles (symmetric).
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
__attribute__((always_inline)) INLINE static void runner_iact_sidm_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct sipart *restrict sipi, struct sipart *restrict sipj, const float a,
    const float H, const int with_cosmology, const struct cosmology *cosmo,
    const struct sidm_props *sidm_props, const integertime_t ti_current,
    const double time_base) {

  message("running sym SIDM force");

  /* Pair separation. */
  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mi = sipi->mass;
  const float mj = sipj->mass;

  /* Compute pair relative velocity. */
  double dv[3];
  dv[0] = sipi->v[0] - sipj->v[0];
  dv[1] = sipi->v[1] - sipj->v[1];
  dv[2] = sipi->v[2] - sipj->v[2];
  const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
  const double vij = sqrt(v2) * cosmo->a_inv;

  /* Get time-step for sipi and sipj */
  double dt_sipi, dt_sipj;

  if (with_cosmology) {

    const integertime_t ti_begin =
        get_integer_time_begin(ti_current - 1, sipi->time_bin);
    const integertime_t ti_step = get_integer_timestep(sipi->time_bin);
    dt_sipi = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);

    const integertime_t tj_begin =
        get_integer_time_begin(ti_current - 1, sipj->time_bin);
    const integertime_t tj_step = get_integer_timestep(sipj->time_bin);
    dt_sipj = cosmology_get_delta_time(cosmo, tj_begin, tj_begin + tj_step);

  } else {
    dt_sipi = get_timestep(sipi->time_bin, time_base);
    dt_sipj = get_timestep(sipj->time_bin, time_base);
  }

  float lambda_ij = sidm_kernel_overlap_tophat(r, hi, hj) * cosmo->a3_inv;

  /* Scattering rates */
  float SIDM_rate_i = mj * sidm_props->sigma_over_m * vij * lambda_ij;
  float SIDM_rate_j = mi * sidm_props->sigma_over_m * vij * lambda_ij;

  sipi->SIDM_rate += SIDM_rate_i;
  sipj->SIDM_rate += SIDM_rate_j;
}

/**
 * @brief Force interaction between two particles (non-symmetric).
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
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_sidm_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct sipart *restrict sipi, struct sipart *restrict sipj, const float a,
    const float H, const int with_cosmology, const struct cosmology *cosmo,
    const struct sidm_props *sidm_props, const integertime_t ti_current,
    const double time_base) {

  /* Pair separation. */
  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mj = sipj->mass;

  /* Compute pair relative velocity. */
  double dv[3];
  dv[0] = sipi->v[0] - sipj->v[0];
  dv[1] = sipi->v[1] - sipj->v[1];
  dv[2] = sipi->v[2] - sipj->v[2];
  const double v2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
  const double vij = sqrt(v2) * cosmo->a_inv;

  /* Get time-step for sipi and sipj */
  double dt_sipi;

  if (with_cosmology) {

    const integertime_t ti_begin =
        get_integer_time_begin(ti_current - 1, sipi->time_bin);
    const integertime_t ti_step = get_integer_timestep(sipi->time_bin);
    dt_sipi = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);

  } else {
    dt_sipi = get_timestep(sipi->time_bin, time_base);
  }

  float lambda_ij = sidm_kernel_overlap_tophat(r, hi, hj) * cosmo->a3_inv;

  /* Scattering rates */
  float SIDM_rate_i = mj * sidm_props->sigma_over_m * vij * lambda_ij;
  sipi->SIDM_rate += SIDM_rate_i;
}

#endif /* SWIFT_BASIC_SIDM_IACT_H */