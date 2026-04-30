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
#include "random.h"
#include "sidm_kernel.h"
#include "sidm_properties.h"
#include "timeline.h"

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
 * @brief Apply SIDM velocity kick to particles (elastic, isotropic scattering).
 *
 * @param sipi First part*icle.
 * @param sipj Second part*icle.
 * @param a Current scale factor.
 * @param vij Relative speed.
 */
__attribute__((always_inline)) INLINE static void kick_siparts(
    struct sipart *restrict sipi, struct sipart *restrict sipj, const float a,
    const double vij, const integertime_t ti_current) {

  /* Get COM velocity. */
  const float mi = sipi->mass;
  const float mj = sipj->mass;
  const float mi_plus_mj = mi + mj;

  const float v_com[3] = {(mi * sipi->v[0] + mj * sipj->v[0]) / mi_plus_mj,
                          (mi * sipi->v[1] + mj * sipj->v[1]) / mi_plus_mj,
                          (mi * sipi->v[2] + mj * sipj->v[2]) / mi_plus_mj};

  /* Get scattering direction. */
  /* For isotropic sampling:
   * Draw cos(theta) uniformly in [-1, 1] and phi uniformly in [0, 2*pi). */
  const float cos_theta =
      2.f * random_unit_interval_two_IDs(sipi->id, sipj->id, ti_current,
                                         random_number_sidm_polar_angle) -
      1.f;
  const float sin_theta = sqrtf(max(0.f, 1.f - cos_theta * cos_theta));
  const float phi =
      2.f * M_PI *
      random_unit_interval_two_IDs(sipi->id, sipj->id, ti_current,
                                   random_number_sidm_azimuthal_angle);

  /* New velocity direction, magnitude of relative veolcity is conserved. */
  const float v_new[3] = {vij * sin_theta * cosf(phi),
                          vij * sin_theta * sinf(phi), vij * cos_theta};

  /* Transform back to original reference frame. */
  /* Mass fractions that distribute new velocity kick. */
  const float fi = mj / mi_plus_mj;
  const float fj = mi / mi_plus_mj;

  /* By momentum conservation sipj gets the opposite COM-frame kick. */
  sipi->v[0] = v_com[0] + fi * v_new[0];
  sipi->v[1] = v_com[1] + fi * v_new[1];
  sipi->v[2] = v_com[2] + fi * v_new[2];

  sipj->v[0] = v_com[0] - fj * v_new[0];
  sipj->v[1] = v_com[1] - fj * v_new[1];
  sipj->v[2] = v_com[2] - fj * v_new[2];
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

  const double hubble_flow = cosmo->a * cosmo->H * r;
  const double vij = sqrt(v2) * cosmo->a_inv + hubble_flow;

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

  /* Interaction probability for sipi and sipj */
  float pij = SIDM_rate_i * dt_sipi;
  float pji = SIDM_rate_j * dt_sipj;
  float p = (pij + pji) * 0.5;

  /* Pair interaction? */  // need two ids here
  const float x = random_unit_interval_two_IDs(sipi->id, sipj->id, ti_current,
                                               random_number_sidm_scattering);

  if (x < p) {
    kick_siparts(sipi, sipj, a, vij, ti_current);
  }
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

  /* Interaction probability */
  float pij = SIDM_rate_i * dt_sipi;

  /* Pair interaction? */
  const float x = random_unit_interval_two_IDs(sipi->id, sipj->id, ti_current,
                                               random_number_sidm_scattering);

  if (x < pij) {
    kick_siparts(sipi, sipj, a, vij, ti_current);
  }
}

#endif /* SWIFT_BASIC_SIDM_IACT_H */
