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
#include "kernel_dark_matter.h"
#include "random.h"
#include "sidm_debug.h"
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
 * @brief Returns a random number in [0, 1[ for an SIDM pair.
 *
 * The IDs are ordered before being combined so that the number does not
 * depend on which particle of the pair plays the role of i (e.g. when the
 * pair is evaluated from both sides over MPI).
 *
 * @param id_i ID of the first particle.
 * @param id_j ID of the second particle.
 * @param ti_current The current integer time.
 * @param type The #random_number_type to generate.
 */
__attribute__((always_inline)) INLINE static double sidm_random_unit_interval(
    const long long id_i, const long long id_j, const integertime_t ti_current,
    const enum random_number_type type) {

  const long long id_lo = min(id_i, id_j);
  const long long id_hi = max(id_i, id_j);

  return random_unit_interval_two_IDs(id_lo, id_hi, ti_current, type);
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

  const float v_com[3] = {
      (mi * sipi->v[0] + mj * sipj->v[0]) / (mi_plus_mj * a),
      (mi * sipi->v[1] + mj * sipj->v[1]) / (mi_plus_mj * a),
      (mi * sipi->v[2] + mj * sipj->v[2]) / (mi_plus_mj * a)};

  /* Get scattering direction. */
  /* For isotropic sampling:
   * Draw cos(theta) uniformly in [-1, 1] and phi uniformly in [0, 2*pi). */
  const float cos_theta =
      2.f * sidm_random_unit_interval(sipi->id, sipj->id, ti_current,
                                      random_number_sidm_polar_angle) -
      1.f;
  const float sin_theta = sqrtf(max(0.f, 1.f - cos_theta * cos_theta));
  const float phi =
      2.f * M_PI *
      sidm_random_unit_interval(sipi->id, sipj->id, ti_current,
                                random_number_sidm_azimuthal_angle);

  /* New velocity direction, magnitude of relative veolcity is conserved. */
  const float v_new[3] = {vij * sin_theta * cosf(phi),
                          vij * sin_theta * sinf(phi), vij * cos_theta};

  /* Transform back to original reference frame. */
  /* Mass fractions that distribute new velocity kick. */
  const float fi = mj / mi_plus_mj;
  const float fj = mi / mi_plus_mj;

  /* The particle with the lower ID always gets +v_new. */
  const float sign = (sipi->id < sipj->id) ? 1.f : -1.f;

  /* By momentum conservation sipj gets the opposite COM-frame kick.
  Also convert velocities back to comoving. */
  sipi->v[0] = a * (v_com[0] + sign * fi * v_new[0]);
  sipi->v[1] = a * (v_com[1] + sign * fi * v_new[1]);
  sipi->v[2] = a * (v_com[2] + sign * fi * v_new[2]);

  sipj->v[0] = a * (v_com[0] - sign * fj * v_new[0]);
  sipj->v[1] = a * (v_com[1] - sign * fj * v_new[1]);
  sipj->v[2] = a * (v_com[2] - sign * fj * v_new[2]);
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

  /* Pair scattering rate. The mean mass makes it symmetric under i <-> j
   * and matches the mean momentum exchange between the two particles. */
  const float m_pair = 0.5f * (mi + mj);
  const float pair_rate = m_pair * sidm_props->sigma_over_m * vij * lambda_ij;

  sipi->SIDM_rate += pair_rate;
  sipj->SIDM_rate += pair_rate;

  /* Inertia tensor contributions */
  sipi->inertia_tensor[0] += mj * (r2 - dx[0] * dx[0]);
  sipi->inertia_tensor[1] += mj * (r2 - dx[1] * dx[1]);
  sipi->inertia_tensor[2] += mj * (r2 - dx[2] * dx[2]);
  sipi->inertia_tensor[3] -= mj * dx[0] * dx[1];
  sipi->inertia_tensor[4] -= mj * dx[0] * dx[2];
  sipi->inertia_tensor[5] -= mj * dx[1] * dx[2];

  sipj->inertia_tensor[0] += mi * (r2 - dx[0] * dx[0]);
  sipj->inertia_tensor[1] += mi * (r2 - dx[1] * dx[1]);
  sipj->inertia_tensor[2] += mi * (r2 - dx[2] * dx[2]);
  sipj->inertia_tensor[3] -= mi * dx[0] * dx[1];
  sipj->inertia_tensor[4] -= mi * dx[0] * dx[2];
  sipj->inertia_tensor[5] -= mi * dx[1] * dx[2];

  /* Interaction probability. Use the shorter time-step: the earlier part of
   * the longer step was already covered by the non-symmetric interactions
   * made while only the particle on the shorter step was active. */
  const float p = pair_rate * min(dt_sipi, dt_sipj);

  /* Pair interaction? */
  const float x = sidm_random_unit_interval(sipi->id, sipj->id, ti_current,
                                            random_number_sidm_scattering);

  if (x < p) {
    sidm_debug_interaction(sipi, sipj);
    kick_siparts(sipi, sipj, a, vij, ti_current);
    sidm_debug_interaction(sipi, sipj);
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

  /* Pair scattering rate (same definition as in the symmetric case) */
  const float m_pair = 0.5f * (mi + mj);
  const float pair_rate = m_pair * sidm_props->sigma_over_m * vij * lambda_ij;
  sipi->SIDM_rate += pair_rate;

  /* Inertia tensor contribution */
  sipi->inertia_tensor[0] += mj * (r2 - dx[0] * dx[0]);
  sipi->inertia_tensor[1] += mj * (r2 - dx[1] * dx[1]);
  sipi->inertia_tensor[2] += mj * (r2 - dx[2] * dx[2]);
  sipi->inertia_tensor[3] -= mj * dx[0] * dx[1];
  sipi->inertia_tensor[4] -= mj * dx[0] * dx[2];
  sipi->inertia_tensor[5] -= mj * dx[1] * dx[2];

  /* Interaction probability, as in the symmetric case. sipj is usually
   * inactive (so dt_sipi is the shorter step), but it can also be active and
   * excluded here because it is foreign or outside this task's h range. */
  const float p = pair_rate * min(dt_sipi, dt_sipj);

  /* Pair interaction? */
  const float x = sidm_random_unit_interval(sipi->id, sipj->id, ti_current,
                                            random_number_sidm_scattering);

  if (x < p) {
    sidm_debug_interaction(sipi, sipj);
    kick_siparts(sipi, sipj, a, vij, ti_current);
    sidm_debug_interaction(sipi, sipj);
  }
}

#endif /* SWIFT_BASIC_SIDM_IACT_H */