/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IACT_GEAR_H
#define SWIFT_RT_IACT_GEAR_H

#include "rt_gradients.h"

/**
 * @file src/rt/GEAR/rt_iact.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * particle interactions.
 */

/**
 * @brief Injection step interaction between star and hydro particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si Star particle.
 * @param pj Hydro particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_inject(
    const float r2, float *dx, const float hi, const float hj,
    struct spart *restrict si, struct part *restrict pj, float a, float H) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  si->rt_data.debug_iact_hydro_inject += 1;
  si->rt_data.debug_radiation_emitted_tot += 1ULL;

  pj->rt_data.debug_iact_stars_inject += 1;
  pj->rt_data.debug_radiation_absorbed_tot += 1ULL;
#endif

  /* If the star doesn't have any neighbours, we
   * have nothing to do here. */
  if (si->density.wcount == 0.f) return;

  float wi;
  const float r = sqrtf(r2);
  const float hi_inv = 1.f / hi;
  const float xi = r * hi_inv;
  kernel_eval(xi, &wi);
  const float hi_inv_dim = pow_dimension(hi_inv);
  /* psi(x_star - x_gas, h_star) */
  const float psi = wi * hi_inv_dim / si->density.wcount;

  /* TODO: this is done differently for RT_HYDRO_CONTROLLED_INJECTION */
  for (int g = 0; g < RT_NGROUPS; g++) {
    pj->rt_data.conserved[g].energy += si->rt_data.emission_this_step[g] * psi;
  }

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* Take note how much energy we actually injected */
  for (int g = 0; g < RT_NGROUPS; g++) {
    float res = si->rt_data.emission_this_step[g] * psi;
    si->rt_data.debug_injected_energy[g] += res;
    si->rt_data.debug_injected_energy_tot[g] += res;
  }
#endif
}

/**
 * @brief Flux calculation between particle i and particle j
 *
 * This method calls runner_iact_rt_fluxes_common with mode 1.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param mode 0 if non-symmetric interaction, 1 if symmetric
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_flux_common(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H, int mode) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (pi->rt_data.debug_injection_done != 1)
    error(
        "Trying to do iact transport when "
        "finalise injection count is %d",
        pi->rt_data.debug_injection_done);

  if (pi->rt_data.debug_calls_iact_gradient == 0)
    error(
        "Called iact transport on particle "
        "with iact gradient count 0");

  if (pi->rt_data.debug_gradients_done != 1)
    error(
        "Trying to do iact transport when "
        "rt_finalise_gradient count is %d",
        pi->rt_data.debug_gradients_done);

  pi->rt_data.debug_calls_iact_transport_interaction += 1;

  if (mode == 1) {

    if (pj->rt_data.debug_injection_done != 1)
      error(
          "Trying to do iact transport when "
          "finalise injection count is %d",
          pj->rt_data.debug_injection_done);

    if (pj->rt_data.debug_calls_iact_gradient == 0)
      error(
          "Called iact transport on particle "
          "with iact gradient count 0");

    if (pj->rt_data.debug_gradients_done != 1)
      error(
          "Trying to do iact transport when "
          "rt_finalise_gradient count is %d",
          pj->rt_data.debug_gradients_done);

    pj->rt_data.debug_calls_iact_transport_interaction += 1;
  }
#endif
}

/**
 * @brief Flux calculation between particle i and particle j
 *
 * This method calls runner_iact_rt_fluxes_common with mode 1.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_transport(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  runner_iact_rt_flux_common(r2, dx, hi, hj, pi, pj, a, H, 1);
}

/**
 * @brief Flux calculation between particle i and particle j: non-symmetric
 * version
 *
 * This method calls runner_iact_rt_fluxes_common with mode 0.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_transport(float r2, const float *dx, float hi, float hj,
                                struct part *restrict pi,
                                struct part *restrict pj, float a, float H) {

  runner_iact_rt_flux_common(r2, dx, hi, hj, pi, pj, a, H, 0);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around rt_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_gradient(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  rt_gradients_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around rt_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_gradient(float r2, const float *dx, float hi, float hj,
                               struct part *restrict pi,
                               struct part *restrict pj, float a, float H) {

  rt_gradients_nonsym_collect(r2, dx, hi, hj, pi, pj);
}

#endif /* SWIFT_RT_IACT_GEAR_H */
