//
// Created by yuyttenh on 29/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_FLUX_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_FLUX_H

#include "riemann.h"

/**
 * @brief Reset the fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_part_reset_fluxes(
    struct part* restrict p) {

  p->flux.mass = 0.0f;
  p->flux.momentum[0] = 0.0f;
  p->flux.momentum[1] = 0.0f;
  p->flux.momentum[2] = 0.0f;
  p->flux.energy = 0.0f;
  p->flux.entropy = 0.0f;
  p->flux.dt = -1.f;

  p->gravity.mflux[0] = 0.0f;
  p->gravity.mflux[1] = 0.0f;
  p->gravity.mflux[2] = 0.0f;
}

/**
 * @brief Get the fluxes for the given particle.
 *
 * @param p Particle.
 * @param flux Fluxes for the particle (array of size 5 or more).
 */
__attribute__((always_inline)) INLINE static void hydro_part_get_fluxes(
    const struct part* restrict p, float* flux) {

  flux[0] = p->flux.mass;
  flux[1] = p->flux.momentum[0];
  flux[2] = p->flux.momentum[1];
  flux[3] = p->flux.momentum[2];
  flux[4] = p->flux.energy;
  flux[5] = p->flux.entropy;
}

/**
 * @brief Compute the flux for the Riemann problem with the given left and right
 * state, and interface normal, surface area and velocity.
 *
 * @param WL Left state variables.
 * @param WR Right state variables.
 * @param n_unit Unit vector of the interface.
 * @param vLR Velocity of the interface.
 * @param Anorm Surface area of the interface.
 * @param fluxes Array to store the result in (of size 5 or more).
 */
__attribute__((always_inline)) INLINE static void hydro_compute_flux(
    const float* WL, const float* WR, const float* n_unit, const float* vLR,
    const float Anorm, float* fluxes) {

  riemann_solve_for_flux(WL, WR, n_unit, vLR, fluxes);
  float entropy_flux;
  if (fluxes[0] > 0.f) {
    entropy_flux = fluxes[0] * WL[5];
  } else {
    entropy_flux = fluxes[0] * WR[5];
  }

  fluxes[0] *= Anorm;
  fluxes[1] *= Anorm;
  fluxes[2] *= Anorm;
  fluxes[3] *= Anorm;
  fluxes[4] *= Anorm;
  fluxes[5] = Anorm * entropy_flux;
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the left of the interparticle interface.
 *
 * @param p Particle.
 * @param fluxes Time integrated fluxes across the interface.
 * @param dx Vector pointing from right particle to left particle.
 */
__attribute__((always_inline)) INLINE static void hydro_part_update_fluxes_left(
    struct part* restrict p, const float* restrict fluxes, const float* restrict dx) {

  p->gravity.mflux[0] += fluxes[0] * dx[0];
  p->gravity.mflux[1] += fluxes[0] * dx[1];
  p->gravity.mflux[2] += fluxes[0] * dx[2];

  p->flux.mass -= fluxes[0];
  p->flux.momentum[0] -= fluxes[1];
  p->flux.momentum[1] -= fluxes[2];
  p->flux.momentum[2] -= fluxes[3];
  p->flux.energy -= fluxes[4];
  p->flux.entropy -= fluxes[5];

  if (dx[0] < 0) {
    p->flux_count -= 1;
  } else {
    p->flux_count += 1;
  }
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the right of the interparticle interface.
 *
 * @param p Particle.
 * @param fluxes Time integrated fluxes across the interface.
 * @param dx Vector pointing from right particle to left particle.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_update_fluxes_right(struct part* restrict p, const float* restrict fluxes,
                               const float* restrict dx) {

  p->gravity.mflux[0] += fluxes[0] * dx[0];
  p->gravity.mflux[1] += fluxes[0] * dx[1];
  p->gravity.mflux[2] += fluxes[0] * dx[2];

  p->flux.mass += fluxes[0];
  p->flux.momentum[0] += fluxes[1];
  p->flux.momentum[1] += fluxes[2];
  p->flux.momentum[2] += fluxes[3];
  p->flux.energy += fluxes[4];
  p->flux.entropy += fluxes[5];

  if (dx[0] < 0) {
    p->flux_count += 1;
  } else {
    p->flux_count -= 1;
  }
}

/**
 * @brief Get the drift term for the density based on the mass flux.
 *
 * @param mass_flux Current mass flux for the particle.
 * @param dt Drift time step (in co-moving units).
 * @param volume Volume of the particle.
 * @return Term that will drift the density.
 */
__attribute__((always_inline)) INLINE static float
hydro_flux_density_drift_term(const float mass_flux, const float dt,
                              const float volume) {

  if (volume > 0.0f) {
    return mass_flux * dt / volume;
  } else {
    return 0.0f;
  }
}

#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_FLUX_H
