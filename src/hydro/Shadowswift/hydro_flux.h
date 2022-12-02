//
// Created by yuyttenh on 29/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_FLUX_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_FLUX_H

#include "hydro_unphysical.h"
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
}

__attribute((always_inline)) INLINE static void hydro_part_apply_fluxes(
    struct part* p, const struct hydro_props* hydro_props,
    const struct cosmology* cosmo) {

  float flux[5];
  hydro_part_get_fluxes(p, flux);

  /* Update conserved variables. */
  p->conserved.mass += flux[0];
  p->conserved.momentum[0] += flux[1];
  p->conserved.momentum[1] += flux[2];
  p->conserved.momentum[2] += flux[3];
#if defined(EOS_ISOTHERMAL_GAS)
  /* We use the EoS equation in a sneaky way here just to get the constant u
   */
  p->conserved.energy =
      p->conserved.mass * gas_internal_energy_from_entropy(0.0f, 0.0f);
#else
  p->conserved.energy += flux[4];
#endif

#ifndef HYDRO_GAMMA_5_3

  const float Pcorr = (dt_hydro - dt_therm) * p->geometry.volume;
  p->conserved.momentum[0] -= Pcorr * p->gradients.P[0];
  p->conserved.momentum[1] -= Pcorr * p->gradients.P[1];
  p->conserved.momentum[2] -= Pcorr * p->gradients.P[2];
#ifdef SHADOWSWIFT_TOTAL_ENERGY
  p->conserved.energy -=
      Pcorr * (p->v[0] * p->gradients.P[0] + p->v[1] * p->gradients.P[1] +
               p->v[2] * p->gradients.P[2]);
#endif
#endif

  /* Apply the minimal energy limit */
  const float min_energy =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;
  if (p->conserved.energy < min_energy * p->conserved.mass) {
    p->conserved.energy = min_energy * p->conserved.mass;
  }

  // MATTHIEU: Apply the entropy floor here.

  /* Check conserved quantities */
  shadowswift_check_physical_quantities(
      "mass", "energy", p->conserved.mass, p->conserved.momentum[0],
      p->conserved.momentum[1], p->conserved.momentum[2], p->conserved.energy);

#ifdef SWIFT_DEBUG_CHECKS
  if (p->conserved.mass < 0.) {
    error(
        "Negative mass after conserved variables update (mass: %g, dmass: "
        "%g)!",
        p->conserved.mass, p->flux.mass);
  }

  if (p->conserved.energy < 0.) {
    error(
        "Negative energy after conserved variables update (energy: %g, "
        "denergy: %g)!",
        p->conserved.energy, p->flux.energy);
  }
#endif

  /* Reset the fluxes so that they do not get used again in the kick1. */
  hydro_part_reset_fluxes(p);
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
 * @param dt Time-step to integrate fluxes over
 * @param fluxes Array to store the result in (of size 5 or more).
 */
__attribute__((always_inline)) INLINE static void hydro_compute_flux(
    const float* WL, const float* WR, const float* n_unit, const float* vLR,
    const float Anorm, const float dt, float* fluxes) {

  riemann_solve_for_flux(WL, WR, n_unit, vLR, fluxes);

  fluxes[0] *= Anorm * dt;
  fluxes[1] *= Anorm * dt;
  fluxes[2] *= Anorm * dt;
  fluxes[3] *= Anorm * dt;
  fluxes[4] *= Anorm * dt;
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the left of the interparticle interface.
 *
 * @param p Particle.
 * @param fluxes Fluxes accross the interface.
 * @param dx Distance between the particles that share the interface.
 */
__attribute__((always_inline)) INLINE static void hydro_part_update_fluxes_left(
    struct part* restrict p, const float* fluxes, const float* dx) {

  p->gravity.mflux[0] += fluxes[0] * dx[0];
  p->gravity.mflux[1] += fluxes[0] * dx[1];
  p->gravity.mflux[2] += fluxes[0] * dx[2];

  p->flux.mass -= fluxes[0];
  p->flux.momentum[0] -= fluxes[1];
  p->flux.momentum[1] -= fluxes[2];
  p->flux.momentum[2] -= fluxes[3];
  p->flux.energy -= fluxes[4];

#ifndef SHADOWSWIFT_TOTAL_ENERGY
  const float ekin =
      0.5f * (p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]);
  p->flux.energy += fluxes[1] * p->v[0];
  p->flux.energy += fluxes[2] * p->v[1];
  p->flux.energy += fluxes[3] * p->v[2];
  p->flux.energy -= fluxes[0] * ekin;
#endif

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
 * @param fluxes Fluxes accross the interface.
 * @param dx Distance between the particles that share the interface.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_update_fluxes_right(struct part* restrict p, const float* fluxes,
                               const float* dx) {

  p->gravity.mflux[0] += fluxes[0] * dx[0];
  p->gravity.mflux[1] += fluxes[0] * dx[1];
  p->gravity.mflux[2] += fluxes[0] * dx[2];

  p->flux.mass += fluxes[0];
  p->flux.momentum[0] += fluxes[1];
  p->flux.momentum[1] += fluxes[2];
  p->flux.momentum[2] += fluxes[3];
  p->flux.energy += fluxes[4];

#ifndef SHADOWSWIFT_TOTAL_ENERGY
  const float ekin =
      0.5f * (p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]);
  p->flux.energy -= fluxes[1] * p->v[0];
  p->flux.energy -= fluxes[2] * p->v[1];
  p->flux.energy -= fluxes[3] * p->v[2];
  p->flux.energy += fluxes[0] * ekin;
#endif

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
