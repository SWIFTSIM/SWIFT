/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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
#ifndef SWIFT_SHADOWSWIFT_HYDRO_H
#define SWIFT_SHADOWSWIFT_HYDRO_H

#include "adiabatic_index.h"
#include "approx_math.h"
#include "cosmology.h"
#include "dimension.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "hydro_flux.h"
#include "hydro_getters.h"
#include "hydro_gradients.h"
#include "hydro_gravity.h"
#include "hydro_parameters.h"
#include "hydro_properties.h"
#include "hydro_setters.h"
#include "hydro_slope_limiters.h"
#include "hydro_space.h"
#include "hydro_unphysical.h"
#include "hydro_velocities.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "pressure_floor.h"
#include "space.h"

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * This function returns the time-step of a particle given its hydro-dynamical
 * state. A typical time-step calculation would be the use of the CFL condition.
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 * @param hydro_properties The SPH parameters
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct hydro_props *restrict hydro_properties,
    const struct cosmology *restrict cosmo) {

  const float CFL_condition = hydro_properties->CFL_condition;

  /* skip the time step calculation if we are using Lloyd's algorithm */
  /* TODO */

#ifdef SWIFT_BOUNDARY_PARTICLES
  if (p->id < space_boundary_parts_interior) {
    return FLT_MAX;
  }
#endif

  float W[6];
  hydro_part_get_primitive_variables(p, W);

  /* v_full is the actual velocity of the particle, v is its
     hydrodynamical velocity. The time step depends on the relative difference
     of the two. */
  float v_rel[3];
  hydro_part_get_relative_fluid_velocity(p, v_rel);
  float vmax =
      sqrtf(v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] + v_rel[2] * v_rel[2]) +
      hydro_get_comoving_soundspeed(p);
  vmax = max(vmax, p->timestepvars.vmax);

  /* Get the comoving psize, since we will compare with another comoving
   * geometric property below */
  float psize = hydro_get_comoving_psize(p);
  /* If the particle shows large deviations from a sphere, better use the
   * minimal distance to any of its faces to compute the timestep */
  if (p->geometry.min_face_dist < 0.25 * psize &&
      p->geometry.min_face_dist > 0.) {
    psize = p->geometry.min_face_dist;
  }

  /* NOTE (yuyttenh, 06/25): To compute the (physical) dt we want to divide the
   * physical particle size (a * psize) by the physical/peculiar velocity
   * (v / a). Hence, the two a factors in front. */
  float dt = cosmo->a * cosmo->a * psize / (vmax + FLT_MIN);

#ifdef SWIFT_DEBUG_CHECKS
  if (dt == 0.f) error("Hydro part wants dt=0!");
#endif

  return CFL_condition * dt;
}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_first_init_part(
    struct part *restrict p, struct xpart *restrict xp) {

  /* Convert internal energy to thermal energy */
#ifdef EOS_ISOTHERMAL_GAS
  p->thermal_energy = Q[0] * gas_internal_energy_from_entropy(0.0f, 0.0f);
#else
  /* This needs to be corrected by a scale factor for cosmological runs, but
   * that happens in `hydro_convert_quantities(...)` */
  p->thermal_energy = p->thermal_energy * p->conserved.mass;
#endif

  /* Compute momentum from velocity and mass */
  p->conserved.momentum[0] = p->conserved.mass * p->v[0];
  p->conserved.momentum[1] = p->conserved.mass * p->v[1];
  p->conserved.momentum[2] = p->conserved.mass * p->v[2];

  /* overwrite all hydro variables if we are using Lloyd's algorithm */
  /* TODO */

  /* initialize the generator velocity based on the primitive fluid velocity */
  hydro_velocities_init(p, xp);

  /* Ignore hydro accelerations present in the initial condition */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* Set some other quantities related to geometry and time-integration to valid
   * initial values*/
  p->time_bin = 0;
  p->flux_count = 0;
  p->geometry.delaunay_flags = 0;
  p->geometry.search_radius = p->h;
  hydro_reset_timestep_vars(p);
}

/**
 * @brief Does some extra hydro operations once the actual physical time step
 * for the particle is known.
 *
 * Not used right now.
 *
 * @param p The particle to act upon.
 * @param dt Physical time step of the particle during the next step.
 */
__attribute__((always_inline)) INLINE static void hydro_timestep_extra(
    struct part *p, float dt) {}

/**
 * @brief Prepares a particle for the next hydro loop.
 *
 * We just unset the delaunay_flags
 *
 * @param p The particle to act upon
 * @param hs #hydro_space containing hydro specific space information.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part *restrict p, const struct hydro_space *hs) {
  p->geometry.delaunay_flags = 0;
}

/**
 * @brief Finishes the density calculation.
 *
 * No Density calculation for ShadowSWIFT.
 *
 * @param p The particle to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part *restrict p, const struct cosmology *cosmo) {}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * This cannot occur in the ShadowSWIFT scheme
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_part_has_no_neighbours(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo) {}

/**
 * @brief Prepare a particle for the gradient calculation.
 *
 * This function is called after the density loop and before the gradient loop.
 *
 * We use it to set the physical timestep for the particle and to copy the
 * actual velocities, which we need to boost our interfaces during the flux
 * calculation. We also initialize the variables used for the time step
 * calculation.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 * @param hydro_props Hydrodynamic properties.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_gradient(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

  hydro_gradients_init(p);
}

/**
 * @brief Resets the variables that are required for a gradient calculation.
 *
 * This function is called after hydro_prepare_gradient.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_gradient(
    struct part *restrict p) {}

/**
 * @brief Finishes the gradient calculation and prepares the particle for the
 * slope limiting loop.
 *
 * Just a wrapper around hydro_gradients_finalize, which can be an empty method,
 * in which case no gradients are used.
 *
 * This method also initializes the force loop variables.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_end_gradient(
    struct part *p) {

  hydro_gradients_finalize(p);

  /* reset the gradients if we are using Lloyd's algorith; we don't use them */
  /* TODO */

#ifdef SHADOWSWIFT_SLOPE_LIMITER_CELL_WIDE
  /* Prepare the slope limiter for this particle */
  hydro_slope_limiter_prepare(p);
#endif
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * For ShadowSWIFT, this function is called before the flux exchange, and is
 * mainly used to apply the pressure floor if needed.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 * @param hydro_props Hydrodynamic properties.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor, const float dt_alpha,
    const float dt_therm) {

  /* Apply the pressure floor. */
  /* Temporarily apply the time extrapolation to density as it is used in the
   * pressure floor... */
  const float rho = p->rho;
  p->rho += p->dW_time[0];
  const float pressure_with_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor,
                                           p->P + p->dW_time[4], cosmo) -
      p->dW_time[4];
  p->rho = rho;
  if (p->P < pressure_with_floor) {
    p->P = pressure_with_floor;
    /* TODO: Should we update the internal energy and/or entropy as well to
     * reflect this change? */
    /* TODO: Should we update the total energy? */
  }
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking  place in the various force tasks.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_reset_acceleration(
    struct part *restrict p) {

  /* Reset the acceleration. */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values(
    struct part *restrict p, const struct xpart *restrict xp,
    const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor) {}

/**
 * @brief Extra operations to be done during the drift
 *
 * This predicted the primitive variables a half timestep into the future, but
 * this is better done during the flux calculation (in the gradients predict,
 * TODO).
 *
 * @param p Particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part *p, struct xpart *xp, float dt_drift, float dt_therm,
    float dt_kick_grav, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct pressure_floor_props *pressure_floor) {

  /* skip the drift if we are using Lloyd's algorithm */
  /* TODO */

#ifdef SHADOWSWIFT_EXTRAPOLATE_TIME
  /* Extrapolate primitive quantities to the current time */
  float W[6];
  hydro_part_get_primitive_variables(p, W);
  hydro_gradients_extrapolate_in_time(p, W, dt_therm, p->dW_time);
#endif

  /* Reset the delaunay flags after a particle has been drifted */
  p->geometry.delaunay_flags = 0;
}

/**
 * @brief Set the particle acceleration after the flux loop
 *
 * We use the new conserved variables to calculate the new velocity of the
 * particle, and use that to derive the change of the velocity over the particle
 * time step.
 *
 * If the particle time step is zero, we set the accelerations to zero. This
 * should only happen at the start of the simulation.
 *
 * @param p Particle to act upon.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_force(
    struct part *p, const struct cosmology *cosmo) {

  /* Reset force variables if we are using Lloyd's algorithm. */
  /* TODO */
}

/**
 * @brief Convert conserved variables into primitive variables.
 *
 * NOTE: This function *may violate energy conservation* in the case where the
 * pressure is recovered from the Thermal or Entropy directly instead of from
 * the Total energy
 *
 * @param p The particle to act upon.
 * @param xp The @xpart
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme
 * @param floor_props The properties of the entropy floor.
 * @param Q (inout) The vector of conserved quantities
 * @param W (out) The vector of primitive quantities
 * @param thermal_energy (out) The updated thermal energy of the particle
 */
__attribute__((always_inline)) INLINE static void
hydro_convert_conserved_to_primitive(
    struct part *restrict p, const struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props, float *Q, float *W,
    double *thermal_energy) {

  const float m_inv = (Q[0] != 0.0f) ? 1.0f / Q[0] : 0.0f;
  const float volume_inv = 1.f / p->geometry.volume;

  /* Update density and velocity */
  W[0] = Q[0] * volume_inv;
  hydro_set_velocity_from_momentum(&Q[1], m_inv, W[0], &W[1]);

  /* Calculate the updated internal energy and entropic function A.
   * NOTE: This may violate energy conservation. */
  double u;
#ifdef EOS_ISOTHERMAL_GAS
  u = gas_internal_energy_from_pressure(0.f, 0.f);
#else
  double Ekin = 0.5f * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * m_inv;
  *thermal_energy = Q[4] - Ekin;
#ifdef SHADOWSWIFT_THERMAL_ENERGY_SWITCH
#if SHADOWSWIFT_THERMAL_ENERGY_SWITCH == THERMAL_ENERGY_SWITCH_SPRINGEL
  const float *g = xp->a_grav;
  float Egrav = Q[0] * sqrtf(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) *
                hydro_get_comoving_psize(p);
  if (*thermal_energy > 1e-2 * p->timestepvars.Ekin &&
      *thermal_energy > 1e-2 * Egrav) {
    /* Recover thermal energy and entropy from total energy */
    u = *thermal_energy * m_inv;
  } else {
    /* Keep entropy conserved and recover thermal and total energy. */
    double A = p->conserved.entropy * m_inv;
    u = gas_internal_energy_from_entropy(W[0], A);
  }
#elif SHADOWSWIFT_THERMAL_ENERGY_SWITCH == THERMAL_ENERGY_SWITCH_SPRINGEL_MACH
  if (p->timestepvars.mach_number > 1.1) {
    /* Recover thermal energy and entropy from total energy */
    u = thermal_energy * m_inv;
  } else {
    /* Keep entropy conserved and recover thermal and total energy. */
    double A = p->conserved.entropy * m_inv;
    u = gas_internal_energy_from_entropy(W[0], A);
  }
#elif SHADOWSWIFT_THERMAL_ENERGY_SWITCH == THERMAL_ENERGY_SWITCH_ASENSIO
  const float *g = xp->a_grav;
  float Egrav = Q[0] * sqrtf(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) *
                hydro_get_comoving_psize(p);
  if (thermal_energy > 1e-2 * (Ekin + Egrav)) {
    /* Recover thermal energy and entropy from total energy */
    u = thermal_energy * m_inv;
  } else if (thermal_energy < 1e-3 * (p->timestepvars.Ekin + thermal_energy) ||
             thermal_energy < 1e-3 * Egrav) {
    /* Keep entropy conserved and recover thermal and total energy. */
    double A = p->conserved.entropy * m_inv;
    u = gas_internal_energy_from_entropy(W[0], A);
  } else {
    /* Use evolved thermal energy to set total energy and entropy */
    u = p->thermal_energy * m_inv;
  }
#else
  error("Unknown thermal energy switch!");
#endif
#else
  u = thermal_energy * m_inv;
#endif

  /* Apply entropy floor.
   * Note that this may also violate energy conservation! */
  /* Explicitly update density here already, since the entropy floor depends on
   * it... */
  p->rho = W[0];
  /* Check against entropy floor */
  const float floor_A = entropy_floor(p, cosmo, floor_props);
  u = fmax(u, gas_internal_energy_from_entropy(W[0], floor_A));
  /* Check against absolute minimum */
  u = fmax(u, hydro_props->minimal_internal_energy /
                  cosmo->a_factor_internal_energy);
#endif

  /* Finally update the pressure and update conserved quantities that are
   * affected by our choice of thermal energy. */
  W[4] = gas_pressure_from_internal_energy(W[0], u);
  W[5] = gas_entropy_from_internal_energy(W[0], u);
  *thermal_energy = Q[0] * u;
  Q[4] = Ekin + *thermal_energy;
  Q[5] = Q[0] * W[5];

  /* reset the primitive variables if we are using Lloyd's algorithm */
  /* TODO */
}

/**
 * @brief Carries out the remaining conversion of the hydrodynamic variables
 * from the initial condition file, now that the volume is known and we have
 * access to the cosmology struct.
 *
 * Mass, velocity, momentum and thermal energy (not comoving) have already been
 * set in `hydro_first_init_part(...)`. We still need to:
 *
 * - Convert thermal energy to comoving thermal energy.
 * - Compute total energy from kinetic energy and thermal energy
 * - Convert mass to density
 * - Compute pressure from density and internal energy
 * - Compute entropy from density and interal energy
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Get some quantities that we'll be using */
  float Q[6], W[6];
  hydro_part_get_conserved_variables(p, Q);
  hydro_part_get_primitive_variables(p, W);
  const float m_inv = Q[0] > 0.f ? 1.f / Q[0] : 0.f;
  const float volume_inv = 1.f / p->geometry.volume;

  /* First of all, convert thermal energy to comoving thermal energy. */
  const float thermal_energy =
      p->thermal_energy / cosmo->a_factor_internal_energy;

  /* Total energy */
  Q[4] =
      thermal_energy + 0.5f * m_inv * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]);

  /* Density */
  W[0] = Q[0] * volume_inv;

  /* Pressure and entropy */
#ifdef EOS_ISOTHERMAL_GAS
  /* although the pressure is not formally used anywhere if an isothermal eos
     has been selected, we still make sure it is set to the correct value */
  W[4] = gas_pressure_from_internal_energy(W[0], 0.0f);
  W[5] = gas_entropy_from_pressure(W[0], W[4]);
  Q[5] = Q[0] * W[5];
#else
  const float u = thermal_energy * m_inv;
  W[4] = gas_pressure_from_internal_energy(W[0], u);
  Q[5] = Q[0] * gas_entropy_from_internal_energy(W[0], u);
  W[5] = Q[5] * m_inv;
#endif

  /* Update quantities */
  hydro_part_set_conserved_variables(p, Q);
  hydro_part_set_primitive_variables(p, W);
  p->thermal_energy = thermal_energy;
}

/**
 * @brief Extra operations done during the kick.
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 * @param dt_grav The time-step for this kick (for gravity quantities).
 * @param dt_hydro The time-step for this kick (for hydro quantities).
 * @param dt_kick_corr The time-step for this kick (for gravity corrections).
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part *restrict p, struct xpart *restrict xp, float dt_therm,
    float dt_grav, float dt_grav_mesh, float dt_hydro, float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  if (dt_therm < 0.0f) {
    /* We are reversing a kick1 due to the timestep limiter */
    /* Note on the fluxes: Since a particle can only receive time integrated
     * fluxes over time steps smaller than or equal to its own time step, we do
     * not need any rescaling or other special care. */

#ifdef SWIFT_DEBUG_CHECKS
    assert(p->timestepvars.last_kick == KICK1);
#endif
    /* Kick generator (undo the last half kick) */
    hydro_generator_velocity_half_kick(p, xp, dt_therm);

    /* Signal that we just did a rollback */
    p->timestepvars.last_kick = ROLLBACK;

    /* Reset the flux.dt */
    p->flux.dt = -1.0f;

    /* Nothing else to do here. */
    return;
  }

  if (p->timestepvars.last_kick == KICK1) {
    /* I.e. we are in kick2 (end of timestep), since the dt_therm > 0. */

    /* Add gravity. We only do this if we have gravity activated. */
    if (p->gpart) {
      /* Retrieve the current value of the gravitational acceleration from the
         gpart. We are only allowed to do this because this is the kick. We
         still need to check whether gpart exists though.*/
      float a_grav[3], grav_kick[3];

      a_grav[0] = p->gpart->a_grav[0] + p->gpart->a_grav_mesh[0];
      a_grav[1] = p->gpart->a_grav[1] + p->gpart->a_grav_mesh[1];
      a_grav[2] = p->gpart->a_grav[2] + p->gpart->a_grav_mesh[2];

      float mdt1 = p->gravity.dt * p->conserved.mass;
      float mdt2 = dt_grav * (p->conserved.mass + p->flux.mass);
      grav_kick[0] = mdt2 * a_grav[0] + mdt1 * xp->a_grav[0];
      grav_kick[1] = mdt2 * a_grav[1] + mdt1 * xp->a_grav[1];
      grav_kick[2] = mdt2 * a_grav[2] + mdt1 * xp->a_grav[2];

      /* apply both half kicks to the momentum */
      /* Note that this also affects the particle movement, as the velocity for
         the particles is set after this. */
      p->conserved.momentum[0] += grav_kick[0];
      p->conserved.momentum[1] += grav_kick[1];
      p->conserved.momentum[2] += grav_kick[2];

      /* Extra hydrodynamic energy due to gravity kick, see eq. 94 in Springel
       * (2010) or eq. 62 in theory/Cosmology/cosmology.pdf. */
      /* Divide total integrated mass flux by the timestep for hydrodynamical
       * quantities. We will later multiply with the correct timestep (these
       * differ for cosmological simulations). */
      float dt_grav_corr1 = 0.f;
      float dt_grav_corr2 = 0.f;
      if (p->flux.dt > 0.) {
        dt_grav_corr1 = p->gravity.dt;
        dt_grav_corr2 = dt_grav;
      }
      float dE_springel = hydro_gravity_energy_update_term(
          dt_grav_corr1, dt_grav_corr2, xp->a_grav, a_grav, p->gravity.mflux,
          p->v_full, grav_kick);

      p->conserved.energy += dE_springel;
    }

#ifdef SWIFT_DEBUG_CHECKS
    assert(p->flux.dt >= 0.0f);
#endif

    float Q[6];
    hydro_part_get_conserved_variables(p, Q);
    if (p->flux.dt > 0.0f) {
      /* Apply the fluxes */
      /* We are in kick2 of a normal timestep (not the very beginning of the
       * simulation) */
      float flux[6];
      hydro_part_get_fluxes(p, flux);

      /* Update conserved variables. */
      Q[0] += flux[0];
      Q[1] += flux[1];
      Q[2] += flux[2];
      Q[3] += flux[3];
#if defined(EOS_ISOTHERMAL_GAS)
      /* We use the EoS equation in a sneaky way here just to get the constant
       * internal energy */
      float u = gas_internal_energy_from_entropy(0.0f, 0.0f);
      float rho = Q[0] / p->geometry.volume;
      p->thermal_energy = p->conserved.mass * u;
      float m_inv = Q[0] > 0.f ? 1.f / Q[0] : 0.f;
      Q[4] = p->thermal_energy +
             0.5f * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * m_inv;
      Q[5] = Q[0] * gas_entropy_from_internal_energy(rho, u);
#else
      Q[4] += flux[4];
      Q[5] += flux[5];
      // See eq. 24 in Alonso Asensio et al. (preprint 2023)
      p->thermal_energy +=
          flux[4] -
          (p->v[0] * flux[1] + p->v[1] * flux[2] + p->v[2] * flux[3]) +
          0.5f * (p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]) *
              flux[0];
#endif

#ifndef HYDRO_GAMMA_5_3
      const float Pcorr = (dt_hydro - dt_therm) * p->geometry.volume;
      Q[1] -= Pcorr * p->gradients.P[0];
      Q[2] -= Pcorr * p->gradients.P[1];
      Q[3] -= Pcorr * p->gradients.P[2];
      Q[4] -=
          Pcorr * (p->v[0] * p->gradients.P[0] + p->v[1] * p->gradients.P[1] +
                   p->v[2] * p->gradients.P[2]);
#endif
    }

    /* Compute primitive quantities. Note that this may also modify the vector
     * of conserved quantities (e.g. for cold flows).
     * This also applies entropy floor/minimal internal energy limit. */
    float W[6];
    double thermal_energy;
    hydro_convert_conserved_to_primitive(p, xp, cosmo, hydro_props, floor_props,
                                         Q, W, &thermal_energy);

    /* Update the particle */
    hydro_part_set_conserved_variables(p, Q);
    hydro_part_set_primitive_variables(p, W);
    p->thermal_energy = thermal_energy;
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

#ifndef SHADOWSWIFT_FIX_PARTICLES
    /* Set the generator velocity to the fluid velocity. Steering will be
     * applied after the first half kick */
    xp->v_full[0] = p->v[0];
    xp->v_full[1] = p->v[1];
    xp->v_full[2] = p->v[2];
#endif

    /* Reset the fluxes so that they do not get used again in the kick1. */
    hydro_part_reset_fluxes(p);

#ifdef SHADOWSWIFT_EXTRAPOLATE_TIME
    /* We are at the end of the particles timestep: reset time extrapolations of
     * primitive quantities */
    p->dW_time[0] = 0.0f;
    p->dW_time[1] = 0.0f;
    p->dW_time[2] = 0.0f;
    p->dW_time[3] = 0.0f;
    p->dW_time[4] = 0.0f;
    p->dW_time[5] = 0.0f;
#endif

    /* Signal we just did kick2 */
    p->timestepvars.last_kick = KICK2;

  } else if (p->timestepvars.last_kick == ROLLBACK) {
#ifdef SWIFT_DEBUG_CHECKS
    assert(p->flux.dt == -1.0f);
#endif
    /* Kick generator */
    hydro_generator_velocity_half_kick(p, xp, dt_therm);

    /* Update the flux.dt */
    p->flux.dt = dt_therm;

    /* Update the gravitational dt */
    p->gravity.dt = dt_grav;
    p->gravity.dt_corr = dt_kick_corr;

    /* Signal we just did a restore */
    p->timestepvars.last_kick = RESTORE_AFTER_ROLLBACK;

  } else if (p->timestepvars.last_kick == RESTORE_AFTER_ROLLBACK) {
    /* We are in kick1 after a rollback. */
#ifdef SWIFT_DEBUG_CHECKS
    assert(p->flux.dt >= 0.0f);
#endif
    /* Kick generator */
    hydro_generator_velocity_half_kick(p, xp, dt_therm);

    /* Add the remainder of this particle's timestep to flux.dt */
    p->flux.dt += 2.f * dt_therm;

    /* Add the remainder of the first half kick to gravity timesteps */
    p->gravity.dt += dt_grav;
    p->gravity.dt_corr += dt_kick_corr;

    /* Now that we have received both half kicks, we can set the actual
     * velocity of the ShadowSWIFT particle (!= fluid velocity) */
    hydro_velocities_set(p, xp, p->flux.dt);

    /* Signal we just did a kick1 */
    p->timestepvars.last_kick = KICK1;

  } else if (p->timestepvars.last_kick == KICK2) {
    /* We are in kick1 after a kick2 (normal scenario). */
#ifdef SWIFT_DEBUG_CHECKS
    assert(p->flux.dt == -1.0f);
#endif
    /* Kick generator */
    hydro_generator_velocity_half_kick(p, xp, dt_therm);

    /* Update the time step used in the flux calculation */
    p->flux.dt = 2.f * dt_therm;

    /* Set the timestep for the first half kick (gravity) */
    p->gravity.dt = dt_grav;
    p->gravity.dt_corr = dt_kick_corr;

    /* Reset the timestep_vars for the next timestep */
    hydro_reset_timestep_vars(p);

    /* Now that we have received both half kicks, we can set the actual
     * velocity of the ShadowSWIFT particle (!= fluid velocity) */
    hydro_velocities_set(p, xp, p->flux.dt);

    /* Signal we just did a kick1 */
    p->timestepvars.last_kick = KICK1;

  } else {
    error("Impossible scenario!");
  }

  /* undo the flux exchange and kick the particles towards their centroid */
  /* TODO Lloyd */
}

/**
 * @brief Split additional hydrodynamic quantities when splitting particles.
 *
 * Note that the mass is already handled in
 * @engine_split_gas_particle_split_mapper, so this will be a no-op for most
 * hydro schemes.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @param particle_split_factor Over how many particles are we splitting?
 */
__attribute__((always_inline)) INLINE static void hydro_split_part(
    struct part *p, struct xpart *xp, const int particle_split_factor) {
  const float fraction = 1.f / (float)particle_split_factor;
  /* Conserved quantities (without mass) */
  p->conserved.momentum[0] *= fraction;
  p->conserved.momentum[2] *= fraction;
  p->conserved.momentum[1] *= fraction;
  p->conserved.energy *= fraction;
  p->conserved.entropy *= fraction;
  /* Any not applied fluxes */
  p->flux.mass *= fraction;
  p->flux.momentum[0] *= fraction;
  p->flux.momentum[1] *= fraction;
  p->flux.momentum[2] *= fraction;
  p->flux.energy *= fraction;
  p->flux.entropy *= fraction;
  /* Volume */
  p->geometry.volume *= fraction;
  /* Temporarily set centroid to position.
   * NOTE: This is incorrect, but at least ensures that the centroid will be in
   * inside the particles new Voronoi cell.
   * NOTE: In shadowswift, the centroid field stores an offset of the particle's
   * position */
  p->geometry.centroid[0] = 0.f;
  p->geometry.centroid[1] = 0.f;
  p->geometry.centroid[2] = 0.f;
}

/**
 * @brief Update given random displacement vector if needed.
 *
 * Might be necessary for schemes with asymmetric cells.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @param displacement (in-out) initial random displacement vector.
 */
__attribute__((always_inline)) INLINE static void hydro_split_part_displacement(
    struct part *p, struct xpart *xp, double *displacement) {
  /* Set some minimal distance */
  double displacement_nrm = sqrt(displacement[0] * displacement[0] +
                                 displacement[1] * displacement[1] +
                                 displacement[2] * displacement[2]);
#ifdef SWIFT_DEBUG_CHECKS
  if (displacement_nrm == 0) {
    error("Displacement vector cannot be 0!");
  }
#endif
  /* Rescale to make displacement smaller */
  double fac = fmax(1e-4, 1e-8 / displacement_nrm);
  displacement[0] *= fac;
  displacement[1] *= fac;
  displacement[2] *= fac;

  /* make perpendicular to cell axis */
  const float *axis = p->geometry.centroid;
  double axis_nrm2 = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
  if (axis_nrm2 < 1e-16) {
    /* Centroidal cell is symmetric */
    return;
  }
  double delta_dot_axis = displacement[0] * axis[0] +
                          displacement[1] * axis[1] + displacement[2] * axis[2];
  double delta_x = displacement[0] - axis[0] * delta_dot_axis / axis_nrm2;
  double delta_y = displacement[1] - axis[1] * delta_dot_axis / axis_nrm2;
  double delta_z = displacement[2] - axis[2] * delta_dot_axis / axis_nrm2;

  if (delta_x * delta_x + delta_y * delta_y + delta_z * delta_z < 1e-17) {
    /* try again with cyclic permutation */
    delta_dot_axis = displacement[1] * axis[0] + displacement[2] * axis[1] +
                     displacement[0] * axis[2];
    delta_x = displacement[0] - axis[0] * delta_dot_axis / axis_nrm2;
    delta_y = displacement[1] - axis[1] * delta_dot_axis / axis_nrm2;
    delta_z = displacement[2] - axis[2] * delta_dot_axis / axis_nrm2;
  }

  displacement[0] = delta_x;
  displacement[1] = delta_y;
  displacement[2] = delta_z;
}

/**
 * @brief Operations performed when a particle gets removed from the
 * simulation volume.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @param time The simulation time.
 */
__attribute__((always_inline)) INLINE static void hydro_remove_part(
    const struct part *p, const struct xpart *xp, const double time) {}

#endif /* SWIFT_SHADOWSWIFT_HYDRO_H */
