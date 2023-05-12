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
      sqrtf(hydro_gamma * W[4] / W[0]);
  vmax = max(vmax, p->timestepvars.vmax);

  float psize = hydro_get_physical_psize(p, cosmo);
  if (p->geometry.min_face_dist < 0.25 * psize)
    psize = p->geometry.min_face_dist;

  float dt = FLT_MAX;
  if (vmax > 0.0f) {
    dt = cosmo->a * psize / vmax;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (dt == 0.f) error("Part wants dt=0!");
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

  float W[6], Q[6];

  W[0] = 0.0f;
  W[1] = p->v[0];
  W[2] = p->v[1];
  W[3] = p->v[2];
  W[4] = 0.0f;
  W[5] = 0.0f;

#ifdef EOS_ISOTHERMAL_GAS
  p->thermal_energy = Q[0] * gas_internal_energy_from_entropy(0.0f, 0.0f);
#else
  p->thermal_energy = p->thermal_energy * p->conserved.mass;
#endif

  Q[0] = p->conserved.mass;
  Q[1] = Q[0] * W[1];
  Q[2] = Q[0] * W[2];
  Q[3] = Q[0] * W[3];
  Q[4] = 0.0f;  // We need the comoving internal energy to compute the total
                // comoving energy, so we cannot do this here...
  Q[5] = 0.0f;  // We need the volume to be able to compute the entropy...?

  /* overwrite all hydro variables if we are using Lloyd's algorithm */
  /* TODO */

  p->time_bin = 0;

  hydro_part_set_primitive_variables(p, W);
  hydro_part_set_conserved_variables(p, Q);

  /* initialize the particle velocity based on the primitive fluid velocity */
  hydro_velocities_init(p, xp);

  /* ignore accelerations present in the initial condition */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  p->flux_count = 0;
  p->geometry.delaunay_flags = 0;
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

  /* Prepare the slope limiter for this particle */
  hydro_slope_limiter_prepare(p);
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * This function is called in the ghost task to convert some quantities coming
 * from the density loop over neighbours into quantities ready to be used in the
 * force loop over neighbours. Quantities are typically read from the density
 * sub-structure and written to the force sub-structure.
 * Examples of calculations done here include the calculation of viscosity term
 * constants, thermal conduction terms, hydro conversions, etc.
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
  hydro_part_reset_fluxes(p);
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
  /* Extrapolate primitive quantities in time */
  float W[6];
  hydro_part_get_primitive_variables(p, W);
  hydro_gradients_extrapolate_in_time(p, W, 0.5f * dt_therm, p->dW_time);

  // MATTHIEU: Apply the entropy floor here.
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
 * @param volume The volume of the particle's associated voronoi cell
 */
__attribute__((always_inline)) INLINE static void
hydro_convert_conserved_to_primitive(struct part *p, struct xpart *xp,
                                     const struct cosmology *cosmo) {
  float W[6], Q[6];
  hydro_part_get_conserved_variables(p, Q);
  const float m_inv = (Q[0] != 0.0f) ? 1.0f / Q[0] : 0.0f;
  const float volume_inv = 1.f / p->geometry.volume;

  W[0] = Q[0] * volume_inv;
  hydro_set_velocity_from_momentum(&Q[1], m_inv, W[0], &W[1]);

#ifdef EOS_ISOTHERMAL_GAS
  /* although the pressure is not formally used anywhere if an isothermal eos
     has been selected, we still make sure it is set to the correct value */
  W[4] = gas_pressure_from_internal_energy(W[0], 0.0f);
#else

  /* Calculate the pressure from the internal energy, make sure that the entropy
   * and total energy stay consistent with our choice of thermal energy.
   * NOTE: This may violate energy conservation. */
  float Ekin = 0.5f * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * m_inv;
  float *g = xp->a_grav;
  float Egrav = Q[0] * sqrtf(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) *
                hydro_get_comoving_psize(p);
  float thermal_energy = Q[4] - Ekin;
  if (thermal_energy > 1e-2 * (Ekin + Egrav)) {
    /* Recover thermal energy and entropy from total energy */
    p->thermal_energy = thermal_energy;
    W[5] = gas_entropy_from_internal_energy(W[0], thermal_energy * m_inv);
    p->conserved.entropy = Q[0] * W[5];
  } else if (thermal_energy < 1e-3 * p->limiter.Ekin ||
             thermal_energy < 1e-3 * Egrav) {
    /* Keep entropy conserved and recover thermal and total energy. */
    W[5] = p->conserved.entropy * m_inv;
    p->thermal_energy = Q[0] * gas_internal_energy_from_entropy(W[0], W[5]);
    p->conserved.energy = Ekin + p->thermal_energy;
  } else {
    /* Use evolved thermal energy to set entropy and total energy */
    p->conserved.energy = Ekin + p->thermal_energy;
    W[5] = gas_entropy_from_internal_energy(W[0], p->thermal_energy * m_inv);
    p->conserved.entropy = Q[0] * W[5];
  }
  /* Calculate pressure from thermal energy */
  W[4] = gas_pressure_from_internal_energy(W[0], p->thermal_energy * m_inv);
#endif
  /* reset the primitive variables if we are using Lloyd's algorithm */
  /* TODO */

  hydro_part_set_primitive_variables(p, W);

  if (m_inv == 0. && (p->v[0] != 0. || p->v[1] != 0. || p->v[2] != 0.)) {
    error("Nonzero v for particle with zero mass!");
  }

  if (p->rho < 0.) {
    error("Negative density!");
  }

  if (p->P < 0.) {
    error("Negative pressure!");
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (p->rho != p->rho) error("NaN density!");
  if (p->v[0] != p->v[0]) error("NaN vx!");
  if (p->v[1] != p->v[1]) error("NaN vy!");
  if (p->v[2] != p->v[2]) error("NaN vz!");
  if (p->P != p->P) error("NaN pressure!");
  if (p->A != p->A) error("NaN entropic function!");
#endif
}

/**
 * @brief Converts the hydrodynamic variables from the initial condition file to
 * conserved variables that can be used during the integration (flux
 * calculation).
 *
 * Requires the volume to be known.
 *
 * The initial condition file contains a mixture of primitive and conserved
 * variables. Mass is a conserved variable, and we just copy the particle
 * mass into the corresponding conserved quantity. We need the volume to
 * also derive a density, which is then used to convert the internal energy
 * to a pressure. However, we do not actually use these variables anymore.
 * We do need to initialize the linear momentum, based on the mass and the
 * velocity of the particle.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Convert thermal energy to comoving thermal energy, we have to do this here
   * because we do not have access to the cosmology struct in
   * hydro_first_init_part()... */
  p->thermal_energy /= cosmo->a_factor_internal_energy;
  float Q[6];
  hydro_part_get_conserved_variables(p, Q);
  float m_inv = Q[0] > 0.f ? 1.f / Q[0] : 0.f;
  p->conserved.energy =
      p->thermal_energy +
      0.5f * m_inv * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]);

  shadowswift_check_physical_quantities(
      "mass", "energy", p->conserved.mass, p->conserved.momentum[0],
      p->conserved.momentum[1], p->conserved.momentum[2], p->conserved.energy,
      p->conserved.entropy);

  hydro_convert_conserved_to_primitive(p, xp, cosmo);
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

  /* Add gravity. We only do this if we have gravity activated. */
  if (p->gpart) {
    /* Retrieve the current value of the gravitational acceleration from the
       gpart. We are only allowed to do this because this is the kick. We still
       need to check whether gpart exists though.*/
    float a_grav[3], grav_kick_factor[3];

    a_grav[0] = p->gpart->a_grav[0] + p->gpart->a_grav_mesh[0];
    a_grav[1] = p->gpart->a_grav[1] + p->gpart->a_grav_mesh[1];
    a_grav[2] = p->gpart->a_grav[2] + p->gpart->a_grav_mesh[2];

    grav_kick_factor[0] = dt_grav * p->gpart->a_grav[0];
    grav_kick_factor[1] = dt_grav * p->gpart->a_grav[1];
    grav_kick_factor[2] = dt_grav * p->gpart->a_grav[2];
    if (dt_grav_mesh != 0) {
      grav_kick_factor[0] += dt_grav_mesh * p->gpart->a_grav_mesh[0];
      grav_kick_factor[1] += dt_grav_mesh * p->gpart->a_grav_mesh[1];
      grav_kick_factor[2] += dt_grav_mesh * p->gpart->a_grav_mesh[2];
    }

    /* Kick the momentum for half a time step */
    /* Note that this also affects the particle movement, as the velocity for
       the particles is set after this. */
    p->conserved.momentum[0] += p->conserved.mass * grav_kick_factor[0];
    p->conserved.momentum[1] += p->conserved.mass * grav_kick_factor[1];
    p->conserved.momentum[2] += p->conserved.mass * grav_kick_factor[2];

    /* Extra *kinetic* energy due to gravity kick */
    float gravity_work_therm = hydro_gravity_energy_update_term(
        dt_kick_corr, p, p->conserved.momentum, a_grav, grav_kick_factor);
    p->conserved.energy += gravity_work_therm;
  }

  if (dt_therm < 0.0f) {
    /* We are reversing a kick1 due to the timestep limiter */
    /* Note on the fluxes: Since a particle can only receive time integrated
     * fluxes over time steps smaller than or equal to its own time step, we do
     * not need any rescaling or other special care. */

#ifdef SWIFT_DEBUG_CHECKS
    assert(p->timestepvars.last_kick == KICK1);
#endif

    /* Signal that we just did a rollback */
    p->timestepvars.last_kick = ROLLBACK;

    /* Reset the flux.dt */
    p->flux.dt = -1.0f;

    /* Nothing else to do here. */
    return;
  }

  if (p->timestepvars.last_kick == KICK1) {
    /* I.e. we are in kick2 (end of timestep), since the dt_therm > 0. */

#ifdef SWIFT_DEBUG_CHECKS
    assert(p->flux.dt >= 0.0f);
#endif

    if (p->flux.dt > 0.0f) {
      /* We are in kick2 of a normal timestep (not the very beginning of the
       * simulation) */
      float flux[6];
      hydro_part_get_fluxes(p, flux);

      /* Update conserved variables. */
      p->conserved.mass += flux[0];
      p->conserved.momentum[0] += flux[1];
      p->conserved.momentum[1] += flux[2];
      p->conserved.momentum[2] += flux[3];
#if defined(EOS_ISOTHERMAL_GAS)
      /* We use the EoS equation in a sneaky way here just to get the constant u
       */
      float u = gas_internal_energy_from_entropy(0.0f, 0.0f);
      p->thermal_energy = p->conserved.mass * u;
      p->conserved.energy =
          p->thermal_energy +
          0.5f *
              (p->conserved.momentum[0] * p->conserved.momentum[0] +
               p->conserved.momentum[1] * p->conserved.momentum[1] +
               p->conserved.momentum[2] * p->conserved.momentum[2]) /
              p->conserved.mass;
      p->conserved.entropy =
          p.conserved.mass * gas_entropy_from_internal_energy(
                                 p->conserved.mass / p->geometry.volume, u);
#else
      p->conserved.energy += flux[4];
      p->conserved.entropy += flux[5];
      // See eq. 24 in Alonso Asensio et al. (preprint 2023)
      p->thermal_energy +=
          flux[4] -
          (p->v[0] * flux[1] + p->v[1] * flux[2] + p->v[2] * flux[3]) +
          0.5f * (p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]) *
              flux[0];
#endif

#ifndef HYDRO_GAMMA_5_3

      const float Pcorr = (dt_hydro - dt_therm) * p->geometry.volume;
      p->conserved.momentum[0] -= Pcorr * p->gradients.P[0];
      p->conserved.momentum[1] -= Pcorr * p->gradients.P[1];
      p->conserved.momentum[2] -= Pcorr * p->gradients.P[2];
      p->conserved.energy -=
          Pcorr * (p->v[0] * p->gradients.P[0] + p->v[1] * p->gradients.P[1] +
                   p->v[2] * p->gradients.P[2]);
#endif

      /* Check conserved quantities */
      shadowswift_check_physical_quantities(
          "mass", "energy", p->conserved.mass, p->conserved.momentum[0],
          p->conserved.momentum[1], p->conserved.momentum[2],
          p->conserved.energy, p->conserved.entropy);

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

      /* Now that the mass is updated, update gpart mass accordingly */
      hydro_gravity_update_gpart_mass(p);
    }

    /* Update primitive quantities. Note that this also updates the fluid
     * velocity p->v. */
    hydro_convert_conserved_to_primitive(p, xp, cosmo);

    /* Reset the fluxes so that they do not get used again in the kick1. */
    hydro_part_reset_fluxes(p);

#ifdef SHADOWSWIFT_EXTRAPOLATE_TIME
    /* Reset time extrapolations of primitive quantities */
    p->dW_time[0] = 0.0f;
    p->dW_time[1] = 0.0f;
    p->dW_time[2] = 0.0f;
    p->dW_time[3] = 0.0f;
    p->dW_time[4] = 0.0f;
    p->dW_time[5] = 0.0f;
#endif

    /* Apply the minimal energy limit */
    const float min_energy =
        hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;
    if (p->thermal_energy < min_energy * p->conserved.mass) {
      hydro_set_comoving_internal_energy(p, min_energy);
    }

    // MATTHIEU: Apply the entropy floor here.

    /* Signal we just did kick2 */
    p->timestepvars.last_kick = KICK2;

  } else if (p->timestepvars.last_kick == ROLLBACK) {
#ifdef SWIFT_DEBUG_CHECKS
    assert(p->flux.dt == -1.0f);
#endif
    /* Update the flux.dt */
    p->flux.dt = dt_therm;

    /* Signal we just did a restore */
    p->timestepvars.last_kick = RESTORE_AFTER_ROLLBACK;

  } else if (p->timestepvars.last_kick == RESTORE_AFTER_ROLLBACK) {
    /* We are in kick1 after a rollback. */
#ifdef SWIFT_DEBUG_CHECKS
    assert(p->flux.dt >= 0.0f);
#endif

    /* Add the remainder of this particle's timestep to flux.dt */
    p->flux.dt += 2.f * dt_therm;

    /* Now that we have received both half kicks, we can set the actual
     * velocity of the ShadowSWIFT particle (!= fluid velocity) */
    hydro_velocities_set(p, xp);

    /* Signal we just did a kick1 */
    p->timestepvars.last_kick = KICK1;

  } else if (p->timestepvars.last_kick == KICK2) {
    /* We are in kick1 after a kick2 (normal scenario). */
#ifdef SWIFT_DEBUG_CHECKS
    assert(p->flux.dt == -1.0f);
#endif

    /* Update the time step used in the flux calculation */
    p->flux.dt = 2.f * dt_therm;

    /* Reset v_max */
    p->timestepvars.vmax = 0.f;

    /* Reset Ekin */
    p->limiter.Ekin = -INFINITY;

    /* Now that we have received both half kicks, we can set the actual
     * velocity of the ShadowSWIFT particle (!= fluid velocity) */
    hydro_velocities_set(p, xp);

    /* Signal we just did a kick1 */
    p->timestepvars.last_kick = KICK1;

  } else {
    error("Impossible scenario!");
  }

  /* undo the flux exchange and kick the particles towards their centroid */
  /* TODO Lloyd */
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
