/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/**
 * @file Minimal/hydro.h
 * @brief Minimal conservative implementation of SPH (Non-neighbour loop
 * equations)
 *
 * The thermal variable is the internal energy (u). Simple constant
 * viscosity term with the Balsara (1995) switch. No thermal conduction
 * term is implemented.
 *
 * This corresponds to equations (43), (44), (45), (101), (103)  and (104) with
 * \f$\beta=3\f$ and \f$\alpha_u=0\f$ of Price, D., Journal of Computational
 * Physics, 2012, Volume 231, Issue 3, pp. 759-794.
 */

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
#include "hydro_velocities.h"
#include "kernel_hydro.h"
#include "minmax.h"

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

  float W[5];
  hydro_part_get_primitive_variables(p, W);

  /* v_full is the actual velocity of the particle, v is its
     hydrodynamical velocity. The time step depends on the relative difference
     of the two. */
  float vrel[3];
  vrel[0] = W[1] - xp->v_full[0];
  vrel[1] = W[2] - xp->v_full[1];
  vrel[2] = W[3] - xp->v_full[2];
  float vmax =
      sqrtf(vrel[0] * vrel[0] + vrel[1] * vrel[1] + vrel[2] * vrel[2]) +
      sqrtf(hydro_gamma * W[4] / W[0]);
  vmax = max(vmax, p->timestepvars.vmax);

  const float psize = cosmo->a * cosmo->a *
                      powf(p->geometry.volume / hydro_dimension_unit_sphere,
                           hydro_dimension_inv);
  float dt = FLT_MAX;
  if (vmax > 0.0f) {
    dt = psize / vmax;
  }
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

  float W[5], Q[5];

  W[0] = 0.0f;
  W[1] = p->v[0];
  W[2] = p->v[1];
  W[3] = p->v[2];
  W[4] = 0.0f;

  Q[0] = p->conserved.mass;
  Q[1] = Q[0] * W[1];
  Q[2] = Q[0] * W[2];
  Q[3] = Q[0] * W[3];
#if defined(EOS_ISOTHERMAL_GAS)
  Q[4] = Q[0] * gas_internal_energy_from_entropy(0.0f, 0.0f);
#else
  Q[4] = p->conserved.energy * Q[0];
#endif

#ifdef SHADOWSWIFT_TOTAL_ENERGY
  Q[4] += 0.5f * (Q[1] * W[1] + Q[2] * W[2] + Q[3] * W[3]);
#endif

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
 * We use this to set the timestep used in the flux calculation.
 *
 * @param p The particle to act upon.
 * @param dt Physical time step of the particle during the next step.
 */
__attribute__((always_inline)) INLINE static void hydro_timestep_extra(
    struct part *p, float dt) {}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Not used in the ShadowSWIFT scheme.
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
    const struct cosmology *cosmo, const struct hydro_props *hydro_props) {

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
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const float dt_alpha) {
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
    const struct cosmology *cosmo) {
  // MATTHIEU: Apply the entropy floor here.
}

/**
 * @brief Converts the hydrodynamic variables from the initial condition file to
 * conserved variables that can be used during the integration
 *
 * We no longer do this, as the mass needs to be provided in the initial
 * condition file, and the mass alone is enough to initialize all conserved
 * variables. This is now done in hydro_first_init_part.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props) {

  p->conserved.energy /= cosmo->a_factor_internal_energy;
}

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
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* skip the drift if we are using Lloyd's algorithm */
  /* TODO */

  float W[5];
  hydro_part_get_primitive_variables(p, W);

  // MATTHIEU: Apply the entropy floor here.

  /* add the gravitational contribution to the fluid velocity drift */
  /* The fluid velocity is now drifted directly by drift_part() */
  //  hydro_gravity_extra_velocity_drift(&W[1], p->v, xp->v_full);

  /* TODO Predict primitive variables forward? Would be with outdated gradients,
   * so maybe better in flux? */

  hydro_part_set_primitive_variables(p, W);

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
 * @param p The particle to act upon.
 * @param volume The volume of the particle's associated voronoi cell
 */
__attribute__((always_inline)) INLINE static void
hydro_convert_conserved_to_primitive(struct part *restrict p,
                                     struct xpart *restrict xp) {

  float W[5], Q[5];
  hydro_part_get_conserved_variables(p, Q);
  const float m_inv = (Q[0] != 0.0f) ? 1.0f / Q[0] : 0.0f;
  const float volume_inv = 1.f / p->geometry.volume;

  W[0] = Q[0] * volume_inv;
  hydro_velocities_from_momentum(&Q[1], m_inv, W[0], &W[1]);

#ifdef EOS_ISOTHERMAL_GAS
  /* although the pressure is not formally used anywhere if an isothermal eos
     has been selected, we still make sure it is set to the correct value */
  W[4] = gas_pressure_from_internal_energy(W[0], 0.0f);
#else

#ifdef SHADOWSWIFT_TOTAL_ENERGY
  /* subtract the kinetic energy; we want the thermal energy */
  Q[4] -= 0.5f * (Q[1] * W[1] + Q[2] * W[2] + Q[3] * W[3]);
#endif

  W[4] = gas_pressure_from_internal_energy(W[0], Q[4] * m_inv);
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
    float dt_grav, float dt_hydro, float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Add gravity. We only do this if we have gravity activated. */
  if (p->gpart) {
    /* Retrieve the current value of the gravitational acceleration from the
       gpart. We are only allowed to do this because this is the kick. We still
       need to check whether gpart exists though.*/
    float a_grav[3];

    /* TODO also add mesh acceleration? */
    a_grav[0] = p->gpart->a_grav[0];
    a_grav[1] = p->gpart->a_grav[1];
    a_grav[2] = p->gpart->a_grav[2];

    p->conserved.energy += hydro_gravity_energy_update_term(
        dt_kick_corr, dt_grav, p, p->conserved.momentum, a_grav);

    /* Kick the momentum for half a time step */
    /* Note that this also affects the particle movement, as the velocity for
       the particles is set after this. */
    p->conserved.momentum[0] += p->conserved.mass * a_grav[0] * dt_grav;
    p->conserved.momentum[1] += p->conserved.mass * a_grav[1] * dt_grav;
    p->conserved.momentum[2] += p->conserved.mass * a_grav[2] * dt_grav;
  }

  if (p->flux.dt > 0.0f) {
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

    /* Update primitive quantities. Note that this also updates the fluid
     * velocity p->v. */
    hydro_convert_conserved_to_primitive(p, xp);

    /* Update gpart mass */
    hydro_gravity_update_gpart_mass(p);

  } else if (p->flux.dt == 0.0f) {
    /* This can only occur in kick2 the beginning of the simulation. We simply
     * need to reset the flux.dt to -1.0f and calculate the primitive
     * quantities. */

    /* Reset the fluxes so that they do not get used again in the kick1. */
    hydro_part_reset_fluxes(p);

    /* Update primitive quantities */
    hydro_convert_conserved_to_primitive(p, xp);

  } else {
    /* flux.dt < 0 implies that we are in kick1. */

    /* Update the timestep used in the flux calculation */
    p->flux.dt = 2.f * dt_therm;

    /* Reset v_max */
    p->timestepvars.vmax = 0.f;

    /* Now that we have received both half kicks, we can set the actual
     * velocity of the ShadowSWIFT particle (!= fluid velocity) */
    hydro_velocities_set(p, xp);
  }

  /* TODO: check for negative dt_therm (implying that we are rolling back a kick
   * due to the timestep limiter or the timestep sync). This means that we must
   * rescale the time integrated fluxes this particle has received. */

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
