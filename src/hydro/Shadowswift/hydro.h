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

#include <float.h>
#include "adiabatic_index.h"
#include "approx_math.h"
#include "cosmology.h"
#include "equation_of_state.h"
#include "hydro_gradients.h"
#include "hydro_properties.h"
#include "hydro_space.h"
#include "voronoi_algorithm.h"

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param hydro_properties Pointer to the hydro parameters.
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    const struct part* restrict p, const struct xpart* restrict xp,
    const struct hydro_props* restrict hydro_properties,
    const struct cosmology* restrict cosmo) {

  const float CFL_condition = hydro_properties->CFL_condition;

  float vrel[3];
  vrel[0] = p->primitives.v[0] - xp->v_full[0];
  vrel[1] = p->primitives.v[1] - xp->v_full[1];
  vrel[2] = p->primitives.v[2] - xp->v_full[2];
  float vmax =
      sqrtf(vrel[0] * vrel[0] + vrel[1] * vrel[1] + vrel[2] * vrel[2]) +
      sqrtf(hydro_gamma * p->primitives.P / p->primitives.rho);
  vmax = max(vmax, p->timestepvars.vmax);

  const float psize =
      cosmo->a *
      powf(p->cell.volume / hydro_dimension_unit_sphere, hydro_dimension_inv);
  float dt = FLT_MAX;
  if (vmax > 0.) {
    dt = psize / vmax;
  }
  return CFL_condition * dt;
}

/**
 * @brief Does some extra hydro operations once the actual physical time step
 * for the particle is known.
 *
 * We use this to store the physical time step, since it is used for the flux
 * exchange during the force loop.
 *
 * We also set the active flag of the particle to inactive. It will be set to
 * active in hydro_init_part, which is called the next time the particle becomes
 * active.
 *
 * @param p The particle to act upon.
 * @param dt Physical time step of the particle during the next step.
 */
__attribute__((always_inline)) INLINE static void hydro_timestep_extra(
    struct part* p, float dt) {

  p->force.dt = dt;
  p->force.active = 0;
}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * In this case, we copy the particle velocities into the corresponding
 * primitive variable field. We do this because the particle velocities in GIZMO
 * can be independent of the actual fluid velocity. The latter is stored as a
 * primitive variable and integrated using the linear momentum, a conserved
 * variable.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_first_init_part(
    struct part* p, struct xpart* xp) {

  const float mass = p->conserved.mass;

  p->time_bin = 0;
  p->wakeup = time_bin_not_awake;

  p->primitives.v[0] = p->v[0];
  p->primitives.v[1] = p->v[1];
  p->primitives.v[2] = p->v[2];

  p->conserved.momentum[0] = mass * p->primitives.v[0];
  p->conserved.momentum[1] = mass * p->primitives.v[1];
  p->conserved.momentum[2] = mass * p->primitives.v[2];

#ifdef EOS_ISOTHERMAL_GAS
  p->conserved.energy = mass * gas_internal_energy_from_entropy(0.f, 0.f);
#else
  p->conserved.energy *= mass;
#endif

#ifdef SHADOWFAX_TOTAL_ENERGY
  p->conserved.energy += 0.5f * (p->conserved.momentum[0] * p->primitives.v[0] +
                                 p->conserved.momentum[1] * p->primitives.v[1] +
                                 p->conserved.momentum[2] * p->primitives.v[2]);
#endif

#if defined(SHADOWFAX_FIX_CELLS)
  p->v[0] = 0.;
  p->v[1] = 0.;
  p->v[2] = 0.;
#else
  p->v[0] = p->primitives.v[0];
  p->v[1] = p->primitives.v[1];
  p->v[2] = p->primitives.v[2];
#endif

  /* set the initial velocity of the cells */
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];

  /* ignore accelerations present in the initial condition */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;
}

/**
 * @brief Prepares a particle for the volume calculation.
 *
 * Simply makes sure all necessary variables are initialized to zero.
 * Initializes the Voronoi cell.
 *
 * @param p The particle to act upon
 * @param hs #hydro_space containing extra information about the space.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part* p, const struct hydro_space* hs) {

  /* make sure we don't enter the no neighbour case in runner.c */
  p->density.wcount = 1.0f;
  p->density.wcount_dh = 0.0f;

  voronoi_cell_init(&p->cell, p->x, hs->anchor, hs->side);

  /* Set the active flag to active. */
  p->force.active = 1;
}

/**
 * @brief Finishes the volume calculation.
 *
 * Calls the finalize method on the Voronoi cell, which calculates the volume
 * and centroid of the cell. We use the return value of this function to set
 * a new value for the smoothing length and possibly force another iteration
 * of the volume calculation for this particle. We then use the volume to
 * convert conserved variables into primitive variables.
 *
 * This method also initializes the gradient variables (if gradients are used).
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part* restrict p, const struct cosmology* cosmo) {

  float volume;
  float m, momentum[3], energy;

  hydro_gradients_init(p);

  float hnew = voronoi_cell_finalize(&p->cell);
  /* Enforce hnew as new smoothing length in the iteration
     This is annoyingly difficult, as we do not have access to the variables
     that govern the loop...
     So here's an idea: let's force in some method somewhere that makes sure
     r->e->hydro_properties->target_neighbours is 1, and
     r->e->hydro_properties->delta_neighbours is 0.
     This way, we can accept an iteration by setting p->density.wcount to 1.
     To get the right correction for h, we set wcount to something else
     (say 0), and then set p->density.wcount_dh to 1/(hnew-h). */
  if (hnew < p->h) {
    /* Iteration succesful: we accept, but manually set h to a smaller value
       for the next time step */
    const float hinvdim = pow_dimension(1.0f / p->h);
    p->density.wcount = hinvdim;
    p->h = 1.1f * hnew;
  } else {
    /* Iteration not succesful: we force h to become 1.1*hnew */
    p->density.wcount = 0.0f;
    p->density.wcount_dh = p->h / (1.1f * hnew - p->h);
    return;
  }
  volume = p->cell.volume;

#ifdef SWIFT_DEBUG_CHECKS
  /* the last condition checks for NaN: a NaN value always evaluates to false,
     even when checked against itself */
  if (volume == 0. || volume == INFINITY || volume != volume) {
    error("Invalid value for cell volume (%g)!", volume);
  }
#endif

  /* compute primitive variables */
  /* eqns (3)-(5) */
  m = p->conserved.mass;
  if (m > 0.) {
    momentum[0] = p->conserved.momentum[0];
    momentum[1] = p->conserved.momentum[1];
    momentum[2] = p->conserved.momentum[2];
    p->primitives.rho = m / volume;
    p->primitives.v[0] = momentum[0] / m;
    p->primitives.v[1] = momentum[1] / m;
    p->primitives.v[2] = momentum[2] / m;

    energy = p->conserved.energy;

#ifdef SHADOWFAX_TOTAL_ENERGY
    energy -= 0.5f * (momentum[0] * p->primitives.v[0] +
                      momentum[1] * p->primitives.v[1] +
                      momentum[2] * p->primitives.v[2]);
#endif

    energy /= m;

    p->primitives.P =
        gas_pressure_from_internal_energy(p->primitives.rho, energy);
  } else {
    p->primitives.rho = 0.;
    p->primitives.v[0] = 0.;
    p->primitives.v[1] = 0.;
    p->primitives.v[2] = 0.;
    p->primitives.P = 0.;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (p->primitives.rho < 0.) {
    error("Negative density!");
  }

  if (p->primitives.P < 0.) {
    error("Negative pressure!");
  }
#endif
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_part_has_no_neighbours(
    struct part* restrict p, struct xpart* restrict xp,
    const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  p->density.wcount = kernel_root * h_inv_dim;
  p->density.wcount_dh = 0.f;
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
    struct part* restrict p, struct xpart* restrict xp,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const float dt_alpha) {

  /* Initialize time step criterion variables */
  p->timestepvars.vmax = 0.0f;

  /* Set the actual velocity of the particle */
  p->force.v_full[0] = xp->v_full[0];
  p->force.v_full[1] = xp->v_full[1];
  p->force.v_full[2] = xp->v_full[2];

  p->conserved.flux.mass = 0.0f;
  p->conserved.flux.momentum[0] = 0.0f;
  p->conserved.flux.momentum[1] = 0.0f;
  p->conserved.flux.momentum[2] = 0.0f;
  p->conserved.flux.energy = 0.0f;
}

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
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_gradient(
    struct part* restrict p, struct xpart* restrict xp,
    const struct cosmology* cosmo) {

  /* Initialize time step criterion variables */
  p->timestepvars.vmax = 0.;
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
    struct part* restrict p) {}

/**
 * @brief Finishes the gradient calculation.
 *
 * Just a wrapper around hydro_gradients_finalize, which can be an empty method,
 * in which case no gradients are used.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_end_gradient(
    struct part* p) {

  hydro_gradients_finalize(p);
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is actually not necessary for Shadowswift, since we just set the
 * accelerations after the flux calculation.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_acceleration(
    struct part* p) {

  /* Reset the acceleration. */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* Reset the time derivatives. */
  p->force.h_dt = 0.0f;
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values(
    struct part* restrict p, const struct xpart* restrict xp) {}

/**
 * @brief Converts the hydrodynamic variables from the initial condition file to
 * conserved variables that can be used during the integration
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
    struct part* p, struct xpart* xp, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props) {}

/**
 * @brief Extra operations to be done during the drift
 *
 * Not used for Shadowswift.
 *
 * @param p Particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param dt The drift time-step.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part* p, struct xpart* xp, float dt_drift, float dt_therm) {}

/**
 * @brief Set the particle acceleration after the flux loop.
 *
 * @param p Particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_end_force(
    struct part* p, const struct cosmology* cosmo) {}

/**
 * @brief Extra operations done during the kick
 *
 * Not used for Shadowswift.
 *
 * @param p Particle to act upon.
 * @param xp Extended particle data to act upon.
 * @param dt Physical time step.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part* p, struct xpart* xp, float dt, float dt_grav, float dt_hydro,
    float dt_kick_corr, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props) {

  /* Update the conserved variables. We do this here and not in the kick,
     since we need the updated variables below. */
  p->conserved.mass += p->conserved.flux.mass * dt;
  p->conserved.momentum[0] += p->conserved.flux.momentum[0] * dt;
  p->conserved.momentum[1] += p->conserved.flux.momentum[1] * dt;
  p->conserved.momentum[2] += p->conserved.flux.momentum[2] * dt;

#ifdef EOS_ISOTHERMAL_GAS
  /* reset the thermal energy */
  p->conserved.energy =
      p->conserved.mass * gas_internal_energy_from_entropy(0.f, 0.f);
#else
  p->conserved.energy += p->conserved.flux.energy * dt;
#endif

#if defined(SHADOWFAX_FIX_CELLS)
  p->v[0] = 0.0f;
  p->v[1] = 0.0f;
  p->v[2] = 0.0f;
#else
  if (p->conserved.mass > 0.0f && p->primitives.rho > 0.0f) {

    const float inverse_mass = 1.f / p->conserved.mass;

    /* Normal case: set particle velocity to fluid velocity. */
    p->v[0] = p->conserved.momentum[0] * inverse_mass;
    p->v[1] = p->conserved.momentum[1] * inverse_mass;
    p->v[2] = p->conserved.momentum[2] * inverse_mass;

#ifdef SHADOWFAX_STEER_CELL_MOTION
    float centroid[3], d[3];
    float volume, csnd, R, vfac, fac, dnrm;
    voronoi_get_centroid(&p->cell, centroid);
    d[0] = centroid[0] - p->x[0];
    d[1] = centroid[1] - p->x[1];
    d[2] = centroid[2] - p->x[2];
    dnrm = sqrtf(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    csnd = sqrtf(hydro_gamma * p->primitives.P / p->primitives.rho);
    volume = p->cell.volume;
    R = get_radius_dimension_sphere(volume);
    fac = 4.0f * dnrm / R;
    if (fac > 0.9f) {
      if (fac < 1.1f) {
        vfac = csnd * (dnrm - 0.225f * R) / dnrm / (0.05f * R);
      } else {
        vfac = csnd / dnrm;
      }
      p->v[0] += vfac * d[0];
      p->v[1] += vfac * d[1];
      p->v[2] += vfac * d[2];
    }
#endif

  } else {
    p->v[0] = 0.;
    p->v[1] = 0.;
    p->v[2] = 0.;
  }
#endif

  /* Now make sure all velocity variables are up to date. */
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];

  if (p->gpart) {
    p->gpart->v_full[0] = p->v[0];
    p->gpart->v_full[1] = p->v[1];
    p->gpart->v_full[2] = p->v[2];
  }
}

/**
 * @brief Returns the internal energy of a particle
 *
 * @param p The particle of interest.
 * @return Internal energy of the particle.
 */
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy(
    const struct part* restrict p) {

  if (p->primitives.rho > 0.) {
    return gas_internal_energy_from_pressure(p->primitives.rho,
                                             p->primitives.P);
  } else {
    return 0.;
  }
}

/**
 * @brief Returns the entropy of a particle
 *
 * @param p The particle of interest.
 * @return Entropy of the particle.
 */
__attribute__((always_inline)) INLINE static float hydro_get_entropy(
    const struct part* restrict p) {

  if (p->primitives.rho > 0.) {
    return gas_entropy_from_pressure(p->primitives.rho, p->primitives.P);
  } else {
    return 0.;
  }
}

/**
 * @brief Returns the sound speed of a particle
 *
 * @param p The particle of interest.
 * @param Sound speed of the particle.
 */
__attribute__((always_inline)) INLINE static float hydro_get_soundspeed(
    const struct part* restrict p) {

  if (p->primitives.rho > 0.) {
    return gas_soundspeed_from_pressure(p->primitives.rho, p->primitives.P);
  } else {
    return 0.;
  }
}

/**
 * @brief Returns the pressure of a particle
 *
 * @param p The particle of interest
 * @param Pressure of the particle.
 */
__attribute__((always_inline)) INLINE static float hydro_get_pressure(
    const struct part* restrict p) {

  return p->primitives.P;
}

/**
 * @brief Returns the mass of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_mass(
    const struct part* restrict p) {

  return p->conserved.mass;
}

/**
 * @brief Returns the velocities drifted to the current time of a particle.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle.
 * @param dt The time since the last kick.
 * @param v (return) The velocities at the current time.
 */
__attribute__((always_inline)) INLINE static void hydro_get_drifted_velocities(
    const struct part* restrict p, const struct xpart* xp, float dt_kick_hydro,
    float dt_kick_grav, float v[3]) {

  v[0] = p->v[0];
  v[1] = p->v[1];
  v[2] = p->v[2];
}

/**
 * @brief Returns the density of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_density(
    const struct part* restrict p) {

  return p->primitives.rho;
}

/**
 * @brief Modifies the thermal state of a particle to the imposed internal
 * energy
 *
 * This overrides the current state of the particle but does *not* change its
 * time-derivatives
 *
 * @param p The particle
 * @param u The new internal energy
 */
__attribute__((always_inline)) INLINE static void hydro_set_internal_energy(
    struct part* restrict p, float u) {

  if (p->primitives.rho > 0.) {
    p->conserved.energy = u * p->conserved.mass;

#ifdef SHADOWFAX_TOTAL_ENERGY
    p->conserved.energy +=
        0.5f * (p->conserved.momentum[0] * p->primitives.v[0] +
                p->conserved.momentum[1] * p->primitives.v[1] +
                p->conserved.momentum[2] * p->primitives.v[2]);
#endif

    p->primitives.P = gas_pressure_from_internal_energy(p->primitives.rho, u);
  }
}

/**
 * @brief Modifies the thermal state of a particle to the imposed entropy
 *
 * This overrides the current state of the particle but does *not* change its
 * time-derivatives
 *
 * @param p The particle
 * @param S The new entropy
 */
__attribute__((always_inline)) INLINE static void hydro_set_entropy(
    struct part* restrict p, float S) {

  if (p->primitives.rho > 0.) {
    p->conserved.energy =
        gas_internal_energy_from_entropy(p->primitives.rho, S) *
        p->conserved.mass;

#ifdef SHADOWFAX_TOTAL_ENERGY
    p->conserved.energy +=
        0.5f * (p->conserved.momentum[0] * p->primitives.v[0] +
                p->conserved.momentum[1] * p->primitives.v[1] +
                p->conserved.momentum[2] * p->primitives.v[2]);
#endif

    p->primitives.P = gas_pressure_from_entropy(p->primitives.rho, S);
  }
}

/**
 * @brief Sets the mass of a particle
 *
 * @param p The particle of interest
 * @param m The mass to set.
 */
__attribute__((always_inline)) INLINE static void hydro_set_mass(
    struct part* restrict p, float m) {

  p->conserved.mass = m;
}

/**
 * @brief Overwrite the initial internal energy of a particle.
 *
 * Note that in the cases where the thermodynamic variable is not
 * internal energy but gets converted later, we must overwrite that
 * field. The conversion to the actual variable happens later after
 * the initial fake time-step.
 *
 * @param p The #part to write to.
 * @param u_init The new initial internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_init_internal_energy(struct part* p, float u_init) {

  p->conserved.energy = u_init * p->conserved.mass;
#ifdef GIZMO_TOTAL_ENERGY
  /* add the kinetic energy */
  p->conserved.energy += 0.5f * p->conserved.mass *
                         (p->conserved.momentum[0] * p->primitives.v[0] +
                          p->conserved.momentum[1] * p->primitives.v[1] +
                          p->conserved.momentum[2] * p->primitives.v[2]);
#endif
  p->primitives.P = hydro_gamma_minus_one * p->primitives.rho * u_init;
}

/**
 * @brief Returns the comoving internal energy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy(const struct part* restrict p) {

  if (p->primitives.rho > 0.)
    return gas_internal_energy_from_pressure(p->primitives.rho,
                                             p->primitives.P);
  else
    return 0.;
}

/**
 * @brief Returns the comoving entropy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_entropy(
    const struct part* restrict p) {

  if (p->primitives.rho > 0.) {
    return gas_entropy_from_pressure(p->primitives.rho, p->primitives.P);
  } else {
    return 0.;
  }
}

/**
 * @brief Returns the sound speed of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_soundspeed(const struct part* restrict p) {

  if (p->primitives.rho > 0.)
    return gas_soundspeed_from_pressure(p->primitives.rho, p->primitives.P);
  else
    return 0.;
}

/**
 * @brief Returns the comoving pressure of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_pressure(
    const struct part* restrict p) {

  return p->primitives.P;
}

/**
 * @brief Returns the comoving density of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_density(
    const struct part* restrict p) {

  return p->primitives.rho;
}

/**
 * @brief Returns the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy(const struct part* restrict p,
                                   const struct cosmology* cosmo) {

  return cosmo->a_factor_internal_energy *
         hydro_get_comoving_internal_energy(p);
}

/**
 * @brief Returns the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_entropy(
    const struct part* restrict p, const struct cosmology* cosmo) {

  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return hydro_get_comoving_entropy(p);
}

/**
 * @brief Returns the physical sound speed of a particle
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_soundspeed(const struct part* restrict p,
                              const struct cosmology* cosmo) {

  return cosmo->a_factor_sound_speed * hydro_get_comoving_soundspeed(p);
}

/**
 * @brief Sets the physical entropy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param entropy The physical entropy
 */
__attribute__((always_inline)) INLINE static void hydro_set_physical_entropy(
    struct part* p, struct xpart* xp, const struct cosmology* cosmo,
    const float entropy) {

  error("Needs implementing");
}

/**
 * @brief Sets the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy(struct part* p, struct xpart* xp,
                                   const struct cosmology* cosmo,
                                   const float u) {
  error("Need implementing");
}

/**
 * @brief Sets the drifted physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(struct part* p,
                                           const struct cosmology* cosmo,
                                           const float u) {
  error("Need implementing");
}

/**
 * @brief Returns the comoving pressure of a particle
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_pressure(
    const struct part* restrict p, const struct cosmology* cosmo) {

  return cosmo->a_factor_pressure * p->primitives.P;
}

/**
 * @brief Returns the physical density of a particle
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_density(
    const struct part* restrict p, const struct cosmology* cosmo) {

  return cosmo->a3_inv * p->primitives.rho;
}

/**
 * @brief Operations performed when a particle gets removed from the
 * simulation volume.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
__attribute__((always_inline)) INLINE static void hydro_remove_part(
    const struct part* p, const struct xpart* xp) {}

#endif /* SWIFT_SHADOWSWIFT_HYDRO_H */
