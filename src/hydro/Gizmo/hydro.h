
/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2016, 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#include <float.h>
#include "adiabatic_index.h"
#include "approx_math.h"
#include "equation_of_state.h"
#include "hydro_gradients.h"
#include "hydro_space.h"
#include "hydro_unphysical.h"
#include "hydro_velocities.h"
#include "minmax.h"
#include "riemann.h"

//#define GIZMO_LLOYD_ITERATION

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param hydro_properties Pointer to the hydro parameters.
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    const struct part* restrict p, const struct xpart* restrict xp,
    const struct hydro_props* restrict hydro_properties) {

  const float CFL_condition = hydro_properties->CFL_condition;

#ifdef GIZMO_LLOYD_ITERATION
  return CFL_condition;
#endif

  float vrel[3];
  vrel[0] = p->primitives.v[0] - xp->v_full[0];
  vrel[1] = p->primitives.v[1] - xp->v_full[1];
  vrel[2] = p->primitives.v[2] - xp->v_full[2];
  float vmax =
      sqrtf(vrel[0] * vrel[0] + vrel[1] * vrel[1] + vrel[2] * vrel[2]) +
      sqrtf(hydro_gamma * p->primitives.P / p->primitives.rho);
  vmax = max(vmax, p->timestepvars.vmax);
  const float psize = powf(p->geometry.volume / hydro_dimension_unit_sphere,
                           hydro_dimension_inv);
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

#ifdef SWIFT_DEBUG_CHECKS
  if (dt == 0.) {
    error("Zero time step assigned to particle!");
  }

  if (dt != dt) {
    error("NaN time step assigned to particle!");
  }
#endif

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

  p->primitives.v[0] = p->v[0];
  p->primitives.v[1] = p->v[1];
  p->primitives.v[2] = p->v[2];

  /* we can already initialize the momentum */
  p->conserved.momentum[0] = mass * p->primitives.v[0];
  p->conserved.momentum[1] = mass * p->primitives.v[1];
  p->conserved.momentum[2] = mass * p->primitives.v[2];

/* and the thermal energy */
/* remember that we store the total thermal energy, not the specific thermal
   energy (as in Gadget) */
#if defined(EOS_ISOTHERMAL_GAS)
  /* this overwrites the internal energy from the initial condition file */
  p->conserved.energy = mass * const_isothermal_internal_energy;
#else
  p->conserved.energy *= mass;
#endif

#ifdef GIZMO_TOTAL_ENERGY
  /* add the total kinetic energy */
  p->conserved.energy += 0.5f * (p->conserved.momentum[0] * p->primitives.v[0] +
                                 p->conserved.momentum[1] * p->primitives.v[1] +
                                 p->conserved.momentum[2] * p->primitives.v[2]);
#endif

#ifdef GIZMO_LLOYD_ITERATION
  /* overwrite all variables to make sure they have safe values */
  p->primitives.rho = 1.;
  p->primitives.v[0] = 0.;
  p->primitives.v[1] = 0.;
  p->primitives.v[2] = 0.;
  p->primitives.P = 1.;

  p->conserved.mass = 1.;
  p->conserved.momentum[0] = 0.;
  p->conserved.momentum[1] = 0.;
  p->conserved.momentum[2] = 0.;
  p->conserved.energy = 1.;

  p->v[0] = 0.;
  p->v[1] = 0.;
  p->v[2] = 0.;
#endif

  /* initialize the particle velocity based on the primitive fluid velocity */
  hydro_velocities_init(p, xp);

  /* we cannot initialize wcorr in init_part, as init_part gets called every
     time the density loop is repeated, and the whole point of storing wcorr
     is to have a way of remembering that we need more neighbours for this
     particle */
  p->density.wcorr = 1.0f;
}

/**
 * @brief Prepares a particle for the volume calculation.
 *
 * Simply makes sure all necessary variables are initialized to zero.
 *
 * @param p The particle to act upon
 * @param hs #hydro_space containing hydro specific space information.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part* p, const struct hydro_space* hs) {

  p->density.wcount = 0.0f;
  p->density.wcount_dh = 0.0f;
  p->geometry.volume = 0.0f;
  p->geometry.matrix_E[0][0] = 0.0f;
  p->geometry.matrix_E[0][1] = 0.0f;
  p->geometry.matrix_E[0][2] = 0.0f;
  p->geometry.matrix_E[1][0] = 0.0f;
  p->geometry.matrix_E[1][1] = 0.0f;
  p->geometry.matrix_E[1][2] = 0.0f;
  p->geometry.matrix_E[2][0] = 0.0f;
  p->geometry.matrix_E[2][1] = 0.0f;
  p->geometry.matrix_E[2][2] = 0.0f;
  p->geometry.centroid[0] = 0.0f;
  p->geometry.centroid[1] = 0.0f;
  p->geometry.centroid[2] = 0.0f;
  p->geometry.Atot = 0.0f;

  /* Set the active flag to active. */
  p->force.active = 1;
}

/**
 * @brief Finishes the volume calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and adds the self-contribution term. Calculates the volume and uses it to
 * update the primitive variables (based on the conserved variables). The latter
 * should only be done for active particles. This is okay, since this method is
 * only called for active particles.
 *
 * Multiplies the components of the matrix E with the appropriate constants and
 * inverts it. Initializes the variables used during the gradient loop. This
 * cannot be done in hydro_prepare_force, since that method is called for all
 * particles, and not just the active ones. If we would initialize the
 * variables there, gradients for passive particles would be zero, while we
 * actually use the old gradients in the flux calculation between active and
 * passive particles.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part* restrict p) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float ih = 1.0f / h;
  const float ihdim = pow_dimension(ih);
  const float ihdim_plus_one = ihdim * ih;

  /* Final operation on the density. */
  p->density.wcount += kernel_root;
  p->density.wcount *= ihdim;

  p->density.wcount_dh -= hydro_dimension * kernel_root;
  p->density.wcount_dh *= ihdim_plus_one;

  /* Final operation on the geometry. */
  /* we multiply with the smoothing kernel normalization ih3 and calculate the
   * volume */
  const float volume = 1.f / (ihdim * (p->geometry.volume + kernel_root));
  p->geometry.volume = volume;

  /* we multiply with the smoothing kernel normalization */
  p->geometry.matrix_E[0][0] = ihdim * p->geometry.matrix_E[0][0];
  p->geometry.matrix_E[0][1] = ihdim * p->geometry.matrix_E[0][1];
  p->geometry.matrix_E[0][2] = ihdim * p->geometry.matrix_E[0][2];
  p->geometry.matrix_E[1][0] = ihdim * p->geometry.matrix_E[1][0];
  p->geometry.matrix_E[1][1] = ihdim * p->geometry.matrix_E[1][1];
  p->geometry.matrix_E[1][2] = ihdim * p->geometry.matrix_E[1][2];
  p->geometry.matrix_E[2][0] = ihdim * p->geometry.matrix_E[2][0];
  p->geometry.matrix_E[2][1] = ihdim * p->geometry.matrix_E[2][1];
  p->geometry.matrix_E[2][2] = ihdim * p->geometry.matrix_E[2][2];

  p->geometry.centroid[0] *= kernel_norm;
  p->geometry.centroid[1] *= kernel_norm;
  p->geometry.centroid[2] *= kernel_norm;

  p->geometry.centroid[0] /= p->density.wcount;
  p->geometry.centroid[1] /= p->density.wcount;
  p->geometry.centroid[2] /= p->density.wcount;

  /* Check the condition number to see if we have a stable geometry. */
  float condition_number_E = 0.0f;
  int i, j;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      condition_number_E +=
          p->geometry.matrix_E[i][j] * p->geometry.matrix_E[i][j];
    }
  }

  invert_dimension_by_dimension_matrix(p->geometry.matrix_E);

  float condition_number_Einv = 0.0f;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      condition_number_Einv +=
          p->geometry.matrix_E[i][j] * p->geometry.matrix_E[i][j];
    }
  }

  float condition_number =
      hydro_dimension_inv * sqrtf(condition_number_E * condition_number_Einv);

  if (condition_number > const_gizmo_max_condition_number &&
      p->density.wcorr > const_gizmo_min_wcorr) {
#ifdef GIZMO_PATHOLOGICAL_ERROR
    error("Condition number larger than %g (%g)!",
          const_gizmo_max_condition_number, condition_number);
#endif
#ifdef GIZMO_PATHOLOGICAL_WARNING
    message("Condition number too large: %g (> %g, p->id: %llu)!",
            condition_number, const_gizmo_max_condition_number, p->id);
#endif
    /* add a correction to the number of neighbours for this particle */
    p->density.wcorr *= const_gizmo_w_correction_factor;
  }

  hydro_gradients_init(p);

  /* compute primitive variables */
  /* eqns (3)-(5) */
  const float m = p->conserved.mass;

#ifdef SWIFT_DEBUG_CHECKS
  if (m < 0.) {
    error("Mass is negative!");
  }

  if (volume == 0.) {
    error("Volume is 0!");
  }
#endif

  float momentum[3];
  momentum[0] = p->conserved.momentum[0];
  momentum[1] = p->conserved.momentum[1];
  momentum[2] = p->conserved.momentum[2];
  p->primitives.rho = m / volume;
  if (m == 0.) {
    p->primitives.v[0] = 0.;
    p->primitives.v[1] = 0.;
    p->primitives.v[2] = 0.;
  } else {
    p->primitives.v[0] = momentum[0] / m;
    p->primitives.v[1] = momentum[1] / m;
    p->primitives.v[2] = momentum[2] / m;
  }

#ifdef EOS_ISOTHERMAL_GAS
  /* although the pressure is not formally used anywhere if an isothermal eos
     has been selected, we still make sure it is set to the correct value */
  p->primitives.P = gas_pressure_from_internal_energy(p->primitives.rho, 0.);
#else

  float energy = p->conserved.energy;

#ifdef GIZMO_TOTAL_ENERGY
  /* subtract the kinetic energy; we want the thermal energy */
  energy -= 0.5f * (momentum[0] * p->primitives.v[0] +
                    momentum[1] * p->primitives.v[1] +
                    momentum[2] * p->primitives.v[2]);
#endif

  /* energy contains the total thermal energy, we want the specific energy.
     this is why we divide by the volume, and not by the density */
  p->primitives.P = hydro_gamma_minus_one * energy / volume;
#endif

  /* sanity checks */
  gizmo_check_physical_quantity("density", p->primitives.rho);
  gizmo_check_physical_quantity("pressure", p->primitives.P);

#ifdef GIZMO_LLOYD_ITERATION
  /* overwrite primitive variables to make sure they still have safe values */
  p->primitives.rho = 1.;
  p->primitives.v[0] = 0.;
  p->primitives.v[1] = 0.;
  p->primitives.v[2] = 0.;
  p->primitives.P = 1.;
#endif

  /* Add a correction factor to wcount (to force a neighbour number increase if
     the geometry matrix is close to singular) */
  p->density.wcount *= p->density.wcorr;
  p->density.wcount_dh *= p->density.wcorr;
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_part_has_no_neighbours(
    struct part* restrict p, struct xpart* restrict xp) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  p->density.wcount = kernel_root * kernel_norm * h_inv_dim;
  p->density.wcount_dh = 0.f;
  p->geometry.volume = 1.0f;
  p->geometry.matrix_E[0][0] = 1.0f;
  p->geometry.matrix_E[0][1] = 0.0f;
  p->geometry.matrix_E[0][2] = 0.0f;
  p->geometry.matrix_E[1][0] = 0.0f;
  p->geometry.matrix_E[1][1] = 1.0f;
  p->geometry.matrix_E[1][2] = 0.0f;
  p->geometry.matrix_E[2][0] = 0.0f;
  p->geometry.matrix_E[2][1] = 0.0f;
  p->geometry.matrix_E[2][2] = 1.0f;
  /* centroid is relative w.r.t. particle position */
  /* by setting the centroid to 0.0f, we make sure no velocity correction is
     applied */
  p->geometry.centroid[0] = 0.0f;
  p->geometry.centroid[1] = 0.0f;
  p->geometry.centroid[2] = 0.0f;
  p->geometry.Atot = 1.0f;
}

/**
 * @brief Prepare a particle for the gradient calculation.
 *
 * The name of this method is confusing, as this method is really called after
 * the density loop and before the gradient loop.
 *
 * We use it to set the physical timestep for the particle and to copy the
 * actual velocities, which we need to boost our interfaces during the flux
 * calculation. We also initialize the variables used for the time step
 * calculation.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part* restrict p, struct xpart* restrict xp) {

  /* Initialize time step criterion variables */
  p->timestepvars.vmax = 0.;

  /* Set the actual velocity of the particle */
  hydro_velocities_prepare_force(p, xp);
}

/**
 * @brief Finishes the gradient calculation.
 *
 * Just a wrapper around hydro_gradients_finalize, which can be an empty method,
 * in which case no gradients are used.
 *
 * This method also initializes the force loop variables.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_end_gradient(
    struct part* p) {

  hydro_gradients_finalize(p);

  p->gravity.mflux[0] = 0.0f;
  p->gravity.mflux[1] = 0.0f;
  p->gravity.mflux[2] = 0.0f;

#ifdef GIZMO_LLOYD_ITERATION
  /* reset the gradients to zero, as we don't want them */
  hydro_gradients_init(p);
#endif
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is actually not necessary for GIZMO, since we just set the accelerations
 * after the flux calculation.
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
 * We no longer do this, as the mass needs to be provided in the initial
 * condition file, and the mass alone is enough to initialize all conserved
 * variables. This is now done in hydro_first_init_part.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part* p, struct xpart* xp) {}

/**
 * @brief Extra operations to be done during the drift
 *
 * @param p Particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param dt The drift time-step.
 * @param t0 Integer start time of the drift interval.
 * @param t1 Integer end time of the drift interval.
 * @param timeBase Conversion factor between integer and physical time.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part* p, struct xpart* xp, float dt) {

#ifdef GIZMO_LLOYD_ITERATION
  return;
#endif

  const float h_inv = 1.0f / p->h;

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt;
  float h_corr;
  if (fabsf(w1) < 0.2f)
    h_corr = approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    h_corr = expf(w1);

  /* Limit the smoothing length correction (and make sure it is always
     positive). */
  if (h_corr < 2.0f && h_corr > 0.) {
    p->h *= h_corr;
  }

/* we temporarily disabled the primitive variable drift.
   This should be reenabled once gravity works, and we have time to check that
   drifting works properly. */
//  const float w2 = -hydro_dimension * w1;
//  if (fabsf(w2) < 0.2f) {
//    p->primitives.rho *= approx_expf(w2);
//  } else {
//    p->primitives.rho *= expf(w2);
//  }

//  p->primitives.v[0] += (p->a_hydro[0] + p->gravity.old_a[0]) * dt;
//  p->primitives.v[1] += (p->a_hydro[1] + p->gravity.old_a[1]) * dt;
//  p->primitives.v[2] += (p->a_hydro[2] + p->gravity.old_a[2]) * dt;

//#if !defined(EOS_ISOTHERMAL_GAS)
//  if (p->conserved.mass > 0.) {
//    const float u = p->conserved.energy + p->du_dt * dt;
//    p->primitives.P =
//        hydro_gamma_minus_one * u * p->primitives.rho / p->conserved.mass;
//  }
//#endif

#ifdef SWIFT_DEBUG_CHECKS
  if (p->h <= 0.) {
    error("Zero or negative smoothing length (%g)!", p->h);
  }
#endif
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
 */
__attribute__((always_inline)) INLINE static void hydro_end_force(
    struct part* p) {

  /* set the variables that are used to drift the primitive variables */

  if (p->force.dt > 0.) {
    p->du_dt = p->conserved.flux.energy / p->force.dt;
  } else {
    p->du_dt = 0.0f;
  }

  hydro_velocities_end_force(p);
}

/**
 * @brief Extra operations done during the kick
 *
 * Not used for GIZMO.
 *
 * @param p Particle to act upon.
 * @param xp Extended particle data to act upon.
 * @param dt Physical time step.
 * @param half_dt Half the physical time step.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part* p, struct xpart* xp, float dt) {

  float a_grav[3];

  /* Update conserved variables. */
  p->conserved.mass += p->conserved.flux.mass;
  p->conserved.momentum[0] += p->conserved.flux.momentum[0];
  p->conserved.momentum[1] += p->conserved.flux.momentum[1];
  p->conserved.momentum[2] += p->conserved.flux.momentum[2];
#if defined(EOS_ISOTHERMAL_GAS)
  p->conserved.energy = p->conserved.mass * const_isothermal_internal_energy;
#else
  p->conserved.energy += p->conserved.flux.energy;
#endif

  gizmo_check_physical_quantity("mass", p->conserved.mass);
  gizmo_check_physical_quantity("energy", p->conserved.energy);

#ifdef SWIFT_DEBUG_CHECKS
  /* Note that this check will only have effect if no GIZMO_UNPHYSICAL option
     was selected. */
  if (p->conserved.mass < 0.) {
    error(
        "Negative mass after conserved variables update (mass: %g, dmass: %g)!",
        p->conserved.mass, p->conserved.flux.mass);
  }

  if (p->conserved.energy < 0.) {
    error(
        "Negative energy after conserved variables update (energy: %g, "
        "denergy: %g)!",
        p->conserved.energy, p->conserved.flux.energy);
  }
#endif

  /* Add gravity. We only do this if we have gravity activated. */
  if (p->gpart) {
    /* Retrieve the current value of the gravitational acceleration from the
       gpart. We are only allowed to do this because this is the kick. We still
       need to check whether gpart exists though.*/
    a_grav[0] = p->gpart->a_grav[0];
    a_grav[1] = p->gpart->a_grav[1];
    a_grav[2] = p->gpart->a_grav[2];

    /* Make sure the gpart knows the mass has changed. */
    p->gpart->mass = p->conserved.mass;

#if !defined(EOS_ISOTHERMAL_GAS)
    /* If the energy needs to be updated, we need to do it before the momentum
       is updated, as the old value of the momentum enters the equations. */
    p->conserved.energy += dt * (p->conserved.momentum[0] * a_grav[0] +
                                 p->conserved.momentum[1] * a_grav[1] +
                                 p->conserved.momentum[2] * a_grav[2]);

    p->conserved.energy += dt * (a_grav[0] * p->gravity.mflux[0] +
                                 a_grav[1] * p->gravity.mflux[1] +
                                 a_grav[2] * p->gravity.mflux[2]);
#endif

    /* Kick the momentum for half a time step */
    /* Note that this also affects the particle movement, as the velocity for
       the particles is set after this. */
    p->conserved.momentum[0] += dt * p->conserved.mass * a_grav[0];
    p->conserved.momentum[1] += dt * p->conserved.mass * a_grav[1];
    p->conserved.momentum[2] += dt * p->conserved.mass * a_grav[2];
  }

  /* reset fluxes */
  /* we can only do this here, since we need to keep the fluxes for inactive
     particles */
  p->conserved.flux.mass = 0.0f;
  p->conserved.flux.momentum[0] = 0.0f;
  p->conserved.flux.momentum[1] = 0.0f;
  p->conserved.flux.momentum[2] = 0.0f;
  p->conserved.flux.energy = 0.0f;

  hydro_velocities_set(p, xp);

#ifdef GIZMO_LLOYD_ITERATION
  /* reset conserved variables to safe values */
  p->conserved.mass = 1.;
  p->conserved.momentum[0] = 0.;
  p->conserved.momentum[1] = 0.;
  p->conserved.momentum[2] = 0.;
  p->conserved.energy = 1.;

  /* set the particle velocities to the Lloyd velocities */
  /* note that centroid is the relative position of the centroid w.r.t. the
     particle position (position - centroid) */
  xp->v_full[0] = -p->geometry.centroid[0] / p->force.dt;
  xp->v_full[1] = -p->geometry.centroid[1] / p->force.dt;
  xp->v_full[2] = -p->geometry.centroid[2] / p->force.dt;
  p->v[0] = xp->v_full[0];
  p->v[1] = xp->v_full[1];
  p->v[2] = xp->v_full[2];
#endif

  /* reset wcorr */
  p->density.wcorr = 1.0f;
}

/**
 * @brief Returns the internal energy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy(
    const struct part* restrict p) {

  if (p->primitives.rho > 0.) {
#ifdef EOS_ISOTHERMAL_GAS
    return p->primitives.P / hydro_gamma_minus_one / p->primitives.rho;
#else
    return const_isothermal_internal_energy;
#endif
  } else {
    return 0.;
  }
}

/**
 * @brief Returns the entropy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_entropy(
    const struct part* restrict p) {

  if (p->primitives.rho > 0.) {
    return p->primitives.P / pow_gamma(p->primitives.rho);
  } else {
    return 0.;
  }
}

/**
 * @brief Returns the sound speed of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_soundspeed(
    const struct part* restrict p) {

  if (p->primitives.rho > 0.) {
#ifdef EOS_ISOTHERMAL_GAS
    return sqrtf(const_isothermal_internal_energy * hydro_gamma *
                 hydro_gamma_minus_one);
#else
    return sqrtf(hydro_gamma * p->primitives.P / p->primitives.rho);
#endif
  } else {
    return 0.;
  }
}

/**
 * @brief Returns the pressure of a particle
 *
 * @param p The particle of interest
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

  /* conserved.energy is NOT the specific energy (u), but the total thermal
     energy (u*m) */
  p->conserved.energy = u * p->conserved.mass;
#ifdef GIZMO_TOTAL_ENERGY
  /* add the kinetic energy */
  p->conserved.energy += 0.5f * p->conserved.mass *
                         (p->conserved.momentum[0] * p->primitives.v[0] +
                          p->conserved.momentum[1] * p->primitives.v[1] +
                          p->conserved.momentum[2] * p->primitives.v[2]);
#endif
  p->primitives.P = hydro_gamma_minus_one * p->primitives.rho * u;
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

  p->conserved.energy = S * pow_gamma_minus_one(p->primitives.rho) *
                        hydro_one_over_gamma_minus_one * p->conserved.mass;
#ifdef GIZMO_TOTAL_ENERGY
  /* add the kinetic energy */
  p->conserved.energy += 0.5f * p->conserved.mass *
                         (p->conserved.momentum[0] * p->primitives.v[0] +
                          p->conserved.momentum[1] * p->primitives.v[1] +
                          p->conserved.momentum[2] * p->primitives.v[2]);
#endif
  p->primitives.P = S * pow_gamma(p->primitives.rho);
}
