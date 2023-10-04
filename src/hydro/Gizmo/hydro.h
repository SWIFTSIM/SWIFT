/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_HYDRO_H
#define SWIFT_GIZMO_HYDRO_H

/**
 * @brief Enable Lloyd's iteration.
 *
 * If you enable the flag below, the code will ignore all hydrodynamical
 * variables and instead run in a mode where
 */
/*#define GIZMO_LLOYD_ITERATION*/

#include "approx_math.h"
#include "entropy_floor.h"
#include "hydro_flux.h"
#include "hydro_getters.h"
#include "hydro_gradients.h"
#include "hydro_lloyd.h"
#include "hydro_parameters.h"
#include "hydro_properties.h"
#include "hydro_setters.h"
#include "hydro_space.h"
#include "hydro_velocities.h"
#include "pressure_floor.h"

#include <float.h>

#if defined(GIZMO_MFV_SPH)
#define SPH_IMPLEMENTATION "GIZMO MFV (Hopkins 2015)"
#elif defined(GIZMO_MFM_SPH)
#define SPH_IMPLEMENTATION "GIZMO MFM (Hopkins 2015)"
#endif

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param hydro_properties Pointer to the hydro parameters.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    const struct part* restrict p, const struct xpart* restrict xp,
    const struct hydro_props* restrict hydro_properties,
    const struct cosmology* restrict cosmo) {

  const float CFL_condition = hydro_properties->CFL_condition;

  /* skip the time step calculation if we are using Lloyd's algorithm */
  hydro_gizmo_lloyd_skip_timestep(CFL_condition);

  float W[5];
  hydro_part_get_primitive_variables(p, W);

  /* The time step depends on the relative difference of the fluid velocity and
   * the particle velocity. */
  float v_rel[3];
  hydro_part_get_relative_fluid_velocity(p, v_rel);
  const float rhoinv = (W[0] > 0.0f) ? 1.0f / W[0] : 0.0f;
  float vmax =
      sqrtf(v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] + v_rel[2] * v_rel[2]) +
      sqrtf(hydro_gamma * W[4] * rhoinv);
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
 * @brief Does some extra hydro operations once the actual physical time step
 * for the particle is known.
 *
 * This method is no longer used, as Gizmo is now unaware of the actual particle
 * time step.
 *
 * @param p The particle to act upon.
 * @param dt Physical time step of the particle during the next step.
 */
__attribute__((always_inline)) INLINE static void hydro_timestep_extra(
    struct part* p, float dt) {

#ifdef SWIFT_DEBUG_CHECKS
  if (dt == 0.0f) {
    error("Zero time step assigned to particle!");
  }

  if (dt != dt) {
    error("NaN time step assigned to particle!");
  }
#endif
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

#ifdef GIZMO_TOTAL_ENERGY
  Q[4] += 0.5f * (Q[1] * W[1] + Q[2] * W[2] + Q[3] * W[3]);
#endif

  /* overwrite all hydro variables if we are using Lloyd's algorithm */
  hydro_gizmo_lloyd_initialize_particle(W, Q, p->v);

  p->time_bin = 0;

  hydro_part_set_primitive_variables(p, W);
  hydro_part_set_conserved_variables(p, Q);

  /* initialize the particle velocity based on the primitive fluid velocity */
  hydro_velocities_init(p, xp);

  /* ignore accelerations present in the initial condition */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* we cannot initialize wcorr in init_part, as init_part gets called every
     time the density loop is repeated, and the whole point of storing wcorr
     is to have a way of remembering that we need more neighbours for this
     particle */
  p->geometry.wcorr = 1.0f;
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

  /* reset the centroid variables used for the velocity correction in MFV */
  hydro_velocities_reset_centroids(p);
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
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part* restrict p, const struct cosmology* cosmo) {

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
  const float volume_inv = ihdim * (p->geometry.volume + kernel_root);
  const float volume = 1.0f / volume_inv;
  p->geometry.volume = volume;

  /* we multiply with the smoothing kernel normalization */
  p->geometry.matrix_E[0][0] *= ihdim;
  p->geometry.matrix_E[0][1] *= ihdim;
  p->geometry.matrix_E[0][2] *= ihdim;
  p->geometry.matrix_E[1][0] *= ihdim;
  p->geometry.matrix_E[1][1] *= ihdim;
  p->geometry.matrix_E[1][2] *= ihdim;
  p->geometry.matrix_E[2][0] *= ihdim;
  p->geometry.matrix_E[2][1] *= ihdim;
  p->geometry.matrix_E[2][2] *= ihdim;

  /* normalise the centroids for MFV */
  hydro_velocities_normalise_centroid(p, p->density.wcount);

  /* Check the condition number to see if we have a stable geometry. */
  const float condition_number_E =
      p->geometry.matrix_E[0][0] * p->geometry.matrix_E[0][0] +
      p->geometry.matrix_E[0][1] * p->geometry.matrix_E[0][1] +
      p->geometry.matrix_E[0][2] * p->geometry.matrix_E[0][2] +
      p->geometry.matrix_E[1][0] * p->geometry.matrix_E[1][0] +
      p->geometry.matrix_E[1][1] * p->geometry.matrix_E[1][1] +
      p->geometry.matrix_E[1][2] * p->geometry.matrix_E[1][2] +
      p->geometry.matrix_E[2][0] * p->geometry.matrix_E[2][0] +
      p->geometry.matrix_E[2][1] * p->geometry.matrix_E[2][1] +
      p->geometry.matrix_E[2][2] * p->geometry.matrix_E[2][2];

  float condition_number = 0.0f;
  if (invert_dimension_by_dimension_matrix(p->geometry.matrix_E) != 0) {
    /* something went wrong in the inversion; force bad condition number */
    condition_number = const_gizmo_max_condition_number + 1.0f;
  } else {
    const float condition_number_Einv =
        p->geometry.matrix_E[0][0] * p->geometry.matrix_E[0][0] +
        p->geometry.matrix_E[0][1] * p->geometry.matrix_E[0][1] +
        p->geometry.matrix_E[0][2] * p->geometry.matrix_E[0][2] +
        p->geometry.matrix_E[1][0] * p->geometry.matrix_E[1][0] +
        p->geometry.matrix_E[1][1] * p->geometry.matrix_E[1][1] +
        p->geometry.matrix_E[1][2] * p->geometry.matrix_E[1][2] +
        p->geometry.matrix_E[2][0] * p->geometry.matrix_E[2][0] +
        p->geometry.matrix_E[2][1] * p->geometry.matrix_E[2][1] +
        p->geometry.matrix_E[2][2] * p->geometry.matrix_E[2][2];

    condition_number =
        hydro_dimension_inv * sqrtf(condition_number_E * condition_number_Einv);
  }

  if (condition_number > const_gizmo_max_condition_number &&
      p->geometry.wcorr > const_gizmo_min_wcorr) {
#ifdef GIZMO_PATHOLOGICAL_ERROR
    error("Condition number larger than %g (%g)!",
          const_gizmo_max_condition_number, condition_number);
#endif
#ifdef GIZMO_PATHOLOGICAL_WARNING
    message("Condition number too large: %g (> %g, p->id: %llu)!",
            condition_number, const_gizmo_max_condition_number, p->id);
#endif
    /* add a correction to the number of neighbours for this particle */
    p->geometry.wcorr = const_gizmo_w_correction_factor * p->geometry.wcorr;
  }

  /* compute primitive variables */
  /* eqns (3)-(5) */
  float Q[5] = {p->conserved.mass, p->conserved.momentum[0],
                p->conserved.momentum[1], p->conserved.momentum[2],
                p->conserved.energy};

#ifdef SWIFT_DEBUG_CHECKS
  if (Q[0] < 0.) {
    error("Mass is negative!");
  }

  if (volume == 0.) {
    error("Volume is 0!");
  }
#endif

  float W[5];

  W[0] = Q[0] * volume_inv;
  const float m_inv = (Q[0] != 0.0f) ? 1.0f / Q[0] : 0.0f;
  W[1] = Q[1] * m_inv;
  W[2] = Q[2] * m_inv;
  W[3] = Q[3] * m_inv;

#ifdef EOS_ISOTHERMAL_GAS
  /* although the pressure is not formally used anywhere if an isothermal eos
     has been selected, we still make sure it is set to the correct value */
  W[4] = gas_pressure_from_internal_energy(W[0], 0.0f);
#else

#ifdef GIZMO_TOTAL_ENERGY
  /* subtract the kinetic energy; we want the thermal energy */
  Q[4] -= 0.5f * (Q[1] * W[1] + Q[2] * W[2] + Q[3] * W[3]);
#endif

  /* energy contains the total thermal energy, we want the specific energy.
     this is why we divide by the volume, and not by the density */
  W[4] = hydro_gamma_minus_one * Q[4] * volume_inv;
#endif

  /* sanity checks */
  gizmo_check_physical_quantities("density", "pressure", W[0], W[1], W[2], W[3],
                                  W[4]);

  /* reset the primitive variables if we are using Lloyd's algorithm */
  hydro_gizmo_lloyd_reset_primitive_variables(W);

  hydro_part_set_primitive_variables(p, W);

  /* Add a correction factor to wcount (to force a neighbour number increase if
     the geometry matrix is close to singular) */
  p->density.wcount *= p->geometry.wcorr;
  p->density.wcount_dh *= p->geometry.wcorr;
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_part_has_no_neighbours(
    struct part* restrict p, struct xpart* restrict xp,
    const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  warning(
      "Gas particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      p->id, h, p->density.wcount);

  /* Re-set problematic values */
  p->density.wcount = kernel_root * h_inv_dim;
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

  /* reset the centroid to disable MFV velocity corrections for this particle */
  hydro_velocities_reset_centroids(p);
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
 * @param hydro_props Hydrodynamic properties.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_gradient(
    struct part* restrict p, struct xpart* restrict xp,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct pressure_floor_props* pressure_floor) {

  /* Initialize time step criterion variables */
  p->timestepvars.vmax = 0.;

  hydro_gradients_init(p);

  hydro_velocities_prepare_force(p, xp);
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
 * This method also initializes the force loop variables.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_end_gradient(
    struct part* p) {

  hydro_gradients_finalize(p);

  /* reset the gradients if we are using Lloyd's algorith; we don't use them */
  hydro_gizmo_lloyd_reset_gradients(p);
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
    struct part* restrict p, struct xpart* restrict xp,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct pressure_floor_props* pressure_floor, const float dt_alpha,
    const float dt_therm) {

  hydro_part_reset_gravity_fluxes(p);
  p->flux.dt = dt_therm;
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
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values(
    struct part* restrict p, const struct xpart* restrict xp,
    const struct cosmology* cosmo,
    const struct pressure_floor_props* pressure_floor) {
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
    struct part* p, struct xpart* xp, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props,
    const struct pressure_floor_props* pressure_floor) {

  p->conserved.energy /= cosmo->a_factor_internal_energy;
}

/**
 * @brief Extra operations to be done during the drift
 *
 * @param p Particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part* p, struct xpart* xp, float dt_drift, float dt_therm,
    float dt_kick_grav, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* floor_props,
    const struct pressure_floor_props* pressure_floor) {

  /* skip the drift if we are using Lloyd's algorithm */
  hydro_gizmo_lloyd_skip_drift();

  const float h_inv = 1.0f / p->h;

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt_drift;
  float h_corr;
  if (fabsf(w1) < 0.2f) {
    h_corr = approx_expf(w1); /* 4th order expansion of exp(w) */
  } else {
    h_corr = expf(w1);
  }

  /* Limit the smoothing length correction (and make sure it is always
     positive). */
  if (h_corr < 2.0f && h_corr > 0.0f) {
    p->h *= h_corr;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (p->h <= 0.) {
    error("Zero or negative smoothing length (%g)!", p->h);
  }
#endif

  /* Reset the particle velocity. (undo the drift) */
  hydro_set_particle_velocity(p, xp->v_full);

  float W[5];
  hydro_part_get_primitive_variables(p, W);

  /* Use the fluid velocity in the rest frame of the particle for the time
   * extrapolation to preserve Galilean invariance. */
  float v_rel[3];
  hydro_part_get_relative_fluid_velocity(p, v_rel);

  float gradrho[3], gradvx[3], gradvy[3], gradvz[3], gradP[3];
  hydro_part_get_gradients(p, gradrho, gradvx, gradvy, gradvz, gradP);

  const float divv = gradvx[0] + gradvy[1] + gradvz[2];

  float Wprime[5];
  Wprime[0] = W[0] - dt_therm * (W[0] * divv + v_rel[0] * gradrho[0] +
                                 v_rel[1] * gradrho[1] + v_rel[2] * gradrho[2]);
  if (W[0] != 0.0f) {
    const float rhoinv = 1.0f / W[0];
    Wprime[1] = W[1] - dt_therm * (v_rel[0] * divv + rhoinv * gradP[0]);
    Wprime[2] = W[2] - dt_therm * (v_rel[1] * divv + rhoinv * gradP[1]);
    Wprime[3] = W[3] - dt_therm * (v_rel[2] * divv + rhoinv * gradP[2]);
  } else {
    Wprime[1] = 0.0f;
    Wprime[2] = 0.0f;
    Wprime[3] = 0.0f;
  }
  Wprime[4] =
      W[4] - dt_therm * (hydro_gamma * W[4] * divv + v_rel[0] * gradP[0] +
                         v_rel[1] * gradP[1] + v_rel[2] * gradP[2]);

  W[0] = Wprime[0];
  W[1] = Wprime[1];
  W[2] = Wprime[2];
  W[3] = Wprime[3];
  W[4] = Wprime[4];

  // MATTHIEU: Apply the entropy floor here.

  /* add the gravitational contribution to the fluid velocity drift */
  /* (MFV only) */
  hydro_gizmo_mfv_extra_velocity_drift(p->v, p->fluid_v, xp->v_full, xp->a_grav,
                                       dt_kick_grav);

  gizmo_check_physical_quantities("density", "pressure", W[0], W[1], W[2], W[3],
                                  W[4]);

  hydro_part_set_primitive_variables(p, W);
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
    struct part* p, const struct cosmology* cosmo) {

  hydro_velocities_end_force(p);

  /* Reset force variables if we are using Lloyd's algorithm. */
  hydro_gizmo_lloyd_end_force(p);
}

/**
 * @brief Extra operations done during the kick
 *
 * @param p Particle to act upon.
 * @param xp Extended particle data to act upon.
 * @param dt_therm Thermal energy time-step @f$\frac{dt}{a^2}@f$.
 * @param dt_grav Gravity time-step @f$\frac{dt}{a}@f$.
 * @param dt_hydro Hydro acceleration time-step
 * @f$\frac{dt}{a^{3(\gamma{}-1)}}@f$.
 * @param dt_kick_corr Gravity correction time-step @f$adt@f$.
 * @param cosmo Cosmology.
 * @param hydro_props Additional hydro properties.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part* p, struct xpart* xp, float dt_therm, float dt_grav,
    float dt_grav_mesh, float dt_hydro, float dt_kick_corr,
    const struct cosmology* cosmo, const struct hydro_props* hydro_props,
    const struct entropy_floor_properties* floor_props) {

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

    p->conserved.energy += hydro_gizmo_mfv_gravity_energy_update_term(
        dt_kick_corr, p, p->conserved.momentum, a_grav, grav_kick_factor);
  }

  if (p->flux.dt > 0.0f) {
    float flux[5];
    hydro_part_get_fluxes(p, flux);

    /* Update conserved variables. */
    p->conserved.mass += hydro_gizmo_mfv_mass_update_term(flux[0], dt_therm);
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

    /* reset the fluxes, so that they do not get used again in kick1 */
    hydro_part_reset_hydro_fluxes(p);
    /* invalidate the particle time step. It is considered to be inactive until
       dt is set again in hydro_prepare_force() */
    p->flux.dt = -1.0f;
  } else if (p->flux.dt == 0.0f) {
    /* something tricky happens at the beginning of the simulation: the flux
       exchange is done for all particles, but using a time step of 0. This
       in itself is not a problem. However, it causes some issues with the
       initialisation of flux.dt for inactive particles, since this value will
       remain 0 until the particle is active again, and its flux.dt is set to
       the actual time step in hydro_prepare_force(). We have to make sure it
       is properly set to -1 here, so that inactive particles are indeed found
       to be inactive during the flux loop. */
    p->flux.dt = -1.0f;
  }

#ifndef HYDRO_GAMMA_5_3

  const float Pcorr = (dt_hydro - dt_therm) * p->geometry.volume;
  p->conserved.momentum[0] -= Pcorr * p->gradients.P[0];
  p->conserved.momentum[1] -= Pcorr * p->gradients.P[1];
  p->conserved.momentum[2] -= Pcorr * p->gradients.P[2];
#ifdef GIZMO_TOTAL_ENERGY
  p->conserved.energy -= Pcorr * (p->fluid_v[0] * p->gradients.P[0] +
                                  p->fluid_v[1] * p->gradients.P[1] +
                                  p->fluid_v[2] * p->gradients.P[2]);
#endif
#endif

  /* Apply the minimal energy limit */
  const float min_energy =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;
  if (p->conserved.energy < min_energy * p->conserved.mass) {
    p->conserved.energy = min_energy * p->conserved.mass;
    p->flux.energy = 0.0f;
  }

  // MATTHIEU: Apply the entropy floor here.

  gizmo_check_physical_quantities(
      "mass", "energy", p->conserved.mass, p->conserved.momentum[0],
      p->conserved.momentum[1], p->conserved.momentum[2], p->conserved.energy);

#ifdef SWIFT_DEBUG_CHECKS
  /* Note that this check will only have effect if no GIZMO_UNPHYSICAL option
     was selected. */
#ifdef GIZMO_MFV_SPH
  if (p->conserved.mass < 0.) {
    error(
        "Negative mass after conserved variables update (mass: %g, dmass: %g)!",
        p->conserved.mass, p->flux.mass);
  }
#endif

  if (p->conserved.energy < 0.) {
    error(
        "Negative energy after conserved variables update (energy: %g, "
        "denergy: %g)!",
        p->conserved.energy, p->flux.energy);
  }
#endif

  hydro_gizmo_mfv_update_gpart_mass(p);
  hydro_velocities_set(p, xp);

  /* undo the flux exchange and kick the particles towards their centroid */
  hydro_gizmo_lloyd_kick(p, xp, dt_therm);

  /* reset wcorr */
  p->geometry.wcorr = 1.0f;
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
    const struct part* p, const struct xpart* xp, const double time) {}

#endif /* SWIFT_GIZMO_HYDRO_H */
