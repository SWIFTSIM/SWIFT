/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

//#define GIZMO_LLOYD_ITERATION

#include "hydro_getters.h"
#include "hydro_gradients.h"
#include "hydro_setters.h"
#include "hydro_space.h"

#if defined(GIZMO_MFV_SPH)
#include "MFV/hydro_velocities.h"
#endif

#include <float.h>

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

#ifdef GIZMO_LLOYD_ITERATION
  return CFL_condition;
#endif

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

#ifdef GIZMO_LLOYD_ITERATION
  /* overwrite all variables to make sure they have safe values */
  W[0] = 1.0f;
  W[1] = 0.0f;
  W[2] = 0.0f;
  W[3] = 0.0f;
  W[4] = 1.0f;

  Q[0] = 1.0f;
  Q[1] = 0.0f;
  Q[2] = 0.0f;
  Q[3] = 0.0f;
  Q[4] = 1.0f;

  p->v[0] = 0.0f;
  p->v[1] = 0.0f;
  p->v[2] = 0.0f;
#endif

  p->time_bin = 0;

  hydro_part_set_primitive_variables(p, W);
  hydro_part_set_conserved_variables(p, Q);

#if defined(GIZMO_MFV_SPH)
  /* initialize the particle velocity based on the primitive fluid velocity */
  hydro_velocities_init(p, xp);
#elif defined(GIZMO_MFM_SPH)
  /* initialize the particle velocity based on the primitive fluid velocity */
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
#endif

  /* ignore accelerations present in the initial condition */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* we cannot initialize wcorr in init_part, as init_part gets called every
     time the density loop is repeated, and the whole point of storing wcorr
     is to have a way of remembering that we need more neighbours for this
     particle */
  hydro_part_set_wcorr(p, 1.0f);
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

  const float wcount_inv = 1.0f / p->density.wcount;
  p->geometry.centroid[0] *= wcount_inv;
  p->geometry.centroid[1] *= wcount_inv;
  p->geometry.centroid[2] *= wcount_inv;

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

  const float condition_number =
      hydro_dimension_inv * sqrtf(condition_number_E * condition_number_Einv);

  if (condition_number > const_gizmo_max_condition_number &&
      hydro_part_get_wcorr(p) > const_gizmo_min_wcorr) {
#ifdef GIZMO_PATHOLOGICAL_ERROR
    error("Condition number larger than %g (%g)!",
          const_gizmo_max_condition_number, condition_number);
#endif
#ifdef GIZMO_PATHOLOGICAL_WARNING
    message("Condition number too large: %g (> %g, p->id: %llu)!",
            condition_number, const_gizmo_max_condition_number, p->id);
#endif
    /* add a correction to the number of neighbours for this particle */
    hydro_part_set_wcorr(
        p, const_gizmo_w_correction_factor * hydro_part_get_wcorr(p));
  }

  /* need to figure out whether or not to move this to hydro_prepare_gradient */
  hydro_gradients_init(p);

  /* compute primitive variables */
  /* eqns (3)-(5) */
  const float Q[5] = {p->conserved.mass, p->conserved.momentum[0],
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
  if (Q[0] == 0.0f) {
    W[1] = 0.;
    W[2] = 0.;
    W[3] = 0.;
  } else {
    const float m_inv = 1.0f / Q[0];
    W[1] = Q[1] * m_inv;
    W[2] = Q[2] * m_inv;
    W[3] = Q[3] * m_inv;
  }

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

#ifdef GIZMO_LLOYD_ITERATION
  /* overwrite primitive variables to make sure they still have safe values */
  W[0] = 1.0f;
  W[1] = 0.0f;
  W[2] = 0.0f;
  W[3] = 0.0f;
  W[4] = 1.0f;
#endif

  hydro_part_set_primitive_variables(p, W);

  /* Add a correction factor to wcount (to force a neighbour number increase if
     the geometry matrix is close to singular) */
  p->density.wcount *= hydro_part_get_wcorr(p);
  p->density.wcount_dh *= hydro_part_get_wcorr(p);
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
  /* centroid is relative w.r.t. particle position */
  /* by setting the centroid to 0.0f, we make sure no velocity correction is
     applied */
  p->geometry.centroid[0] = 0.0f;
  p->geometry.centroid[1] = 0.0f;
  p->geometry.centroid[2] = 0.0f;
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
    const struct cosmology* cosmo, const struct hydro_props* hydro_props) {

  /* Initialize time step criterion variables */
  p->timestepvars.vmax = 0.;

  /* not sure if we need to do this here or in end_density() */
  hydro_gradients_init(p);

#if defined(GIZMO_MFV_SPH)
  /* Set the actual velocity of the particle */
  hydro_velocities_prepare_force(p, xp);
#endif
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

#ifdef GIZMO_LLOYD_ITERATION
  /* reset the gradients to zero, as we don't want them */
  hydro_gradients_init(p);
#endif
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

  hydro_part_reset_fluxes(p);
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
    const struct cosmology* cosmo) {
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
    const struct hydro_props* hydro_props) {

  p->conserved.energy /= cosmo->a_factor_internal_energy;
}

#if defined(GIZMO_MFV_SPH)
#include "MFV/hydro.h"
#define SPH_IMPLEMENTATION "GIZMO MFV (Hopkins 2015)"
#elif defined(GIZMO_MFM_SPH)
#include "MFM/hydro.h"
#define SPH_IMPLEMENTATION "GIZMO MFM (Hopkins 2015)"
#endif

#endif /* SWIFT_GIZMO_HYDRO_H */
