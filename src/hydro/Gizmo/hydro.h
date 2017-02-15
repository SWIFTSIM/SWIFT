/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "minmax.h"

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

  return CFL_condition * p->h / fabsf(p->timestepvars.vmax);
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

  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];

  p->primitives.v[0] = p->v[0];
  p->primitives.v[1] = p->v[1];
  p->primitives.v[2] = p->v[2];
}

/**
 * @brief Prepares a particle for the volume calculation.
 *
 * Simply makes sure all necessary variables are initialized to zero.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part* p) {

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

  /* Final operation on the density. */
  p->density.wcount += kernel_root;
  p->density.wcount *= kernel_norm;

  p->density.wcount_dh *= ih * kernel_gamma * kernel_norm;

  const float ihdim = pow_dimension(ih);

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

  invert_dimension_by_dimension_matrix(p->geometry.matrix_E);

  hydro_gradients_init(p);

  /* compute primitive variables */
  /* eqns (3)-(5) */
  const float m = p->conserved.mass;
  float momentum[3];
  momentum[0] = p->conserved.momentum[0];
  momentum[1] = p->conserved.momentum[1];
  momentum[2] = p->conserved.momentum[2];
  p->primitives.rho = m / volume;
  p->primitives.v[0] = momentum[0] / m;
  p->primitives.v[1] = momentum[1] / m;
  p->primitives.v[2] = momentum[2] / m;
  const float energy = p->conserved.energy;
  p->primitives.P = hydro_gamma_minus_one * energy / volume;
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
 * @param ti_current Current integer time.
 * @param timeBase Conversion factor between integer time and physical time.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part* restrict p, struct xpart* restrict xp) {

  /* Set the physical time step */
  p->force.dt = get_timestep(p->time_bin, 0.);  // MATTHIEU 0

  /* Initialize time step criterion variables */
  p->timestepvars.vmax = 0.0f;

  /* Set the actual velocity of the particle */
  p->force.v_full[0] = xp->v_full[0];
  p->force.v_full[1] = xp->v_full[1];
  p->force.v_full[2] = xp->v_full[2];
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
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part* p, struct xpart* xp) {

  const float volume = p->geometry.volume;
  const float m = p->conserved.mass;
  p->primitives.rho = m / volume;

  /* first get the initial velocities, as they were overwritten in end_density
   */
  p->primitives.v[0] = p->v[0];
  p->primitives.v[1] = p->v[1];
  p->primitives.v[2] = p->v[2];

  p->conserved.momentum[0] = m * p->primitives.v[0];
  p->conserved.momentum[1] = m * p->primitives.v[1];
  p->conserved.momentum[2] = m * p->primitives.v[2];

  p->primitives.P =
      hydro_gamma_minus_one * p->conserved.energy * p->primitives.rho;

  p->conserved.energy *= m;
}

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

  const float h_inv = 1.0f / p->h;

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt;
  if (fabsf(w1) < 0.2f)
    p->h *= approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    p->h *= expf(w1);

  const float w2 = -hydro_dimension * w1;
  if (fabsf(w2) < 0.2f) {
    p->primitives.rho *= approx_expf(w2);
  } else {
    p->primitives.rho *= expf(w2);
  }

  p->primitives.v[0] += (p->a_hydro[0] + p->gravity.old_a[0]) * dt;
  p->primitives.v[1] += (p->a_hydro[1] + p->gravity.old_a[1]) * dt;
  p->primitives.v[2] += (p->a_hydro[2] + p->gravity.old_a[2]) * dt;
  const float u = p->conserved.energy + p->du_dt * dt;
  p->primitives.P =
      hydro_gamma_minus_one * u * p->primitives.rho / p->conserved.mass;
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

  /* Add normalization to h_dt. */
  p->force.h_dt *= p->h * hydro_dimension_inv;

  /* Set the hydro acceleration, based on the new momentum and mass */
  /* NOTE: the momentum and mass are only correct for active particles, since
           only active particles have received flux contributions from all their
           neighbours. Since this method is only called for active particles,
           this is indeed the case. */
  if (p->force.dt) {
    float mnew;
    float vnew[3];

    mnew = p->conserved.mass + p->conserved.flux.mass;
    vnew[0] = (p->conserved.momentum[0] + p->conserved.flux.momentum[0]) / mnew;
    vnew[1] = (p->conserved.momentum[1] + p->conserved.flux.momentum[1]) / mnew;
    vnew[2] = (p->conserved.momentum[2] + p->conserved.flux.momentum[2]) / mnew;

    p->a_hydro[0] = (vnew[0] - p->force.v_full[0]) / p->force.dt;
    p->a_hydro[1] = (vnew[1] - p->force.v_full[1]) / p->force.dt;
    p->a_hydro[2] = (vnew[2] - p->force.v_full[2]) / p->force.dt;

    p->du_dt = p->conserved.flux.energy / p->force.dt;
  } else {
    p->a_hydro[0] = 0.0f;
    p->a_hydro[1] = 0.0f;
    p->a_hydro[2] = 0.0f;

    p->du_dt = 0.0f;
  }
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

  float oldm, oldp[3], anew[3];
  const float half_dt = 0.5f * dt;  // MATTHIEU

  /* Retrieve the current value of the gravitational acceleration from the
     gpart. We are only allowed to do this because this is the kick. We still
     need to check whether gpart exists though.*/
  if (p->gpart) {
    anew[0] = p->gpart->a_grav[0];
    anew[1] = p->gpart->a_grav[1];
    anew[2] = p->gpart->a_grav[2];

    /* Copy the old mass and momentum before updating the conserved variables */
    oldm = p->conserved.mass;
    oldp[0] = p->conserved.momentum[0];
    oldp[1] = p->conserved.momentum[1];
    oldp[2] = p->conserved.momentum[2];
  }

  /* Update conserved variables. */
  p->conserved.mass += p->conserved.flux.mass;
  p->conserved.momentum[0] += p->conserved.flux.momentum[0];
  p->conserved.momentum[1] += p->conserved.flux.momentum[1];
  p->conserved.momentum[2] += p->conserved.flux.momentum[2];
  p->conserved.energy += p->conserved.flux.energy;

  /* Add gravity. We only do this if we have gravity activated. */
  if (p->gpart) {
    p->conserved.momentum[0] +=
        half_dt * (oldm * p->gravity.old_a[0] + p->conserved.mass * anew[0]);
    p->conserved.momentum[1] +=
        half_dt * (oldm * p->gravity.old_a[1] + p->conserved.mass * anew[1]);
    p->conserved.momentum[2] +=
        half_dt * (oldm * p->gravity.old_a[2] + p->conserved.mass * anew[2]);

    float paold, panew;
    paold = oldp[0] * p->gravity.old_a[0] + oldp[1] * p->gravity.old_a[1] +
            oldp[2] * p->gravity.old_a[2];
    panew = p->conserved.momentum[0] * anew[0] +
            p->conserved.momentum[1] * anew[1] +
            p->conserved.momentum[2] * anew[2];
    p->conserved.energy += half_dt * (paold + panew);

    float fluxaold, fluxanew;
    fluxaold = p->gravity.old_a[0] * p->gravity.old_mflux[0] +
               p->gravity.old_a[1] * p->gravity.old_mflux[1] +
               p->gravity.old_a[2] * p->gravity.old_mflux[2];
    fluxanew = anew[0] * p->gravity.mflux[0] + anew[1] * p->gravity.mflux[1] +
               anew[2] * p->gravity.mflux[2];
    p->conserved.energy += half_dt * (fluxaold + fluxanew);

    /* Store gravitational acceleration and mass flux for next step */
    p->gravity.old_a[0] = anew[0];
    p->gravity.old_a[1] = anew[1];
    p->gravity.old_a[2] = anew[2];
    p->gravity.old_mflux[0] = p->gravity.mflux[0];
    p->gravity.old_mflux[1] = p->gravity.mflux[1];
    p->gravity.old_mflux[2] = p->gravity.mflux[2];
  }

  /* reset fluxes */
  /* we can only do this here, since we need to keep the fluxes for inactive
     particles */
  p->conserved.flux.mass = 0.0f;
  p->conserved.flux.momentum[0] = 0.0f;
  p->conserved.flux.momentum[1] = 0.0f;
  p->conserved.flux.momentum[2] = 0.0f;
  p->conserved.flux.energy = 0.0f;
}

/**
 * @brief Returns the internal energy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy(
    const struct part* restrict p) {

  return p->primitives.P / hydro_gamma_minus_one / p->primitives.rho;
}

/**
 * @brief Returns the entropy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_entropy(
    const struct part* restrict p) {

  return p->primitives.P / pow_gamma(p->primitives.rho);
}

/**
 * @brief Returns the sound speed of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_soundspeed(
    const struct part* restrict p) {

  return sqrtf(hydro_gamma * p->primitives.P / p->primitives.rho);
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

  p->conserved.energy = gas_internal_energy_from_entropy(p->primitives.rho, S) *
                        p->conserved.mass;
  p->primitives.P = gas_pressure_from_entropy(p->primitives.rho, S);
}
