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
#ifndef SWIFT_DEFAULT_HYDRO_H
#define SWIFT_DEFAULT_HYDRO_H

#include "adiabatic_index.h"
#include "approx_math.h"
#include "equation_of_state.h"
#include "hydro_space.h"
#include "minmax.h"

#include <float.h>

/**
 * @brief Returns the internal energy of a particle
 *
 * @param p The particle of interest
 * @param dt Time since the last kick
 */
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy(
    const struct part *restrict p) {

  return p->u;
}

/**
 * @brief Returns the pressure of a particle
 *
 * @param p The particle of interest
 * @param dt Time since the last kick
 */
__attribute__((always_inline)) INLINE static float hydro_get_pressure(
    const struct part *restrict p) {

  return gas_pressure_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the entropy of a particle
 *
 * @param p The particle of interest
 * @param dt Time since the last kick
 */
__attribute__((always_inline)) INLINE static float hydro_get_entropy(
    const struct part *restrict p) {

  return gas_entropy_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the sound speed of a particle
 *
 * @param p The particle of interest
 * @param dt Time since the last kick
 */
__attribute__((always_inline)) INLINE static float hydro_get_soundspeed(
    const struct part *restrict p) {

  return p->force.soundspeed;
}

/**
 * @brief Returns the density of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_density(
    const struct part *restrict p) {

  return p->rho;
}

/**
 * @brief Returns the mass of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_mass(
    const struct part *restrict p) {

  return p->mass;
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy_dt(
    const struct part *restrict p) {

  return p->force.u_dt;
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest.
 * @param du_dt The new time derivative of the internal energy.
 */
__attribute__((always_inline)) INLINE static void hydro_set_internal_energy_dt(
    struct part *restrict p, float du_dt) {

  p->force.u_dt = du_dt;
}

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 *
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct hydro_props *restrict hydro_properties) {

  const float CFL_condition = hydro_properties->CFL_condition;

  /* CFL condition */
  const float dt_cfl =
      2.f * kernel_gamma * CFL_condition * p->h / p->force.v_sig;

  /* Limit change in u */
  const float dt_u_change =
      (p->force.u_dt != 0.0f) ? fabsf(const_max_u_change * p->u / p->force.u_dt)
                              : FLT_MAX;

  return min(dt_cfl, dt_u_change);
}

/**
 * @brief Does some extra hydro operations once the actual physical time step
 * for the particle is known.
 *
 * @param p The particle to act upon.
 * @param dt Physical time step of the particle during the next step.
 */
__attribute__((always_inline)) INLINE static void hydro_timestep_extra(
    struct part *p, float dt) {}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the variaous density tasks
 *
 * @param p The particle to act upon
 * @param hs #hydro_space containing hydro specific space information.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part *restrict p, const struct hydro_space *hs) {
  p->density.wcount = 0.f;
  p->density.wcount_dh = 0.f;
  p->rho = 0.f;
  p->rho_dh = 0.f;
  p->density.div_v = 0.f;
  p->density.rot_v[0] = 0.f;
  p->density.rot_v[1] = 0.f;
  p->density.rot_v[2] = 0.f;
}

/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part *restrict p) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Final operation on the density (add self-contribution). */
  p->rho += p->mass * kernel_root;
  p->rho_dh -= hydro_dimension * p->mass * kernel_root;
  p->density.wcount += kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  p->rho *= h_inv_dim;
  p->rho_dh *= h_inv_dim_plus_one;
  p->density.wcount *= kernel_norm;
  p->density.wcount_dh *= h_inv * kernel_gamma * kernel_norm;

  const float irho = 1.f / p->rho;

  /* Finish calculation of the velocity curl components */
  p->density.rot_v[0] *= h_inv_dim_plus_one * irho;
  p->density.rot_v[1] *= h_inv_dim_plus_one * irho;
  p->density.rot_v[2] *= h_inv_dim_plus_one * irho;

  /* Finish calculation of the velocity divergence */
  p->density.div_v *= h_inv_dim_plus_one * irho;
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_part_has_no_neighbours(
    struct part *restrict p, struct xpart *restrict xp) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  p->rho = p->mass * kernel_root * h_inv_dim;
  p->density.wcount = kernel_root * kernel_norm * h_inv_dim;
  p->rho_dh = 0.f;
  p->density.wcount_dh = 0.f;
  p->density.div_v = 0.f;
  p->density.rot_v[0] = 0.f;
  p->density.rot_v[1] = 0.f;
  p->density.rot_v[2] = 0.f;
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * Computes viscosity term, conduction term and smoothing length gradient terms.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param time The current time
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part *restrict p, struct xpart *restrict xp) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;

  /* Pre-compute some stuff for the balsara switch. */
  const float normDiv_v = fabs(p->density.div_v);
  const float normRot_v = sqrtf(p->density.rot_v[0] * p->density.rot_v[0] +
                                p->density.rot_v[1] * p->density.rot_v[1] +
                                p->density.rot_v[2] * p->density.rot_v[2]);

  /* Compute this particle's sound speed. */
  const float u = p->u;
  const float fc = p->force.soundspeed =
      sqrtf(hydro_gamma * hydro_gamma_minus_one * u);

  /* Compute the derivative term */
  p->rho_dh = 1.f / (1.f + hydro_dimension_inv * p->h * p->rho_dh / p->rho);

  /* Compute the P/Omega/rho2. */
  xp->omega = 1.0f + hydro_dimension_inv * h * p->rho_dh / p->rho;
  p->force.P_over_rho2 = u * hydro_gamma_minus_one / (p->rho * xp->omega);

  /* Balsara switch */
  p->force.balsara = normDiv_v / (normDiv_v + normRot_v + 0.0001f * fc * h_inv);

  /* Viscosity parameter decay time */
  /* const float tau = h / (2.f * const_viscosity_length * p->force.soundspeed);
   */

  /* Viscosity source term */
  /* const float S = max(-normDiv_v, 0.f); */

  /* Compute the particle's viscosity parameter time derivative */
  /* const float alpha_dot = (const_viscosity_alpha_min - p->alpha) / tau + */
  /*                         (const_viscosity_alpha_max - p->alpha) * S; */

  /* Update particle's viscosity paramter */
  /* p->alpha += alpha_dot * (p->ti_end - p->ti_begin) * timeBase; */  // MATTHIEU
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking place in the variaous force tasks
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_reset_acceleration(
    struct part *restrict p) {

  /* Reset the acceleration. */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* Reset the time derivatives. */
  p->force.u_dt = 0.0f;
  p->force.h_dt = 0.0f;
  p->force.v_sig = 0.0f;
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values(
    struct part *restrict p, const struct xpart *restrict xp) {

  /* Re-set the predicted velocities */
  p->v[0] = xp->v_full[0];
  p->v[1] = xp->v_full[1];
  p->v[2] = xp->v_full[2];
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param p The particle
 * @param xp The extended data of the particle
 * @param dt The drift time-step.
 * @param t0 The time at the start of the drift
 * @param t1 The time at the end of the drift
 * @param timeBase The minimal time-step size
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part *restrict p, struct xpart *restrict xp, float dt) {
  float u, w;

  const float h_inv = 1.f / p->h;

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt;
  if (fabsf(w1) < 0.2f)
    p->h *= approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    p->h *= expf(w1);

  /* Predict density */
  const float w2 = -hydro_dimension * w1;
  if (fabsf(w2) < 0.2f)
    p->rho *= approx_expf(w2); /* 4th order expansion of exp(w) */
  else
    p->rho *= expf(w2);

  /* Predict internal energy */
  w = p->force.u_dt / p->u * dt;
  if (fabsf(w) < 0.2f)
    u = p->u *= approx_expf(w);
  else
    u = p->u *= expf(w);

  /* Predict gradient term */
  p->force.P_over_rho2 = u * hydro_gamma_minus_one / (p->rho * xp->omega);
}

/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the forces and accelerationsby the appropiate constants
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_end_force(
    struct part *restrict p) {
  p->force.h_dt *= p->h * hydro_dimension_inv;
}

/**
 * @brief Kick the additional variables
 *
 * @param p The particle to act upon
 * @param xp The particle extended data to act upon
 * @param dt The time-step for this kick
 * @param half_dt The half time-step for this kick
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part *restrict p, struct xpart *restrict xp, float dt) {}

/**
 *  @brief Converts hydro quantity of a particle at the start of a run
 *
 * Requires the density to be known
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *restrict p, struct xpart *restrict xp) {}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_first_init_part(
    struct part *restrict p, struct xpart *restrict xp) {

  p->time_bin = 0;
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
  xp->u_full = p->u;

  hydro_reset_acceleration(p);
  hydro_init_part(p, NULL);
}

#endif /* SWIFT_DEFAULT_HYDRO_H */
