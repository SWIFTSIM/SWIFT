/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PRESSURE_ENERGY_HYDRO_H
#define SWIFT_PRESSURE_ENERGY_HYDRO_H

/**
 * @file PressureEnergy/hydro.h
 * @brief Pressure-Energy implementation of SPH (Non-neighbour loop
 * equations)
 *
 * The thermal variable is the energy (u) and the pressure is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows equations (16), (17) and (18) of Hopkins, P., MNRAS, 2013,
 * Volume 428, Issue 4, pp. 2840-2856 with a simple Balsara viscosity term.
 */

#include "adiabatic_index.h"
#include "approx_math.h"
#include "dimension.h"
#include "equation_of_state.h"
#include "hydro_properties.h"
#include "kernel_hydro.h"
#include "minmax.h"

#include <float.h>

/**
 * @brief Returns the internal energy of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy(
    const struct part *restrict p) {

  return p->u;
}

/**
 * @brief Returns the pressure of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_pressure(
    const struct part *restrict p) {

  return gas_pressure_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the entropy of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_entropy(
    const struct part *restrict p) {

  return gas_entropy_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the sound speed of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_soundspeed(
    const struct part *restrict p) {

  return p->force.soundspeed;
}

/**
 * @brief Returns the physical density of a particle
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
 * @brief Returns the velocities drifted to the current time of a particle.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle.
 * @param dt The time since the last kick.
 * @param v (return) The velocities at the current time.
 */
__attribute__((always_inline)) INLINE static void hydro_get_drifted_velocities(
    const struct part *restrict p, const struct xpart *xp, float dt,
    float v[3]) {

  v[0] = xp->v_full[0] + p->a_hydro[0] * dt;
  v[1] = xp->v_full[1] + p->a_hydro[1] * dt;
  v[2] = xp->v_full[2] + p->a_hydro[2] * dt;
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_energy_dt(
    const struct part *restrict p) {
  /* Since we track energy, we can forget about entropy conversions etc. */
  return p->u_dt;
}

/**
 * @brief Sets the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest.
 * @param du_dt The new time derivative of the internal energy.
 */
__attribute__((always_inline)) INLINE static void hydro_set_energy_dt(
    struct part *restrict p, float du_dt) {
  /* Since we track energy, we can forget about entropy conversions etc. */
  p->u_dt = du_dt;
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

  /* U change limited */
  const float dt_u_changes = 
    (p->u_dt != 0.0f) ? fabsf(const_max_u_change * p->u / p->u_dt)
                      : FLT_MAX;

  return min(dt_u_changes, dt_cfl);
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
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part *restrict p, const struct hydro_space *hs) {

  p->rho = 0.f;
  p->pressure_bar = 0.f;
  p->pressure_bar_dh = 0.f;
  p->density.wcount = 0.f;
  p->density.wcount_dh = 0.f;

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
  p->density.wcount += kernel_root;
  p->density.wcount_dh -= hydro_dimension * kernel_root * h_inv;
  p->rho += p->mass * kernel_root;
  p->pressure_bar += p->mass * p->u * kernel_root;
  p->pressure_bar_dh -=
    p->mass * p->u * h_inv * (hydro_dimension * kernel_root);

  /* Finish the calculation by inserting the missing h-factors */
  p->density.wcount *= h_inv_dim;
  p->density.wcount_dh *= h_inv_dim;
  p->pressure_bar *= h_inv_dim;
  p->pressure_bar_dh *= h_inv_dim;

  /* Insert gamma factors */
  p->pressure_bar *= hydro_gamma_minus_one;
  p->pressure_bar_dh *= hydro_gamma_minus_one;

  /* Finish calculation of the velocity curl components */
  const float factor = h_inv_dim_plus_one / p->rho;
  p->density.rot_v[0] *= factor;
  p->density.rot_v[1] *= factor;
  p->density.rot_v[2] *= factor;

  /* Finish calculation of the velocity divergence */
  p->density.div_v *= factor;
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
  p->density.wcount = kernel_root * kernel_norm * h_inv_dim;
  p->density.wcount_dh = 0.f;
  p->pressure_bar =
    kernel_root * p->mass * p->u * hydro_gamma_minus_one * h_inv_dim;
  p->rho = kernel_root * p->mass * h_inv_dim;
  p->pressure_bar_dh = 0;
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
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part *restrict p, struct xpart *restrict xp) {
  
  const float fac_mu = 1.f; /* Will change with cosmological integration */

  /* Compute the norm of the curl */
  const float curl_v = sqrtf(p->density.rot_v[0] * p->density.rot_v[0] +
                             p->density.rot_v[1] * p->density.rot_v[1] +
                             p->density.rot_v[2] * p->density.rot_v[2]);

  /* Compute the norm of div v */
  const float abs_div_v = fabsf(p->density.div_v);

  /* Compute the pressure */
  const float pressure =
    gas_pressure_from_internal_energy(p->rho, p->u);

  /* Compute the sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Compute this `useful' property */
  const float rho_inv = 1.f / p->rho;
  const float P_over_rho2 = pressure * rho_inv * rho_inv;

  /* Compute the Balsara switch */
  const float balsara =
    abs_div_v / (abs_div_v + curl_v + 0.0001f * soundspeed / fac_mu / p->h);


  /* Update variables */
  p->force.soundspeed = soundspeed;
  p->force.balsara = balsara;
  p->force.P_over_rho2 = P_over_rho2;
  p->force.f = 1.f; /* This cannot be computed individually in Pressure-Energy */
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
  p->force.h_dt = 0.0f;
  p->u_dt = 0.0f;

  /* Reset maximal signal velocity */
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

  p->u = xp->u_full;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param p The particle
 * @param xp The extended data of the particle
 * @param dt The drift time-step.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part *restrict p, const struct xpart *restrict xp, float dt) {

  const float h_inv = 1.f / p->h;

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt;
  if (fabsf(w1) < 0.2f)
    p->h *= approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    p->h *= expf(w1);

  /* Predict density */

  const float w2 = -hydro_dimension * w1;
  if (fabsf(w2) < 0.2f) {
    p->rho *= approx_expf(w2); /* 4th order expansion of exp(w) */
    p->pressure_bar *= approx_expf(w2);
  } else {
    p->rho *= expf(w2);
    p->pressure_bar *= expf(w2);
  }

  /* Predict the energy */
  p->u += p->u_dt * dt;

  /* Compute the pressure */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  p->force.soundspeed = soundspeed;
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
    struct part *restrict p, struct xpart *restrict xp, float dt) {

  /* Do not decrease the energy (temperature) by more than a factor of 2*/
  if (dt > 0. && p->u_dt * dt < -0.5f * p->u) {
    p->u_dt = -0.5f * p->u / dt;
  }
  xp->u_full += p->u_dt * dt;

  /* Compute the pressure */
  const float pressure =
      gas_pressure_from_internal_energy(p->rho, xp->u_full);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Update the sound speed */
  p->force.soundspeed = soundspeed;
}

/**
 *  @brief Converts hydro quantity of a particle at the start of a run
 *
 * Requires the density to be known
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *restrict p, struct xpart *restrict xp) {

  p->entropy = hydro_get_entropy(p);

  /* Compute the pressure */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  /* Compute the sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  const float rho_inv = 1.f / p->rho;
  const float P_over_rho2 = pressure * rho_inv * rho_inv;

  /* Update sound speed */
  p->force.soundspeed = soundspeed;
  p->force.P_over_rho2 = P_over_rho2;
}

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
  p->rho = 0.f;
  p->pressure_bar = 0.f;
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];

  xp->u_full = p->u;

  hydro_reset_acceleration(p);
  hydro_init_part(p, NULL);
}

/* Begin actual helper functions for the specific hydro implementation.
 *
 * Here we define:
 *
 * hydro_h_term(*pi, *pj)
 */

/**
 *  @brief Calculates the f_ij factors for a given particle i and j (i.e. the
 *         'h-terms'.
 *
 *  Requires the relevant quantities to be known; here (in Pressure Energy) we
 *  require: 
 *    - pressure_bar_dh
 *    - number density (wcount)
 *    - internal_energy
 *    - number density dh (wcount_dh)
 *
 *  @param pi The particle i
 *  @param pj The particle j
 *  @param hi The smoothing length of particle i
 */
__attribute__((always_inline)) INLINE static float hydro_h_term(
    struct part *restrict pi, struct part *restrict pj, float hi) {
  
  /* First constrct h_i / (n_D n_i) */
  
  const float hi_over_density = hi / (hydro_dimension * pi->density.wcount);

  /* Construct the 'bottom term' */
  const float bottom_term = 1 + hi_over_density * pi->density.wcount_dh;

  /* Construct the 'top term' */
  const float top_term_top = pi->pressure_bar_dh * hi_over_density;
  const float top_term_bottom =
    hydro_gamma_minus_one * pj->mass * pj->u;

  /* Putting it all together... */
  const float f_ij = 1 - (top_term_top / (bottom_term * top_term_bottom));

  return f_ij;
}

#endif /* SWIFT_PRESSURE_ENERGY_HYDRO_H */
