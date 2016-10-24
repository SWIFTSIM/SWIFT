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
#ifndef SWIFT_PRESSURE_ENTROPY_HYDRO_H
#define SWIFT_PRESSURE_ENTROPY_HYDRO_H

/**
 * @file PressureEntropy/hydro.h
 * @brief Pressure-Entropy implementation of SPH (Non-neighbour loop
 * equations)
 *
 * The thermal variable is the entropy (S) and the entropy is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows eqautions (19), (21) and (22) of Hopkins, P., MNRAS, 2013,
 * Volume 428, Issue 4, pp. 2840-2856 with a simple Balsara viscosity term.
 */

#include "adiabatic_index.h"
#include "approx_math.h"
#include "dimension.h"
#include "equation_of_state.h"
#include "hydro_properties.h"
#include "kernel_hydro.h"
#include "minmax.h"

/**
 * @brief Returns the internal energy of a particle
 *
 * @param p The particle of interest
 * @param dt Time since the last kick
 */
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy(
    const struct part *restrict p, float dt) {

  const float entropy = p->entropy + p->entropy_dt * dt;

  return gas_internal_energy_from_entropy(p->rho_bar, entropy);
}

/**
 * @brief Returns the pressure of a particle
 *
 * @param p The particle of interest
 * @param dt Time since the last kick
 */
__attribute__((always_inline)) INLINE static float hydro_get_pressure(
    const struct part *restrict p, float dt) {

  const float entropy = p->entropy + p->entropy_dt * dt;

  return gas_pressure_from_entropy(p->rho_bar, entropy);
}

/**
 * @brief Returns the entropy of a particle
 *
 * @param p The particle of interest
 * @param dt Time since the last kick
 */
__attribute__((always_inline)) INLINE static float hydro_get_entropy(
    const struct part *restrict p, float dt) {

  return p->entropy + p->entropy_dt * dt;
}

/**
 * @brief Returns the sound speed of a particle
 *
 * @param p The particle of interest
 * @param dt Time since the last kick
 */
__attribute__((always_inline)) INLINE static float hydro_get_soundspeed(
    const struct part *restrict p, float dt) {

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
 * @brief Modifies the thermal state of a particle to the imposed internal
 * energy
 *
 * This overwrites the current state of the particle but does *not* change its
 * time-derivatives. Entropy, pressure, sound-speed and signal velocity will be
 * updated.
 *
 * @param p The particle
 * @param u The new internal energy
 */
__attribute__((always_inline)) INLINE static void hydro_set_internal_energy(
    struct part *restrict p, float u) {

  p->entropy = gas_entropy_from_internal_energy(p->rho_bar, u);
  p->entropy_one_over_gamma = pow_one_over_gamma(p->entropy);

  /* Compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho_bar, p->entropy);

  /* Compute the sound speed from the pressure*/
  const float soundspeed = gas_soundspeed_from_pressure(p->rho_bar, pressure);

  /* Update the signal velocity */
  const float v_sig_old = p->force.v_sig;
  const float v_sig_new = p->force.v_sig - p->force.soundspeed + soundspeed;
  const float v_sig = max(v_sig_old, v_sig_new);

  const float rho_bar_inv = 1.f / p->rho_bar;

  p->force.soundspeed = soundspeed;
  p->force.P_over_rho2 = pressure * rho_bar_inv * rho_bar_inv;
  p->force.v_sig = v_sig;
}

/**
 * @brief Modifies the thermal state of a particle to the imposed entropy
 *
 * This overwrites the current state of the particle but does *not* change its
 * time-derivatives. Entropy, pressure, sound-speed and signal velocity will be
 * updated.
 *
 * @param p The particle
 * @param S The new entropy
 */
__attribute__((always_inline)) INLINE static void hydro_set_entropy(
    struct part *restrict p, float S) {

  p->entropy = S;
  p->entropy_one_over_gamma = pow_one_over_gamma(p->entropy);

  /* Compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho_bar, p->entropy);

  /* Compute the sound speed from the pressure*/
  const float soundspeed = gas_soundspeed_from_pressure(p->rho_bar, pressure);

  /* Update the signal velocity */
  const float v_sig_old = p->force.v_sig;
  const float v_sig_new = p->force.v_sig - p->force.soundspeed + soundspeed;
  const float v_sig = max(v_sig_old, v_sig_new);

  const float rho_bar_inv = 1.f / p->rho_bar;

  p->force.soundspeed = soundspeed;
  p->force.P_over_rho2 = pressure * rho_bar_inv * rho_bar_inv;
  p->force.v_sig = v_sig;
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

  return dt_cfl;
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

  p->ti_begin = 0;
  p->ti_end = 0;
  p->rho_bar = 0.f;
  p->entropy_one_over_gamma = pow_one_over_gamma(p->entropy);
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the variaous density tasks
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part *restrict p) {

  p->rho = 0.f;
  p->rho_bar = 0.f;
  p->density.wcount = 0.f;
  p->density.wcount_dh = 0.f;
  p->density.rho_dh = 0.f;
  p->density.pressure_dh = 0.f;

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
 * @param ti_current The current time (on the integer timeline)
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part *restrict p, int ti_current) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Final operation on the density (add self-contribution). */
  p->rho += p->mass * kernel_root;
  p->rho_bar += p->mass * p->entropy_one_over_gamma * kernel_root;
  p->density.rho_dh -= hydro_dimension * p->mass * kernel_root;
  p->density.pressure_dh -=
      hydro_dimension * p->mass * p->entropy_one_over_gamma * kernel_root;
  p->density.wcount += kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  p->rho *= h_inv_dim;
  p->rho_bar *= h_inv_dim;
  p->density.rho_dh *= h_inv_dim_plus_one;
  p->density.pressure_dh *= h_inv_dim_plus_one;
  p->density.wcount *= kernel_norm;
  p->density.wcount_dh *= h_inv * kernel_gamma * kernel_norm;

  const float rho_inv = 1.f / p->rho;
  const float entropy_minus_one_over_gamma = 1.f / p->entropy_one_over_gamma;

  /* Final operation on the weighted density */
  p->rho_bar *= entropy_minus_one_over_gamma;

  /* Finish calculation of the velocity curl components */
  p->density.rot_v[0] *= h_inv_dim_plus_one * rho_inv;
  p->density.rot_v[1] *= h_inv_dim_plus_one * rho_inv;
  p->density.rot_v[2] *= h_inv_dim_plus_one * rho_inv;

  /* Finish calculation of the velocity divergence */
  p->density.div_v *= h_inv_dim_plus_one * rho_inv;
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * Computes viscosity term, conduction term and smoothing length gradient terms.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param ti_current The current time (on the timeline)
 * @param timeBase The minimal time-step size
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part *restrict p, struct xpart *restrict xp, int ti_current,
    double timeBase) {

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  /* Compute the norm of the curl */
  const float curl_v = sqrtf(p->density.rot_v[0] * p->density.rot_v[0] +
                             p->density.rot_v[1] * p->density.rot_v[1] +
                             p->density.rot_v[2] * p->density.rot_v[2]);

  /* Compute the norm of div v */
  const float abs_div_v = fabsf(p->density.div_v);

  /* Compute the pressure */
  const float half_dt = (ti_current - (p->ti_begin + p->ti_end) / 2) * timeBase;
  const float entropy = hydro_get_entropy(p, half_dt);
  const float pressure = gas_pressure_from_entropy(p->rho_bar, entropy);

  /* Compute the sound speed from the pressure*/
  const float soundspeed = gas_soundspeed_from_pressure(p->rho_bar, pressure);

  /* Compute the Balsara switch */
  const float balsara =
      abs_div_v / (abs_div_v + curl_v + 0.0001f * soundspeed / fac_mu / p->h);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_bar_inv = 1.f / p->rho_bar;
  const float P_over_rho2 = pressure * rho_bar_inv * rho_bar_inv;

  /* Compute "grad h" term (note we use rho here and not rho_bar !)*/
  const float rho_inv = 1.f / p->rho;
  const float rho_dh =
      1.f / (1.f + hydro_dimension_inv * p->h * p->density.rho_dh * rho_inv);
  const float pressure_dh =
      p->density.pressure_dh * rho_inv * p->h * hydro_dimension_inv;

  const float grad_h_term = rho_dh * pressure_dh;

  /* Update variables. */
  p->force.soundspeed = soundspeed;
  p->force.P_over_rho2 = P_over_rho2;
  p->force.balsara = balsara;
  p->force.f = grad_h_term;
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
  p->entropy_dt = 0.0f;
  p->force.h_dt = 0.0f;

  /* Reset maximal signal velocity */
  p->force.v_sig = 0.0f;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param p The particle
 * @param xp The extended data of the particle
 * @param dt The drift time-step.
 * @param t0 The time at the start of the drift (on the timeline).
 * @param t1 The time at the end of the drift (on the timeline).
 * @param timeBase The minimal time-step size
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part *restrict p, const struct xpart *restrict xp, float dt, int t0,
    int t1, double timeBase) {

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
    p->rho_bar *= approx_expf(w2);
  } else {
    p->rho *= expf(w2);
    p->rho_bar *= expf(w2);
  }

  /* Drift the entropy */
  const float dt_entr = (t1 - (p->ti_begin + p->ti_end) / 2) * timeBase;
  const float entropy = hydro_get_entropy(p, dt_entr);

  /* Compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho_bar, entropy);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho_bar, pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_bar_inv = 1.f / p->rho_bar;
  const float P_over_rho2 = pressure * rho_bar_inv * rho_bar_inv;

  /* Update the variables */
  p->entropy_one_over_gamma = pow_one_over_gamma(entropy);
  p->force.soundspeed = soundspeed;
  p->force.P_over_rho2 = P_over_rho2;
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

  p->entropy_dt =
      0.5f * gas_entropy_from_internal_energy(p->rho_bar, p->entropy_dt);
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
    struct part *restrict p, struct xpart *restrict xp, float dt,
    float half_dt) {

  /* Do not decrease the entropy (temperature) by more than a factor of 2*/
  const float entropy_change = p->entropy_dt * dt;
  if (entropy_change > -0.5f * p->entropy)
    p->entropy += entropy_change;
  else
    p->entropy *= 0.5f;

  /* Do not 'overcool' when timestep increases */
  if (p->entropy + p->entropy_dt * half_dt < 0.5f * p->entropy)
    p->entropy_dt = -0.5f * p->entropy / half_dt;

  /* Compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho_bar, p->entropy);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho_bar, pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_bar_inv = 1.f / p->rho_bar;
  const float P_over_rho2 = pressure * rho_bar_inv * rho_bar_inv;

  p->entropy_one_over_gamma = pow_one_over_gamma(p->entropy);
  p->force.soundspeed = soundspeed;
  p->force.P_over_rho2 = P_over_rho2;
}

/**
 *  @brief Converts hydro quantity of a particle at the start of a run
 *
 * Requires the density to be known
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *restrict p) {

  /* We read u in the entropy field. We now get S from u */
  p->entropy = gas_entropy_from_internal_energy(p->rho_bar, p->entropy);
  p->entropy_one_over_gamma = pow_one_over_gamma(p->entropy);

  /* Compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho_bar, p->entropy);

  /* Compute the sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho_bar, pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_bar_inv = 1.f / p->rho_bar;
  const float P_over_rho2 = pressure * rho_bar_inv * rho_bar_inv;

  p->force.soundspeed = soundspeed;
  p->force.P_over_rho2 = P_over_rho2;
}

#endif /* SWIFT_PRESSURE_ENTROPY_HYDRO_H */
