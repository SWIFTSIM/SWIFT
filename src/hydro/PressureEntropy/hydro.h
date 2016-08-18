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
 * @file PressureEntropy/hydro_part.h
 * @brief Pressure-Entropy implementation of SPH (Particle definition)
 *
 * The thermal variable is the entropy (S) and the entropy is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows Hopkins, P., MNRAS, 2013, Volume 428, Issue 4, pp. 2840-2856
 */

#include "adiabatic_index.h"
#include "approx_math.h"
#include "dimension.h"
#include "equation_of_state.h"
#include "hydro_properties.h"
#include "kernel_hydro.h"

/**
 * @brief Returns the internal energy of a particle
 *
 * @param p The particle of interest
 * @param dt Time since the last kick
 */
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy(
    const struct part *restrict p, float dt) {

  const float entropy = p->entropy + p->entropy_dt * dt;

  return gas_internal_energy_from_entropy(p->rho, entropy);
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

  return gas_pressure_from_entropy(p->rho, entropy);
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
    struct part *restrict p, float u) {

  p->entropy = gas_entropy_from_internal_energy(p->rho, u);
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
    struct part *restrict p, float S) {

  p->entropy = S;
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
  p->density.wcount = 0.f;
  p->density.wcount_dh = 0.f;
  p->rho = 0.f;
  p->weightedPressure = 0.f;
  p->density.weightedPressure_dh = 0.f;
  // p->rho_dh = 0.f;
  /* p->density.div_v = 0.f; */
  /* p->density.rot_v[0] = 0.f; */
  /* p->density.rot_v[1] = 0.f; */
  /* p->density.rot_v[2] = 0.f; */
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
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  // const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Final operation on the density (add self-contribution). */
  p->rho += p->mass * kernel_root;
  p->density.wcount += kernel_root;
  // p->density.wcount_dh -= hydro_dimension * kernel_root;
  p->weightedPressure += p->mass * pow_one_over_gamma(p->entropy) * kernel_root;
  p->density.weightedPressure_dh -=
      hydro_dimension * p->mass * pow_one_over_gamma(p->entropy) * kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  p->rho *= h_inv_dim;
  p->density.wcount *= kernel_norm;
  p->density.wcount_dh *= h_inv * kernel_gamma * kernel_norm;
  p->weightedPressure *= h_inv_dim;
  p->density.weightedPressure_dh *= h_inv * kernel_gamma * kernel_norm;

  /* Final operation on the weighted pressure */
  p->weightedPressure = pow_gamma(p->weightedPressure);

  /* const float irho = 1.f / p->rho; */

  /* Compute the derivative term */
  // p->rho_dh = 1.f / (1.f + hydro_dimension_inv * p->h * p->rho_dh * irho);

  /* Finish calculation of the velocity curl components */
  // p->density.rot_v[0] *= h_inv_dim_plus_one * irho;
  // p->density.rot_v[1] *= h_inv_dim_plus_one * irho;
  // p->density.rot_v[2] *= h_inv_dim_plus_one * irho;

  /* Finish calculation of the velocity divergence */
  // p->density.div_v *= h_inv_dim_plus_one * irho;
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

  /* Compute the pressure */
  const float half_dt = (ti_current - (p->ti_begin + p->ti_end) / 2) * timeBase;

  /* Compute the sound speed from the actual pressure*/
  const float pressure = hydro_get_pressure(p, half_dt);
  const float irho = 1.f / p->rho;
  const float soundspeed = sqrtf(hydro_gamma * pressure * irho);

  /* Compute the "i" part of the f_ij term */
  const float number_density = p->density.wcount / kernel_norm;
  const float factor = p->h / (hydro_dimension * number_density);
  const float f_ij =
      factor * p->density.weightedPressure_dh / (1.f + factor * number_density);

  /* Comput the pressure term */
  const float pressure_term = pow_one_minus_two_over_gamma(p->weightedPressure);

  /* Update variables. */
  p->force.soundspeed = soundspeed;
  p->force.f_ij = f_ij;
  p->force.pressure_term = pressure_term;
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
 * @param t0 The time at the start of the drift
 * @param t1 The time at the end of the drift
 * @param timeBase The minimal time-step size
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part *restrict p, const struct xpart *restrict xp, int t0, int t1,
    double timeBase) {

  /* Drift the pressure */
  const float dt_entr = (t1 - (p->ti_begin + p->ti_end) / 2) * timeBase;
  const float pressure = hydro_get_pressure(p, dt_entr);

  const float irho = 1.f / p->rho;

  /* Compute the new sound speed */
  const float soundspeed = sqrtf(hydro_gamma * pressure * irho);

  const float w1 = -hydro_dimension * p->force.h_dt * dt_entr / p->h;
  if (fabsf(w1) < 0.2f)
    p->weightedPressure *= approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    p->weightedPressure *= expf(w1);

  const float pressure_term = pow_one_minus_two_over_gamma(p->weightedPressure);

  /* Update variables */
  p->force.soundspeed = soundspeed;
  p->force.pressure_term = pressure_term;
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

  p->entropy_dt *=
      0.5f * hydro_gamma_minus_one * pow_minus_gamma_minus_one(p->rho);
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
  p->entropy = gas_entropy_from_internal_energy(p->rho, p->entropy);
}

#endif /* SWIFT_PRESSURE_ENTROPY_HYDRO_H */
