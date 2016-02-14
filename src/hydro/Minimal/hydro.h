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

#include "approx_math.h"

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 *
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    struct part* p, struct xpart* xp) {

  /* CFL condition */
  const float dt_cfl = 2.f * const_cfl * kernel_gamma * p->h / p->force.v_sig;

  return dt_cfl;
}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the variaous density tasks
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_init_part(struct part* p) {
  p->density.wcount = 0.f;
  p->density.wcount_dh = 0.f;
  p->rho = 0.f;
  p->rho_dh = 0.f;
}

/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 *
 * @param p The particle to act upon
 * @param time The current time
 */
__attribute__((always_inline))
    INLINE static void hydro_end_density(struct part* p, float time) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float ih = 1.0f / h;
  const float ih2 = ih * ih;
  const float ih4 = ih2 * ih2;

  /* Final operation on the density (add self-contribution). */
  p->rho += p->mass * kernel_root;
  p->rho_dh -= 3.0f * p->mass * kernel_root * kernel_igamma;
  p->density.wcount += kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  p->rho *= ih * ih2;
  p->rho_dh *= ih4;
  p->density.wcount *= (4.0f / 3.0f * M_PI * kernel_gamma3);
  p->density.wcount_dh *= ih * (4.0f / 3.0f * M_PI * kernel_gamma3);
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
    struct part* p, struct xpart* xp, float time) {

  p->force.pressure = p->rho * p->u * (const_hydro_gamma - 1.f);
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking place in the variaous force tasks
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_reset_acceleration(struct part* p) {

  /* Reset the acceleration. */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* Reset the time derivatives. */
  p->u_dt = 0.0f;
  p->h_dt = 0.0f;
  p->force.v_sig = 0.0f;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param p The particle
 * @param xp The extended data of the particle
 * @param t0 The time at the start of the drift
 * @param t1 The time at the end of the drift
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part* p, struct xpart* xp, float t0, float t1) {

  const float dt = t1 - t0;

  /* Predict internal energy */
  const float w = p->u_dt / p->u * dt;
  if (fabsf(w) < 0.2f)
    p->u *= approx_expf(w); /* 4th order expansion of exp(w) */
  else
    p->u *= expf(w);

  /* Need to recompute the pressure as well */
  p->force.pressure = p->rho * p->u * (const_hydro_gamma - 1.f);
}

/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the forces and accelerationsby the appropiate constants
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_end_force(struct part* p) {}

/**
 * @brief Kick the additional variables
 *
 * @param p The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline))
    INLINE static void hydro_kick_extra(struct part* p, float dt) {}

/**
 * @brief Converts hydro quantity of a particle
 *
 * Requires the density to be known
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline))
    INLINE static void hydro_convert_quantities(struct part* p) {}
