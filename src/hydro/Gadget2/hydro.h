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

#include "adiabatic_index.h"

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
    struct part *restrict p, struct xpart *restrict xp) {}

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
 * @param ti_current The current time (on the integer timeline)
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part *restrict p, int ti_current) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float ih = 1.0f / h;
  const float ih2 = ih * ih;
  const float ih4 = ih2 * ih2;

  /* Final operation on the density (add self-contribution). */
  p->rho += p->mass * kernel_root;
  p->rho_dh -= 3.0f * p->mass * kernel_root;
  p->density.wcount += kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  p->rho *= ih * ih2;
  p->rho_dh *= ih4;
  p->density.wcount *= kernel_norm;
  p->density.wcount_dh *= ih * kernel_gamma * kernel_norm;

  const float irho = 1.f / p->rho;

  /* Compute the derivative term */
  p->rho_dh = 1.f / (1.f + 0.33333333f * p->h * p->rho_dh * irho);

  /* Finish calculation of the velocity curl components */
  p->density.rot_v[0] *= ih4 * irho;
  p->density.rot_v[1] *= ih4 * irho;
  p->density.rot_v[2] *= ih4 * irho;

  /* Finish calculation of the velocity divergence */
  p->density.div_v *= ih4 * irho;
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
  const float pressure =
      (p->entropy + p->entropy_dt * half_dt) * pow_gamma(p->rho);

  const float irho = 1.f / p->rho;

  /* Divide the pressure by the density and density gradient */
  const float P_over_rho2 = pressure * irho * irho * p->rho_dh;

  /* Compute the sound speed */
  const float soundspeed = sqrtf(hydro_gamma * pressure * irho);

  /* Compute the Balsara switch */
  const float balsara =
      abs_div_v / (abs_div_v + curl_v + 0.0001f * soundspeed / fac_mu / p->h);

  /* Update variables. */
  p->force.P_over_rho2 = P_over_rho2;
  p->force.soundspeed = soundspeed;
  p->force.balsara = balsara;
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
  const float pressure =
      (p->entropy + p->entropy_dt * dt_entr) * pow_gamma(p->rho);

  const float irho = 1.f / p->rho;

  /* Divide the pressure by the density and density gradient */
  const float P_over_rho2 = pressure * irho * irho * p->rho_dh;

  /* Compute the new sound speed */
  const float soundspeed = sqrtf(hydro_gamma * pressure * irho);

  /* Update variables */
  p->force.P_over_rho2 = P_over_rho2;
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

  p->force.h_dt *= p->h * 0.333333333f;

  p->entropy_dt *= hydro_gamma_minus_one * pow_minus_gamma_minus_one(p->rho);
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

  p->entropy =
      hydro_gamma_minus_one * p->entropy * pow_minus_gamma_minus_one(p->rho);
}

/**
 * @brief Returns the internal energy of a particle
 *
 * @param p The particle of interest
 * @param dt Time since the last kick
 */
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy(
    const struct part *restrict p, float dt) {

  const float entropy = p->entropy + p->entropy_dt * dt;

  return entropy * pow_gamma_minus_one(p->rho) * hydro_one_over_gamma_minus_one;
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
