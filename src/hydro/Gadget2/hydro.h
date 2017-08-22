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
#ifndef SWIFT_GADGET2_HYDRO_H
#define SWIFT_GADGET2_HYDRO_H

/**
 * @file Gadget2/hydro.h
 * @brief SPH interaction functions following the Gadget-2 version of SPH.
 *
 * The interactions computed here are the ones presented in the Gadget-2 paper
 * Springel, V., MNRAS, Volume 364, Issue 4, pp. 1105-1134.
 * We use the same numerical coefficients as the Gadget-2 code. When used with
 * the Spline-3 kernel, the results should be equivalent to the ones obtained
 * with Gadget-2 up to the rounding errors and interactions missed by the
 * Gadget-2 tree-code neighbours search.
 */

#include "adiabatic_index.h"
#include "approx_math.h"
#include "dimension.h"
#include "equation_of_state.h"
#include "hydro_properties.h"
#include "hydro_space.h"
#include "kernel_hydro.h"
#include "minmax.h"

/**
 * @brief Returns the internal energy of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_internal_energy(
    const struct part *restrict p) {

  return gas_internal_energy_from_entropy(p->rho, p->entropy);
}

/**
 * @brief Returns the pressure of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_pressure(
    const struct part *restrict p) {

  return gas_pressure_from_entropy(p->rho, p->entropy);
}

/**
 * @brief Returns the entropy of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_entropy(
    const struct part *restrict p) {

  return p->entropy;
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

  return gas_internal_energy_from_entropy(p->rho, p->entropy_dt);
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

  p->entropy_dt = gas_entropy_from_internal_energy(p->rho, du_dt);
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

  p->rho = 0.f;
  p->density.wcount = 0.f;
  p->density.wcount_dh = 0.f;
  p->density.rho_dh = 0.f;
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
  p->density.rho_dh -= hydro_dimension * p->mass * kernel_root;
  p->density.wcount += kernel_root;
  p->density.wcount_dh -= hydro_dimension * kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  p->rho *= h_inv_dim;
  p->density.rho_dh *= h_inv_dim_plus_one;
  p->density.wcount *= h_inv_dim;
  p->density.wcount_dh *= h_inv_dim_plus_one;

  const float rho_inv = 1.f / p->rho;

  /* Finish calculation of the velocity curl components */
  p->density.rot_v[0] *= h_inv_dim_plus_one * rho_inv;
  p->density.rot_v[1] *= h_inv_dim_plus_one * rho_inv;
  p->density.rot_v[2] *= h_inv_dim_plus_one * rho_inv;

  /* Finish calculation of the velocity divergence */
  p->density.div_v *= h_inv_dim_plus_one * rho_inv;
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
  p->density.rho_dh = 0.f;
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
 * @param ti_current The current time (on the timeline)
 * @param timeBase The minimal time-step size
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part *restrict p, struct xpart *restrict xp) {

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  /* Inverse of the physical density */
  const float rho_inv = 1.f / p->rho;

  /* Compute the norm of the curl */
  const float curl_v = sqrtf(p->density.rot_v[0] * p->density.rot_v[0] +
                             p->density.rot_v[1] * p->density.rot_v[1] +
                             p->density.rot_v[2] * p->density.rot_v[2]);

  /* Compute the norm of div v */
  const float abs_div_v = fabsf(p->density.div_v);

  /* Compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho, p->entropy);

  /* Compute the sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float P_over_rho2 = pressure * rho_inv * rho_inv;

  /* Compute the Balsara switch */
  const float balsara =
      abs_div_v / (abs_div_v + curl_v + 0.0001f * soundspeed / fac_mu / p->h);

  /* Compute the "grad h" term */
  const float omega_inv =
      1.f / (1.f + hydro_dimension_inv * p->h * p->density.rho_dh * rho_inv);

  /* Update variables. */
  p->force.f = omega_inv;
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
  p->force.v_sig = p->force.soundspeed;
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

  /* Re-set the entropy */
  p->entropy = xp->entropy_full;
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
  if (fabsf(w2) < 0.2f)
    p->rho *= approx_expf(w2); /* 4th order expansion of exp(w) */
  else
    p->rho *= expf(w2);

  /* Predict the entropy */
  p->entropy += p->entropy_dt * dt;

  /* Re-compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho, p->entropy);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_inv = 1.f / p->rho;
  const float P_over_rho2 = pressure * rho_inv * rho_inv;

  /* Update variables */
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
      0.5f * gas_entropy_from_internal_energy(p->rho, p->entropy_dt);
}

/**
 * @brief Kick the additional variables
 *
 * @param p The particle to act upon
 * @param xp The particle extended data to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part *restrict p, struct xpart *restrict xp, float dt) {

  /* Do not decrease the entropy by more than a factor of 2 */
  if (dt > 0. && p->entropy_dt * dt < -0.5f * xp->entropy_full) {
    /* message("Warning! Limiting entropy_dt. Possible cooling error.\n
     * entropy_full = %g \n entropy_dt * dt =%g \n", */
    /* 	     xp->entropy_full,p->entropy_dt * dt); */
    p->entropy_dt = -0.5f * xp->entropy_full / dt;
  }
  xp->entropy_full += p->entropy_dt * dt;

  /* Compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho, xp->entropy_full);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_inv = 1.f / p->rho;
  const float P_over_rho2 = pressure * rho_inv * rho_inv;

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
    struct part *restrict p, struct xpart *restrict xp) {

  /* We read u in the entropy field. We now get S from u */
  xp->entropy_full = gas_entropy_from_internal_energy(p->rho, p->entropy);
  p->entropy = xp->entropy_full;

  /* Compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho, p->entropy);

  /* Compute the sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_inv = 1.f / p->rho;
  const float P_over_rho2 = pressure * rho_inv * rho_inv;

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
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
  xp->entropy_full = p->entropy;

  hydro_reset_acceleration(p);
  hydro_init_part(p, NULL);
}

#endif /* SWIFT_GADGET2_HYDRO_H */
