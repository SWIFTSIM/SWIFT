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
#include "cosmology.h"
#include "dimension.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "hydro_properties.h"
#include "hydro_space.h"
#include "kernel_hydro.h"
#include "minmax.h"

#include "./hydro_parameters.h"

/**
 * @brief Returns the comoving internal energy of a particle at the last
 * time the particle was kicked.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy(const struct part *restrict p,
                                   const struct xpart *restrict xp) {

  return gas_internal_energy_from_entropy(p->rho_bar, xp->entropy_full);
}

/**
 * @brief Returns the physical internal energy of a particle at the last
 * time the particle was kicked.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy(const struct part *restrict p,
                                   const struct xpart *restrict xp,
                                   const struct cosmology *cosmo) {

  return gas_internal_energy_from_entropy(p->rho_bar * cosmo->a3_inv,
                                          xp->entropy_full);
}
/**
 * @brief Returns the comoving internal energy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_internal_energy(const struct part *restrict p) {

  return gas_internal_energy_from_entropy(p->rho_bar, p->entropy);
}

/**
 * @brief Returns the physical internal energy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_physical_internal_energy(const struct part *restrict p,
                                           const struct cosmology *cosmo) {

  return gas_internal_energy_from_entropy(p->rho_bar * cosmo->a3_inv,
                                          p->entropy);
}

/**
 * @brief Returns the comoving pressure of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_pressure(
    const struct part *restrict p) {

  return gas_pressure_from_entropy(p->rho_bar, p->entropy);
}

/**
 * @brief Returns the physical pressure of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_pressure(
    const struct part *restrict p, const struct cosmology *cosmo) {

  return gas_pressure_from_entropy(p->rho_bar * cosmo->a3_inv, p->entropy);
}

/**
 * @brief Returns the comoving entropy of a particle at the last
 * time the particle was kicked.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_entropy(
    const struct part *restrict p, const struct xpart *restrict xp) {

  return xp->entropy_full;
}

/**
 * @brief Returns the physical entropy of a particl at the last
 * time the particle was kicked.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_entropy(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct cosmology *cosmo) {

  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return xp->entropy_full;
}

/**
 * @brief Sets the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param u The physical entropy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy(struct part *p, struct xpart *xp,
                                   const struct cosmology *cosmo,
                                   const float u) {
  error("To be implemented");
}

/**
 * @brief Sets the drifted physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param u The physical entropy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(struct part *p,
                                           const struct cosmology *cosmo,
                                           const float u) {
  error("To be implemented");
}

/**
 * @brief Update the value of the viscosity alpha for the scheme.
 *
 * @param p the particle of interest
 * @param alpha the new value for the viscosity coefficient.
 */
__attribute__((always_inline)) INLINE static void hydro_set_viscosity_alpha(
    struct part *restrict p, float alpha) {
  /* This scheme has fixed alpha */
}

/**
 * @brief Update the value of the viscosity alpha to the
 *        feedback reset value for the scheme.
 *
 * @param p the particle of interest
 */
__attribute__((always_inline)) INLINE static void
hydro_set_viscosity_alpha_max_feedback(struct part *restrict p) {
  /* This scheme has fixed alpha */
}

/**
 * @brief Returns the comoving entropy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_entropy(const struct part *restrict p) {

  return p->entropy;
}

/**
 * @brief Returns the physical entropy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_physical_entropy(const struct part *restrict p,
                                   const struct cosmology *cosmo) {

  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return p->entropy;
}

/**
 * @brief Returns the comoving sound speed of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_soundspeed(const struct part *restrict p) {

  return p->force.soundspeed;
}

/**
 * @brief Returns the physical sound speed of a particle
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_soundspeed(const struct part *restrict p,
                              const struct cosmology *cosmo) {

  return cosmo->a_factor_sound_speed * p->force.soundspeed;
}

/**
 * @brief Returns the comoving density of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_density(
    const struct part *restrict p) {

  return p->rho;
}

/**
 * @brief Returns the physical density of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_density(
    const struct part *restrict p, const struct cosmology *cosmo) {

  return p->rho * cosmo->a3_inv;
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
 * @brief Sets the mass of a particle
 *
 * @param p The particle of interest
 * @param m The mass to set.
 */
__attribute__((always_inline)) INLINE static void hydro_set_mass(
    struct part *restrict p, float m) {

  p->mass = m;
}

/**
 * @brief Returns the velocities drifted to the current time of a particle.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle.
 * @param dt_kick_hydro The time (for hydro accelerations) since the last kick.
 * @param dt_kick_grav The time (for gravity accelerations) since the last kick.
 * @param v (return) The velocities at the current time.
 */
__attribute__((always_inline)) INLINE static void hydro_get_drifted_velocities(
    const struct part *restrict p, const struct xpart *xp, float dt_kick_hydro,
    float dt_kick_grav, float v[3]) {

  v[0] = xp->v_full[0] + p->a_hydro[0] * dt_kick_hydro +
         xp->a_grav[0] * dt_kick_grav;
  v[1] = xp->v_full[1] + p->a_hydro[1] * dt_kick_hydro +
         xp->a_grav[1] * dt_kick_grav;
  v[2] = xp->v_full[2] + p->a_hydro[2] * dt_kick_hydro +
         xp->a_grav[2] * dt_kick_grav;
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy_dt(const struct part *restrict p) {

  return gas_internal_energy_from_entropy(p->rho_bar, p->entropy_dt);
}

/**
 * @brief Returns the time derivative of physical internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy_dt(const struct part *restrict p,
                                      const struct cosmology *cosmo) {

  return gas_internal_energy_from_entropy(p->rho_bar * cosmo->a3_inv,
                                          p->entropy_dt);
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest.
 * @param du_dt The new time derivative of the internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_comoving_internal_energy_dt(struct part *restrict p, float du_dt) {

  p->entropy_dt = gas_entropy_from_internal_energy(p->rho_bar, du_dt);
}

/**
 * @brief Sets the time derivative of the physical internal energy of a particle
 *
 * We assume a constant density for the conversion to entropy.
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param du_dt The time derivative of the internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy_dt(struct part *restrict p,
                                      const struct cosmology *restrict cosmo,
                                      float du_dt) {
  p->entropy_dt =
      gas_entropy_from_internal_energy(p->rho_bar * cosmo->a3_inv, du_dt);
}
/**
 * @brief Sets the physical entropy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param entropy The physical entropy
 */
__attribute__((always_inline)) INLINE static void hydro_set_physical_entropy(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const float entropy) {

  /* Note there is no conversion from physical to comoving entropy */
  xp->entropy_full = entropy;
}

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param hydro_properties The constants used in the scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct hydro_props *restrict hydro_properties,
    const struct cosmology *restrict cosmo) {

  const float CFL = hydro_properties->CFL_condition;

  /* CFL condition */
  const float dt_cfl = 2.f * kernel_gamma * CFL * cosmo->a * p->h /
                       (cosmo->a_factor_sound_speed * p->force.v_sig);

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
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part *restrict p, const struct cosmology *cosmo) {

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
  p->density.wcount_dh -= hydro_dimension * kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  p->rho *= h_inv_dim;
  p->rho_bar *= h_inv_dim;
  p->density.rho_dh *= h_inv_dim_plus_one;
  p->density.pressure_dh *= h_inv_dim_plus_one;
  p->density.wcount *= h_inv_dim;
  p->density.wcount_dh *= h_inv_dim_plus_one;

  const float rho_inv = 1.f / p->rho;
  const float a_inv2 = cosmo->a2_inv;
  const float entropy_minus_one_over_gamma = 1.f / p->entropy_one_over_gamma;

  /* Final operation on the weighted density */
  p->rho_bar *= entropy_minus_one_over_gamma;

  /* Finish calculation of the velocity curl components */
  p->density.rot_v[0] *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  p->density.rot_v[1] *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  p->density.rot_v[2] *= h_inv_dim_plus_one * a_inv2 * rho_inv;

  /* Finish calculation of the velocity divergence */
  p->density.div_v *= h_inv_dim_plus_one * rho_inv;
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_part_has_no_neighbours(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  p->rho = p->mass * kernel_root * h_inv_dim;
  p->rho_bar = p->mass * kernel_root * h_inv_dim;
  p->density.wcount = kernel_root * h_inv_dim;
  p->density.rho_dh = 0.f;
  p->density.wcount_dh = 0.f;
  p->density.pressure_dh = 0.f;
  p->density.div_v = 0.f;
  p->density.rot_v[0] = 0.f;
  p->density.rot_v[1] = 0.f;
  p->density.rot_v[2] = 0.f;
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
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const float dt_alpha) {

  const float fac_mu = cosmo->a_factor_mu;

  /* Compute the norm of the curl */
  const float curl_v = sqrtf(p->density.rot_v[0] * p->density.rot_v[0] +
                             p->density.rot_v[1] * p->density.rot_v[1] +
                             p->density.rot_v[2] * p->density.rot_v[2]);

  /* Compute the norm of div v */
  const float abs_div_v = fabsf(p->density.div_v);

  /* Compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho_bar, p->entropy);

  /* Compute the sound speed from the pressure*/
  const float soundspeed = gas_soundspeed_from_pressure(p->rho_bar, pressure);

  /* Compute the Balsara switch */
  const float balsara =
      hydro_props->viscosity.alpha * abs_div_v /
      (abs_div_v + curl_v + 0.0001f * soundspeed * fac_mu / p->h);

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
  p->entropy_one_over_gamma = pow_one_over_gamma(p->entropy);
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param p The particle
 * @param xp The extended data of the particle
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part *restrict p, const struct xpart *restrict xp, float dt_drift,
    float dt_therm, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Predict the entropy */
  p->entropy += p->entropy_dt * dt_therm;

  /* Check against entropy floor */
  const float floor_A = entropy_floor(p, cosmo, floor_props);

  /* Check against absolute minimum */
  const float min_u = hydro_props->minimal_internal_energy;
  const float min_A =
      gas_entropy_from_internal_energy(p->rho * cosmo->a3_inv, min_u);

  p->entropy = max(p->entropy, floor_A);
  p->entropy = max(p->entropy, min_A);

  const float h_inv = 1.f / p->h;

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt_drift;
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

  /* Compute the pressure */
  const float pressure = gas_pressure_from_entropy(p->rho_bar, p->entropy);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho_bar, pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_bar_inv = 1.f / p->rho_bar;
  const float P_over_rho2 = pressure * rho_bar_inv * rho_bar_inv;

  /* Update the variables */
  p->entropy_one_over_gamma = pow_one_over_gamma(p->entropy);
  p->force.soundspeed = soundspeed;
  p->force.P_over_rho2 = P_over_rho2;
}

/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the forces and accelerationsby the appropiate constants
 *
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_force(
    struct part *restrict p, const struct cosmology *cosmo) {

  p->force.h_dt *= p->h * hydro_dimension_inv;

  p->entropy_dt = 0.5f * cosmo->a2_inv *
                  gas_entropy_from_internal_energy(p->rho_bar, p->entropy_dt);
}

/**
 * @brief Kick the additional variables
 *
 * @param p The particle to act upon
 * @param xp The particle extended data to act upon
 * @param dt_therm The time-step for this kick (for thermodynamic quantities)
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part *restrict p, struct xpart *restrict xp, float dt_therm,
    float dt_grav, float dt_hydro, float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Integrate the entropy forward in time */
  const float delta_entropy = p->entropy_dt * dt_therm;

  /* Do not decrease the entropy by more than a factor of 2 */
  xp->entropy_full =
      max(xp->entropy_full + delta_entropy, 0.5f * xp->entropy_full);

  /* Check against entropy floor */
  const float floor_A = entropy_floor(p, cosmo, floor_props);

  /* Check against absolute minimum */
  const float min_u = hydro_props->minimal_internal_energy;
  const float min_A =
      gas_entropy_from_internal_energy(p->rho * cosmo->a3_inv, min_u);

  /* Take highest of both limits */
  const float entropy_min = max(min_A, floor_A);

  if (xp->entropy_full < entropy_min) {
    xp->entropy_full = entropy_min;
    p->entropy_dt = 0.f;
  }
}

/**
 *  @brief Converts hydro quantity of a particle at the start of a run
 *
 * Requires the density to be known
 *
 * @param p The particle to act upon
 * @param xp The extended data.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props) {

  /* We read u in the entropy field. We now get S from u */
  xp->entropy_full =
      gas_entropy_from_internal_energy(p->rho_bar * cosmo->a3_inv, p->entropy);
  p->entropy = xp->entropy_full;
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
  p->wakeup = time_bin_not_awake;
  p->rho_bar = 0.f;
  p->entropy_one_over_gamma = pow_one_over_gamma(p->entropy);
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
  xp->a_grav[0] = 0.f;
  xp->a_grav[1] = 0.f;
  xp->a_grav[2] = 0.f;

  hydro_reset_acceleration(p);
  hydro_init_part(p, NULL);
}

/**
 * @brief Overwrite the initial internal energy of a particle.
 *
 * Note that in the cases where the thermodynamic variable is not
 * internal energy but gets converted later, we must overwrite that
 * field. The conversion to the actual variable happens later after
 * the initial fake time-step.
 *
 * @param p The #part to write to.
 * @param u_init The new initial internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_init_internal_energy(struct part *p, float u_init) {

  p->entropy = u_init;
}

/**
 * @brief Operations performed when a particle gets removed from the
 * simulation volume.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
__attribute__((always_inline)) INLINE static void hydro_remove_part(
    const struct part *p, const struct xpart *xp) {}

#endif /* SWIFT_PRESSURE_ENTROPY_HYDRO_H */
