/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl) &
 *                    Josh Borrow (joshua.borrow@durham.ac.uk)
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
#ifndef SWIFT_ANARCHY_PU_HYDRO_H
#define SWIFT_ANARCHY_PU_HYDRO_H

/**
 * @file AnarchyPU/hydro.h
 * @brief P-U conservative implementation of SPH,
 *        with added ANARCHY physics (Cullen & Denhen 2011 AV,
 *        Price 2008 thermal diffusion (Non-neighbour loop
 *        equations)
 *
 * This implementation corresponds to the one presented in the SWIFT
 * documentation and in Hopkins, "A general class of Lagrangian smoothed
 * particle hydrodynamics methods and implications for fluid mixing problems",
 * MNRAS, 2013.
 */

#include "adiabatic_index.h"
#include "approx_math.h"
#include "cosmology.h"
#include "dimension.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "hydro_properties.h"
#include "hydro_space.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "pressure_floor.h"

#include <float.h>

/**
 * @brief Returns the comoving internal energy of a particle at the last
 * time the particle was kicked.
 *
 * For implementations where the main thermodynamic variable
 * is not internal energy, this function computes the internal
 * energy from the thermodynamic variable.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy(size_t pind,
                                   ) {

  return part_get_u_full(pind);
}

/**
 * @brief Returns the physical internal energy of a particle at the last
 * time the particle was kicked.
 *
 * For implementations where the main thermodynamic variable
 * is not internal energy, this function computes the internal
 * energy from the thermodynamic variable and converts it to
 * physical coordinates.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy(size_t pind,
                                   const struct cosmology *cosmo) {

  return part_get_u_full(pind) * cosmo->a_factor_internal_energy;
}

/**
 * @brief Returns the comoving internal energy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_internal_energy(size_t pind) {

  return part_get_u(pind);
}

/**
 * @brief Returns the physical internal energy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_physical_internal_energy(size_t pind,
                                           const struct cosmology *cosmo) {

  return part_get_u(pind) * cosmo->a_factor_internal_energy;
}

/**
 * @brief Returns the comoving pressure of a particle
 *
 * Computes the pressure based on the particle's properties.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_pressure(
    size_t pind) {

  return part_get_pressure_bar(pind);
}

/**
 * @brief Returns the physical pressure of a particle
 *
 * Computes the pressure based on the particle's properties and
 * convert it to physical coordinates.
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_pressure(
    size_t pind, const struct cosmology *cosmo) {

  return cosmo->a_factor_pressure * part_get_pressure_bar(pind);
}

/**
 * @brief Returns the comoving entropy of a particle at the last
 * time the particle was kicked.
 *
 * For implementations where the main thermodynamic variable
 * is not entropy, this function computes the entropy from
 * the thermodynamic variable.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_entropy(
    size_t pind) {

  return gas_entropy_from_internal_energy(part_get_rho(pind), part_get_u_full(pind));
}

/**
 * @brief Returns the physical entropy of a particle at the last
 * time the particle was kicked.
 *
 * For implementations where the main thermodynamic variable
 * is not entropy, this function computes the entropy from
 * the thermodynamic variable and converts it to
 * physical coordinates.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_entropy(
    size_t pind,
    const struct cosmology *cosmo) {

  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return gas_entropy_from_internal_energy(part_get_rho(pind), part_get_u_full(pind));
}

/**
 * @brief Returns the comoving entropy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_entropy(size_t pind) {

  return gas_entropy_from_internal_energy(part_get_rho(pind), part_get_u(pind));
}

/**
 * @brief Returns the physical entropy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_physical_entropy(size_t pind,
                                   const struct cosmology *cosmo) {

  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return gas_entropy_from_internal_energy(part_get_rho(pind), part_get_u(pind));
}

/**
 * @brief Returns the comoving sound speed of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_soundspeed(size_t pind) {

  /* Compute the sound speed -- see theory section for justification */
  /* IDEAL GAS ONLY -- P-U does not work with generic EoS. */

  return gas_soundspeed_from_pressure(part_get_rho(pind), part_get_pressure_bar(pind));
}

/**
 * @brief Returns the physical sound speed of a particle
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_soundspeed(size_t pind,
                              const struct cosmology *cosmo) {

  return cosmo->a_factor_sound_speed * hydro_get_comoving_soundspeed(p);
}

/**
 * @brief Returns the comoving density of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_density(
    size_t pind) {

  return part_get_rho(pind);
}

/**
 * @brief Returns the comoving density of a particle.
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_density(
    size_t pind, const struct cosmology *cosmo) {

  return cosmo->a3_inv * part_get_rho(pind);
}

/**
 * @brief Returns the mass of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_mass(
    size_t pind) {

  return part_get_mass(pind);
}

/**
 * @brief Sets the mass of a particle
 *
 * @param p The particle of interest
 * @param m The mass to set.
 */
__attribute__((always_inline)) INLINE static void hydro_set_mass(
    size_t pind, float m) {

  part_set_mass(pind, m);
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy_dt(size_t pind) {

  return part_get_u_dt(pind);
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 * @param cosmo Cosmology data structure
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy_dt(size_t pind,
                                      const struct cosmology *cosmo) {

  return part_get_u_dt(pind) * cosmo->a_factor_internal_energy;
}

/**
 * @brief Sets the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest.
 * @param du_dt The new time derivative of the internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_comoving_internal_energy_dt(size_t pind, float du_dt) {

  part_set_u_dt(pind, du_dt);
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param du_dt The new time derivative of the internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy_dt(size_t pind,
                                      const struct cosmology *cosmo,
                                      float du_dt) {

  part_set_u_dt(pind, du_dt / cosmo->a_factor_internal_energy);
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
    size_t pind, const struct cosmology *cosmo,
    const float entropy) {

  /* Note there is no conversion from physical to comoving entropy */
  const float comoving_entropy = entropy;
  float u_full = gas_internal_energy_from_entropy(part_get_rho(pind), comoving_entropy);
  part_set_u_full(pind, u_full);
}

/**
 * @brief Sets the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy(size_t pind,
                                   const struct cosmology *cosmo,
                                   const float u) {

  part_set_u_full(pind, u / cosmo->a_factor_internal_energy);
}

/**
 * @brief Sets the drifted physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(size_t pind,
                                           const struct cosmology *cosmo,
                                           const float u) {

  /* Store ratio of new internal energy to old internal energy, as we use this
   * in the drifting of the pressure. */
  float internal_energy_ratio = 1.f / part_get_u(pind);

  /* Update the internal energy */
  part_set_u(pind, u / cosmo->a_factor_internal_energy);
  internal_energy_ratio *= part_get_u(pind);

  /* Now we can use this to 'update' the value of the smoothed pressure. To
   * truly update this variable, we would need another loop over neighbours
   * using the new internal energies of everyone, but that's not feasible. */

  pressure_bar *= internal_energy_ratio;
  part_set_pressure_bar(pind, pressure_bar);

  /* Now recompute the extra quantities */

  /* Compute the sound speed */
  const float soundspeed =
    gas_soundspeed_from_pressure(part_get_rho(pind), part_get_pressure_bar(pind));

  /* Update variables. */
  part_set_soundspeed(pind, soundspeed);

  part_set_v_sig(pind, max(part_get_v_sig(pind), 2.f * soundspeed));
}

/**
 * @brief Correct the signal velocity of the particle partaking in
 * supernova (kinetic) feedback based on the velocity kick the particle receives
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param dv_phys The velocity kick received by the particle expressed in
 * physical units (note that dv_phys must be positive or equal to zero)
 */
__attribute__((always_inline)) INLINE static void
hydro_set_v_sig_based_on_velocity_kick(size_t pind,
                                       const struct cosmology *cosmo,
                                       const float dv_phys) {

  /* Compute the velocity kick in comoving coordinates */
  const float dv = dv_phys / cosmo->a_factor_sound_speed;

  /* Sound speed */
  const float soundspeed = hydro_get_comoving_soundspeed(pind);

  /* Update the signal velocity */
  part_set_v_sig(
		 pind, max(2.f * soundspeed, part_get_v_sig(pind) + const_viscosity_beta * dv);
}

/**
 * @brief Update the value of the viscosity alpha for the scheme.
 *
 * @param p the particle of interest
 * @param alpha the new value for the viscosity coefficient.
 */
__attribute__((always_inline)) INLINE static void hydro_set_viscosity_alpha(
    size_t pind, float alpha) {
    part_set_alpha_av(pind, alpha);
}

/**
 * @brief Update the value of the diffusive coefficients to the
 *        feedback reset value for the scheme.
 *
 * @param p the particle of interest
 */
__attribute__((always_inline)) INLINE static void
hydro_diffusive_feedback_reset(size_t pind) {
  /* Set the viscosity to the max, and the diffusion to the min */
  hydro_set_viscosity_alpha(pind,
                            hydro_props_default_viscosity_alpha_feedback_reset);

  part_set_alpha_diff(pind, hydro_props_default_diffusion_alpha_feedback_reset);
}

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * This function returns the time-step of a particle given its hydro-dynamical
 * state. A typical time-step calculation would be the use of the CFL condition.
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 * @param hydro_properties The SPH parameters
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    size_t pind,
    const struct hydro_props *restrict hydro_properties,
    const struct cosmology *restrict cosmo) {

  const float CFL_condition = hydro_properties->CFL_condition;

  /* CFL condition */
  const float dt_cfl = 2.f * kernel_gamma * CFL_condition * cosmo->a * part_get_h(pind) /
    (cosmo->a_factor_sound_speed * part_get_v_sig(pind));

  return dt_cfl;
}

/**
 * @brief Compute the signal velocity between two gas particles
 *
 * This is eq. (103) of Price D., JCoPh, 2012, Vol. 231, Issue 3.
 *
 * @param dx Comoving vector separating both particles (pi - pj).
 * @brief pi The first #part.
 * @brief pj The second #part.
 * @brief mu_ij The velocity on the axis linking the particles, or zero if the
 * particles are moving away from each other,
 * @brief beta The non-linear viscosity constant.
 */
__attribute__((always_inline)) INLINE static float hydro_signal_velocity(
    const float dx[3], size_t pindi,
    size_t pindj, const float mu_ij, const float beta) {

  const float ci = part_get_soundspeed(pindi);
  const float cj = part_get_soundspeed(pindj);

  return ci + cj - beta * mu_ij;
}

/**
 * @brief Does some extra hydro operations once the actual physical time step
 * for the particle is known.
 *
 * @param p The particle to act upon.
 * @param dt Physical time step of the particle during the next step.
 */
__attribute__((always_inline)) INLINE static void hydro_timestep_extra(
    struct size_t pind, float dt) {}

/**
 * @brief Operations performed when a particle gets removed from the
 * simulation volume.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @param time The simulation time.
 */
__attribute__((always_inline)) INLINE static void hydro_remove_part(
    size_t pind, const double time) {}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various density loop over neighbours. Typically, all fields of the
 * density sub-structure of a particle get zeroed in here.
 *
 * @param p The particle to act upon
 * @param hs #hydro_space containing hydro specific space information.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    size_t pind, const struct hydro_space *hs) {
  
  part_set_wcount(pind, 0.f);
  part_set_wcount_dh(pind, 0.f);
  part_set_rho(pind, 0.f);
  part_set_rho_dh(pind, 0.f);
  part_set_pressure_bar(pind, 0.f);
  part_set_pressure_bar_dh(pind, 0.f);

  part_set_rot_v_ind(pind, 0, 0.f);
  part_set_rot_v_ind(pind, 1, 0.f);
  part_set_rot_v_ind(pind, 2, 0.f);

  part_set_div_v(pind, 0.f);
  part_set_laplace_u(pind, 0.f);

}

/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 * Additional quantities such as velocity gradients will also get the final
 * terms added to them here.
 *
 * Also adds/multiplies the cosmological terms if need be.
 *
 * @param p The particle to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    size_t pind, const struct cosmology *cosmo) {

  /* Some smoothing length multiples. */
  const float h = part_get_h_(pind);
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Final operation on the density (add self-contribution). */
  const float m = part_get_mass(pind);

  float rho = part_get_rho(pind);
  rho += m * kernel_root;

  float rho_dh = part_get_rho_dh(pind);
  rho_dh -= hydro_dimension * m * kernel_root;

  float pressure_bar = part_get_pressure_bar(pind);
  pressure_bar += m * part_get_u(pind) * kernel_root;

  float pressure_bar_dh = part_get_pressure_bar_dh(pind);
  pressure_bar_dh -= hydro_dimension * mass * part_get_u(pind) * kernel_root;

  float wcount = part_get_wcount(pind);
  wcount += kernel_root;

  float wcount_dh = part_get_wcount_dh(pind);
  wcount_dh -= hydro_dimension * kernel_root;

 
  /* Finish the calculation by inserting the missing h-factors */
  rho *= h_inv_dim;
  rho_dh *= h_inv_dim_plus_one;
  pressure_bar *= (h_inv_dim * hydro_gamma_minus_one);
  pressure_bar_dh *= (h_inv_dim_plus_one * hydro_gamma_minus_one);
  wcount *= h_inv_dim;
  wcount_dh *= h_inv_dim_plus_one;

  part_set_rho(pind, rho);
  part_set_rho_dh(pind, rho_pd);
  part_set_pressure_bar(pind, pressure_bar);
  part_set_pressure_bar_dh(pind, pressure_bar_dh);
  part_set_wcount(pind, wcount);
  part_set_wcount_dh(pind, wcount_dh);

  const float rho_inv = 1.f / rho;
  const float a_inv2 = cosmo->a2_inv;

  /* Finish calculation of the velocity curl components */
  float *rot_v = part_get_rot_v(pind);
  rot_v[0] *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  rot_v[1] *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  rot_v[2] *= h_inv_dim_plus_one * a_inv2 * rho_inv;

  /* Finish calculation of the velocity divergence */
  float div_v = part_get_div_v(pind);
  div_v *= h_inv_dim_plus_one * rho_inv * a_inv2;
  div_v += cosmo->H * hydro_dimension;
  part_set_div_v(pind, div_v);
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
    size_t pind,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

  const float fac_B = cosmo->a_factor_Balsara_eps;

  /* Compute the norm of the curl */
  const float *rot_v = part_get_rot_v(pind);
  const float curl_v = sqrtf(rot_v[0] * rot_v[0] +
                             rot_v[1] * rot_v[1] +
                             rot_v[2] * rot_v[2]);

  /* Compute the norm of div v */
  const float abs_div_v = fabsf(part_get_div_v(pind));

  /* Compute the sound speed -- see theory section for justification */
  const float soundspeed = hydro_get_comoving_soundspeed(pind);

  /* Compute the Balsara switch */
  const float h = part_get_h(pind);
  const float wcount = part_get_wcount(pind);
  const float wcount_dh = part_get_wcount_dh(pind);
  const float balsara =
      abs_div_v / (abs_div_v + curl_v + 0.0001f * soundspeed * fac_B / h);

  /* Compute the "grad h" term */
  const float common_factor = h / (hydro_dimension * wcount);
  float grad_h_term;

  /* Ignore changing-kernel effects when h ~= h_max */
  if (h > 0.9999f * hydro_props->h_max) {
    grad_h_term = 0.f;
    warning("h ~ h_max for particle with ID %lld (h: %g)", part_get_id(pind), h);
  } else {
    const float grad_W_term = common_factor * wcount_dh;
    if (grad_W_term < -0.9999f) {
      /* if we get here, we either had very small neighbour contributions
         (which should be treated as a no neighbour case in the ghost) or
         a very weird particle distribution (e.g. particles sitting on
         top of each other). Either way, we cannot use the normal
         expression, since that would lead to overflow or excessive round
         off and cause excessively high accelerations in the force loop */
      grad_h_term = 0.f;
      warning(
          "grad_W_term very small for particle with ID %lld (h: %g, wcount: "
          "%g, wcount_dh: %g).",
          part_get_id(pind), h, wcount, wcount_dh);
    } else {
      grad_h_term = (part_get_pressure_bar_dh(pind) * common_factor *
                     hydro_one_over_gamma_minus_one) /
                    (1.f + grad_W_term);
    }
  }

  /* Update variables. */
  part_set_f_gradh(pind, grad_h_term);
  part_set_soundspeed(pind, soundspeed);
  part_set_balsara(pind, balsara);
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
    size_t pind) {

  part_set_v_sig(pind, 2.f * part_get_soundspeed(pind));
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
    size_t pind) {

  /* Some smoothing length multiples. */
  const float h = part_get_h(pind);
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Include the extra factors in the del^2 u */

  part_set_laplace_u(pind, part_get_laplace_u(pind) * 2.f * h_inv_dim_plus_one);
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * In the desperate case where a particle has no neighbours (likely because
 * of the h_max ceiling), set the particle fields to something sensible to avoid
 * NaNs in the next calculations.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_part_has_no_neighbours(
    size_t pind,
    const struct cosmology *cosmo) {

  /* Some smoothing length multiples. */
  const float h = part_get_h(pind);
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  warning(
      "Gas particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      part_get_id(pind), h, part_get_wcount(pind));

  /* Re-set problematic values */
  part_set_rho(pind, part_get_mass(pind) * kernel_root * h_inv_dim);
  part_set_v_sig(pind, 0.f);
  part_set_pressure_bar(pind, part_get_mass(pind) * part_get_u(pind) * hydro_gamma_minus_one * kernel_root * h_inv_dim);
  part_set_wcount(pind, kernel_root * h_inv_dim);
  part_set_rho_dh(pind, 0.f);
  part_set_wcount_dh(pind, 0.f);
  part_set_pressure_bar_dh(pind, 0.f);

  part_set_rot_v_ind(pind, 0, 0.f);
  part_set_rot_v_ind(pind, 1, 0.f);
  part_set_rot_v_ind(pind, 2, 0.f);
  
  /* Probably not shocking, so this is safe to do */
  part_set_div_v(pind, 0.f);
  part_set_laplace_u(pind, 0.f);
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
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    size_t pind,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor, const float dt_alpha,
    const float dt_therm) {

  /* Here we need to update the artificial viscosity */

  /* We use in this function that h is the radius of support */
  const float kernel_support_physical = part_get_h(pind) * cosmo->a * kernel_gamma;
  const float kernel_support_physical_inv = 1.f / kernel_support_physical;
  const float v_sig_physical = part_get_v_sig(pind) * cosmo->a_factor_sound_speed;
  const float soundspeed_physical = hydro_get_physical_soundspeed(pind, cosmo);

  const float sound_crossing_time_inverse =
      soundspeed_physical * kernel_support_physical_inv;

  /* Construct time differential of div.v implicitly following the ANARCHY spec
   */

  const float div_v = part_get_div_v(pind);
  const float div_v_prev_step = part_get_div_v_previous_step(pind);

  const float div_v_dt =
      dt_alpha == 0.f
          ? 0.f
          : (div_v - div_v_previous_step) / dt_alpha;

  /* Construct the source term for the AV; if shock detected this is _positive_
   * as div_v_dt should be _negative_ before the shock hits */
  const float S = kernel_support_physical * kernel_support_physical *
                  max(0.f, -1.f * div_v_dt);
  /* 0.25 factor comes from our definition of v_sig (sum of soundspeeds rather
   * than mean). */
  /* Note this is v_sig_physical squared, not comoving */
  const float v_sig_square = 0.25 * v_sig_physical * v_sig_physical;

  /* Calculate the current appropriate value of the AV based on the above */
  const float alpha_loc =
      hydro_props->viscosity.alpha_max * S / (v_sig_square + S);

  float alpha_visc = part_get_alpha_av(pind);

  if (alpha_loc > alpha_visc) {
    /* Reset the value of alpha to the appropriate value */
    alpha_visc = alpha_loc;
  } else {
    /* Integrate the alpha forward in time to decay back to alpha = alpha_loc */
    alpha_visc =
        alpha_loc + (alpha_visc - alpha_loc) *
                        expf(-dt_alpha * sound_crossing_time_inverse *
                             hydro_props->viscosity.length);
  }

  /* Check that we did not hit the minimum */
  alpha_visc =
      max(alpha_visc, hydro_props->viscosity.alpha_min);
  part_set_alpha_av(pind, alpha_visc);

  /* Set our old div_v to the one for the next loop */
  part_set_div_v_previous_step(pind, div_v);
  part_set_div_v_dt(pind, div_v_dt);

  /* Now for the diffusive alpha */
  const float alpha_diff = part_get_alpha_diff(pind);

  const float diffusion_timescale_physical_inverse =
      v_sig_physical * kernel_support_physical_inv;

  const float sqrt_u_inv = 1.f / sqrtf(part_get_u(pind));

  /* Calculate initial value of alpha dt before bounding */
  /* Evolution term: following Schaller+ 2015. This is made up of several
     cosmology factors: physical smoothing length, sound speed from laplace(u) /
     sqrt(u), and the 1 / a^2 coming from the laplace operator. */
  float alpha_diff_dt = hydro_props->diffusion.beta * kernel_support_physical *
    part_get_laplace_u(pind) * cosmo->a_factor_sound_speed *
                        sqrt_u_inv * cosmo->a2_inv;

  /* Decay term: not documented in Schaller+ 2015 but was present
   * in the original EAGLE code and in appendix of Schaye+ 2015 */
  alpha_diff_dt -= (alpha_diff - hydro_props->diffusion.alpha_min) *
                   diffusion_timescale_physical_inverse;

  float new_diffusion_alpha = alpha_diff;
  new_diffusion_alpha += alpha_diff_dt * dt_alpha;

  /* Consistency checks to ensure min < alpha < max */
  new_diffusion_alpha =
      min(new_diffusion_alpha, hydro_props->diffusion.alpha_max);
  new_diffusion_alpha =
      max(new_diffusion_alpha, hydro_props->diffusion.alpha_min);

  part_set_alpha_diff(pind, new_diffusion_alpha);
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking  place in the various force tasks.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_reset_acceleration(
    size_t pind) {

  /* Reset the acceleration. */
  part_set_a_hydro_ind(pind, 0, 0.f);
  part_set_a_hydro_ind(pind, 1, 0.f);
  part_set_a_hydro_ind(pind, 2, 0.f);

  /* Reset the time derivatives. */
  part_set_u_dt(pind, 0.0f);
  part_set_h_dt(pind, 0.0f);
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
    size_t pind,
    const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor) {

  /* Re-set the predicted velocities */
  part_set_v_ind(pind, 0, part_get_v_full_ind(pind,0));
  part_set_v_ind(pind, 1, part_get_v_full_ind(pind,1));
  part_set_v_ind(pind, 2, part_get_v_full_ind(pind,2));
  
  /* Re-set the entropy */
  part_set_u(pind, part_get_u_full(pind));

  /* Compute the sound speed */
  const float soundspeed = hydro_get_comoving_soundspeed(pind);

  part_set_soundspeed(pind, soundspeed);
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * Additional hydrodynamic quantites are drifted forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * Note the different time-step sizes used for the different quantities as they
 * include cosmological factors.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 * @param dt_kick_grav The time-step for gravity quantities.
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    size_t pind, float dt_drift,
    float dt_therm, float dt_kick_grav, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Store ratio of new internal energy to old internal energy, as we use this
   * in the drifting of the pressure. */
  float internal_energy_ratio = 1.f / part_get_u(pind);

  /* Predict the internal energy */
  float u = part_get_u(pind);
  u += part_get_u_dt(pind) * dt_therm;

  float h = part_get_h(pind);
  const float h_inv = 1.f / h;

  /* Predict smoothing length */
  const float w1 = part_get_h_dt(pind) * h_inv * dt_drift;
  if (fabsf(w1) < 0.2f)
    h *= approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    h *= expf(w1);
  
  part_set_h(pind, h);

  /* Predict density and weighted pressure */
  float rho = part_get_rho(pind);
  float pressure_bar = part_get_pressure_bar(pind);
  const float w2 = -hydro_dimension * w1;
  if (fabsf(w2) < 0.2f) {
    const float expf_approx =
        approx_expf(w2); /* 4th order expansion of exp(w) */
    rho *= expf_approx;
    pressure_bar *= expf_approx;
  } else {
    const float expf_exact = expf(w2);
    rho *= expf_exact;
    pressure_bar *= expf_exact;
  }
  part_set_rho(pind, rho);
  part_set_pressure_bar(pind, pressure_bar);

  /* Check against entropy floor */
  const float floor_A = entropy_floor(pind, cosmo, floor_props);
  const float floor_u = gas_internal_energy_from_entropy(part_get_rho(pind), floor_A);

  /* Check against absolute minimum */
  const float min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  u = max(u, floor_u);
  u = max(u, min_u);
  part_set_u(pind, u);

  /* Now that p->u has been properly bounded, we can use it to apply the
   * drift for the pressure */
  internal_energy_ratio *= part_get_u(pind);

  /* Now we can use this to 'update' the value of the smoothed pressure. To
   * truly update this variable, we would need another loop over neighbours
   * using the new internal energies of everyone, but that's not feasible. */
  float pressure_bar = part_get_pressure_bar(pind);
  pressure_bar *= internal_energy_ratio;
  part_set_pressure_bar(pind, pressure_bar);

  /* Compute the new sound speed */
  const float soundspeed =
    gas_soundspeed_from_pressure(part_get_rho(pind), part_get_pressure_bar(pind));
  part_set_soundspeed(pind, soundspeed);

  /* Update the signal velocity */
  part_set_v_sig(pind, max(part_get_v_sig(pind), 2.f * soundspeed));
}

/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the force and accelerations by the appropiate constants
 * and add the self-contribution term. In most cases, there is little
 * to do here.
 *
 * Cosmological terms are also added/multiplied here.
 *
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_force(
    size_t pind, const struct cosmology *cosmo) {

  const float h = part_get_h(pind);
  const float h_dt = part_get_h_dt(pind);
  part_set_h_dt(pind, h_dt * h * hydro_dimension_inv);

}

/**
 * @brief Kick the additional variables
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 * @param dt_grav The time-step for this kick (for gravity quantities).
 * @param dt_grav_mesh The time-step for this kick (mesh gravity).
 * @param dt_hydro The time-step for this kick (for hydro quantities).
 * @param dt_kick_corr The time-step for this kick (for gravity corrections).
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    size_t pind, float dt_therm,
    float dt_grav, float dt_grav_mesh, float dt_hydro, float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Integrate the internal energy forward in time */
  const float delta_u = part_get_u_dt(pind) * dt_therm;

  /* Do not decrease the energy by more than a factor of 2*/
  float u_full = part_get_u_full(pind);
  float u_full_new = max(u_full + delta_u, 0.5f * u_full);
  part_set_u_full(pind, u_full_new);
  
  /* Check against entropy floor */
  const float floor_A = entropy_floor(p, cosmo, floor_props);
  const float floor_u = gas_internal_energy_from_entropy(part_get_rho(pind), floor_A);

  /* Check against absolute minimum */
  const float min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  /* Take highest of both limits */
  const float energy_min = max(min_u, floor_u);

  if (part_get_u_full(pind) < energy_min) {
    part_set_u_full(pind, energy_min);
    part_set_u_dt(pind, 0.f);
  }
}

/**
 * @brief Converts hydro quantity of a particle at the start of a run
 *
 * This function is called once at the end of the engine_init_particle()
 * routine (at the start of a calculation) after the densities of
 * particles have been computed.
 * This can be used to convert internal energy into entropy for instance.
 *
 * @param p The particle to act upon
 * @param xp The extended particle to act upon
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    size_t pind,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Convert the physcial internal energy to the comoving one. */
  /* u' = a^(3(g-1)) u */
  const float factor = 1.f / cosmo->a_factor_internal_energy;
  float u = part_get_u(pind);
  u *= factor;
  part_set_u_full(pind, u);

  /* Apply the minimal energy limit */
  const float min_energy =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;
  if (part_get_u_full(pind) < min_energy) {
    part_set_u_full(pind, min_energy);
    u = min_energy;
    part_set_u_dt(pind, 0.f);
  }
  part_set_u(pind, u);

  /* Note that unlike Minimal the pressure and sound speed cannot be calculated
   * here because they are smoothed properties in this scheme. */

  /* Set the initial value of the artificial viscosity based on the non-variable
     schemes for safety */

  part_set_alpha_av(pind, hydro_props->viscosity.alpha);
  /* Initialise this here to keep all the AV variables together */
  part_set_div_v_previous_step = (pind, 0.f);

  /* Set the initial values for the thermal diffusion */
  part_set_alpha_diff(pind, hydro_props->diffusion.alpha);
}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_first_init_part(
    size_t pind) {
  part_set_time_bin(pind, 0);

  part_set_v_full_ind(pind, 0, part_get_v_ind(pind, 0));
  part_set_v_full_ind(pind, 1, part_get_v_ind(pind, 1));
  part_set_v_full_ind(pind, 2, part_get_v_ind(pind, 2));
  part_set_u_full(pind, part_get_u(pind));


  hydro_reset_acceleration(pind);
  hydro_init_part(pind, NULL);
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
hydro_set_init_internal_energy(size_t pind, float u_init) {

  part_set_u(pind, u_init);
}

#endif /* SWIFT_ANARCHY_PU_HYDRO_H */
