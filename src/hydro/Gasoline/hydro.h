/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_GASOLINE_HYDRO_H
#define SWIFT_GASOLINE_HYDRO_H

/**
 * @file Gasoline/hydro.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added Gasoline physics (Wadsley+ 2017) (Non-neighbour loop
 *        equations)
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
hydro_get_comoving_internal_energy(const struct part *restrict p,
                                   const struct xpart *restrict xp) {
  return xp->u_full;
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
hydro_get_physical_internal_energy(const struct part *restrict p,
                                   const struct xpart *restrict xp,
                                   const struct cosmology *cosmo) {
  return xp->u_full * cosmo->a_factor_internal_energy;
}

/**
 * @brief Returns the comoving internal energy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_internal_energy(const struct part *restrict p) {
  return p->u;
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
  return p->u * cosmo->a_factor_internal_energy;
}

/**
 * @brief Returns the comoving pressure of a particle
 *
 * Computes the pressure based on the particle's properties.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_pressure(
    const struct part *restrict p) {
  return gas_pressure_from_internal_energy(p->rho, p->u);
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
    const struct part *restrict p, const struct cosmology *cosmo) {
  return cosmo->a_factor_pressure * hydro_get_comoving_pressure(p);
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
    const struct part *restrict p, const struct xpart *restrict xp) {
  return gas_entropy_from_internal_energy(p->rho, xp->u_full);
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
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct cosmology *cosmo) {
  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return gas_entropy_from_internal_energy(p->rho, xp->u_full);
}

/**
 * @brief Returns the comoving entropy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_entropy(const struct part *restrict p) {
  return gas_entropy_from_internal_energy(p->rho, p->u);
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
  return gas_entropy_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the comoving sound speed of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_soundspeed(const struct part *restrict p) {
  return gas_soundspeed_from_internal_energy(p->rho, p->u);
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
  return cosmo->a_factor_sound_speed * hydro_get_comoving_soundspeed(p);
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
 * @brief Returns the comoving density of a particle.
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_density(
    const struct part *restrict p, const struct cosmology *cosmo) {
  return cosmo->a3_inv * p->rho;
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
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy_dt(const struct part *restrict p) {
  return p->u_dt;
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
hydro_get_physical_internal_energy_dt(const struct part *restrict p,
                                      const struct cosmology *cosmo) {
  return p->u_dt * cosmo->a_factor_internal_energy;
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
hydro_set_comoving_internal_energy_dt(struct part *restrict p, float du_dt) {
  p->u_dt = du_dt;
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
hydro_set_physical_internal_energy_dt(struct part *restrict p,
                                      const struct cosmology *cosmo,
                                      float du_dt) {
  p->u_dt = du_dt / cosmo->a_factor_internal_energy;
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
  const float comoving_entropy = entropy;
  xp->u_full = gas_internal_energy_from_entropy(p->rho, comoving_entropy);
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
hydro_set_physical_internal_energy(struct part *p, struct xpart *xp,
                                   const struct cosmology *cosmo,
                                   const float u) {
  xp->u_full = u / cosmo->a_factor_internal_energy;
}

/**
 * @brief Sets the drifted physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param pressure_floor The properties of the pressure floor.
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(
    struct part *p, const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor, const float u) {
  /* There is no need to use the floor here as this function is called in the
   * feedback, so the new value of the internal energy should be strictly
   * higher than the old value. */

  p->u = u / cosmo->a_factor_internal_energy;

  /* Now recompute the extra quantities */

  /* Compute the sound speed */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  /* Update variables. */
  p->force.soundspeed = soundspeed;
  p->force.pressure = pressure_including_floor;

  p->viscosity.v_sig = max(p->viscosity.v_sig, soundspeed);
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
hydro_set_v_sig_based_on_velocity_kick(struct part *p,
                                       const struct cosmology *cosmo,
                                       const float dv_phys) {
  /* Compute the velocity kick in comoving coordinates */
  const float dv = dv_phys / cosmo->a_factor_sound_speed;

  /* Sound speed */
  const float soundspeed = hydro_get_comoving_soundspeed(p);

  /* Update the signal velocity */
  p->viscosity.v_sig =
      max(soundspeed, p->viscosity.v_sig + const_viscosity_beta * dv);
}

/**
 * @brief Update the value of the viscosity alpha for the scheme.
 *
 * @param p the particle of interest
 * @param alpha the new value for the viscosity coefficient.
 */
__attribute__((always_inline)) INLINE static void hydro_set_viscosity_alpha(
    struct part *restrict p, float alpha) {
  p->viscosity.alpha = alpha;
}

/**
 * @brief Update the value of the diffusive coefficients to the
 *        feedback reset value for the scheme.
 *
 * @param p the particle of interest
 */
__attribute__((always_inline)) INLINE static void
hydro_diffusive_feedback_reset(struct part *restrict p) {
  /* Set the viscosity to the max, and the diffusion to the min */
  hydro_set_viscosity_alpha(p,
                            hydro_props_default_viscosity_alpha_feedback_reset);

  p->diffusion.rate = hydro_props_default_diffusion_coefficient_feedback_reset;
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
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct hydro_props *restrict hydro_properties,
    const struct cosmology *restrict cosmo) {
  const float CFL_condition = hydro_properties->CFL_condition;

  /* CFL condition */
  const float dt_cfl = 2.f * kernel_gamma * CFL_condition * cosmo->a * p->h /
                       (cosmo->a_factor_sound_speed * p->viscosity.v_sig);

  return dt_cfl;
}

/**
 * @brief Compute the signal velocity between two gas particles
 *
 * @param dx Comoving vector separating both particles (pi - pj).
 * @brief pi The first #part.
 * @brief pj The second #part.
 * @brief mu_ij The velocity on the axis linking the particles, or zero if the
 * particles are moving away from each other,
 * @brief beta The non-linear viscosity constant.
 */
__attribute__((always_inline)) INLINE static float hydro_signal_velocity(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const float mu_ij, const float beta) {

  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  return 0.5f * (ci + cj) - mu_ij;
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
 * @brief Operations performed when a particle gets removed from the
 * simulation volume.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @param time The simulation time.
 */
__attribute__((always_inline)) INLINE static void hydro_remove_part(
    const struct part *p, const struct xpart *xp, const double time) {}

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
    struct part *restrict p, const struct hydro_space *hs) {
  p->density.wcount = 0.f;
  p->density.wcount_dh = 0.f;
  p->rho = 0.f;
  p->weighted_wcount = 0.f;
  p->weighted_neighbour_wcount = 0.f;
  p->density.rho_dh = 0.f;

  p->smooth_pressure_gradient[0] = 0.f;
  p->smooth_pressure_gradient[1] = 0.f;
  p->smooth_pressure_gradient[2] = 0.f;

  p->viscosity.shock_indicator = 0.f;
  p->viscosity.shock_limiter = 0.f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      p->viscosity.velocity_gradient[i][j] = 0.f;
    }
  }
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
    struct part *restrict p, const struct cosmology *cosmo) {
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

  /* Finish caclulation of the pressure gradients; note that the
   * kernel derivative is zero at our particle's position. */
  p->smooth_pressure_gradient[0] *= hydro_gamma_minus_one * h_inv_dim_plus_one;
  p->smooth_pressure_gradient[1] *= hydro_gamma_minus_one * h_inv_dim_plus_one;
  p->smooth_pressure_gradient[2] *= hydro_gamma_minus_one * h_inv_dim_plus_one;

  /* Finish calculation of the velocity gradient tensor, and
   * guard against FPEs here. */
  const float velocity_gradient_norm =
      p->weighted_wcount == 0.f ? 0.f
                                : 3.f * cosmo->a2_inv / p->weighted_wcount;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      const float hubble_term = i == j ? hydro_dimension * cosmo->H : 0.f;

      p->viscosity.velocity_gradient[i][j] *= velocity_gradient_norm;

      p->viscosity.velocity_gradient[i][j] += hubble_term;
    }
  }
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
 * @param pressure_floor The properties of the pressure floor.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_gradient(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Compute the sound speed  */
  const float pressure = hydro_get_comoving_pressure(p);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  /* Update variables. */
  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;

  const float mod_pressure_gradient =
      sqrtf(p->smooth_pressure_gradient[0] * p->smooth_pressure_gradient[0] +
            p->smooth_pressure_gradient[1] * p->smooth_pressure_gradient[1] +
            p->smooth_pressure_gradient[2] * p->smooth_pressure_gradient[2]);

  float unit_pressure_gradient[3];

  /* As this is normalised, no cosmology factor is required */
  unit_pressure_gradient[0] =
      p->smooth_pressure_gradient[0] / mod_pressure_gradient;
  unit_pressure_gradient[1] =
      p->smooth_pressure_gradient[1] / mod_pressure_gradient;
  unit_pressure_gradient[2] =
      p->smooth_pressure_gradient[2] / mod_pressure_gradient;

  float dv_dn = 0.f;
  float shear_norm2 = 0.f;
  float traceless_shear_norm2 = 0.f;
  float div_v = 0.f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dv_dn += unit_pressure_gradient[i] *
               p->viscosity.velocity_gradient[i][j] * unit_pressure_gradient[j];

      const float shear_component =
          0.5f * (p->viscosity.velocity_gradient[i][j] +
                  p->viscosity.velocity_gradient[j][i]);
      const float shear_component2 = shear_component * shear_component;

      shear_norm2 += shear_component2;
      traceless_shear_norm2 += i == j ? 0.f : shear_component2;
      div_v += i == j ? (1. / 3.) * p->viscosity.velocity_gradient[i][j] : 0.f;
    }
  }

  const float shock_indicator =
      (3.f / 2.f) * (dv_dn + max(-(1.f / 3.f) * div_v, 0.f));

  /* Now do the conduction coefficient; note that no limiter is used here. */
  /* These square roots are not included in the original documentation */
  const float h_physical = p->h * cosmo->a;

  const float diffusion_rate = hydro_props->diffusion.coefficient *
                               sqrtf(traceless_shear_norm2) * h_physical *
                               h_physical;

  p->diffusion.rate = diffusion_rate;
  p->viscosity.tensor_norm = sqrtf(shear_norm2);
  p->viscosity.shock_indicator = shock_indicator;
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
    struct part *restrict p) {
  p->viscosity.v_sig = p->force.soundspeed;
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
    struct part *p) {
  /* The f_i is calculated explicitly in Gasoline. */
  p->force.f = p->weighted_wcount / (p->weighted_neighbour_wcount * p->rho);

  /* Calculate smoothing length powers */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Apply correct normalisation */
  p->viscosity.shock_limiter *= h_inv_dim;
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
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo) {
  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  warning(
      "Gas particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      p->id, h, p->density.wcount);

  /* Re-set problematic values */
  p->rho = p->mass * kernel_root * h_inv_dim;
  p->viscosity.v_sig = 0.f;
  p->density.wcount = kernel_root * h_inv_dim;
  p->density.rho_dh = 0.f;
  p->density.wcount_dh = 0.f;
  /* Set to 1 as these are only used by taking the ratio */
  p->weighted_wcount = 1.f;
  p->weighted_neighbour_wcount = 1.f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      p->viscosity.velocity_gradient[i][j] = 0.f;
    }
  }
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
 * @param pressure_floor The properties of the pressure floor.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor, const float dt_alpha,
    const float dt_therm) {
  /* Here we need to update the artificial viscosity alpha */
  const float d_shock_indicator_dt =
      dt_alpha == 0.f ? 0.f
                      : (p->viscosity.shock_indicator -
                         p->viscosity.shock_indicator_previous_step) /
                            dt_alpha;

  /* A_i in all of the equations */
  const float v_sig_physical = p->viscosity.v_sig * cosmo->a_factor_sound_speed;
  const float soundspeed_physical =
      gas_soundspeed_from_pressure(hydro_get_physical_density(p, cosmo),
                                   hydro_get_physical_pressure(p, cosmo));
  const float h_physical = p->h * cosmo->a;

  /* Note that viscosity.shock_limiter still includes the cosmology dependence
   * from the density, so what we do here is correct (i.e. include no explicit
   * h-factors). */
  const float shock_limiter_core =
      0.5f * (1.f - p->viscosity.shock_limiter / p->rho);
  const float shock_limiter_core_2 = shock_limiter_core * shock_limiter_core;
  const float shock_limiter = shock_limiter_core_2 * shock_limiter_core_2;

  const float shock_detector = 2.f * h_physical * h_physical * kernel_gamma2 *
                               shock_limiter *
                               max(-1.f * d_shock_indicator_dt, 0.f);

  const float alpha_loc =
      hydro_props->viscosity.alpha_max *
      (shock_detector / (shock_detector + v_sig_physical * v_sig_physical));

  /* TODO: Probably use physical _h_ throughout this function */
  const float d_alpha_dt = (alpha_loc - p->viscosity.alpha) *
                           hydro_props->viscosity.length * soundspeed_physical /
                           h_physical;

  float new_alpha = p->viscosity.alpha;

  if (new_alpha < alpha_loc) {
    new_alpha = alpha_loc;
  } else {
    /* Very basic time integration here */
    new_alpha += d_alpha_dt * dt_alpha;
  }

  new_alpha = max(new_alpha, hydro_props->viscosity.alpha_min);
  new_alpha = min(new_alpha, hydro_props->viscosity.alpha_max);

  p->viscosity.shock_indicator_previous_step = p->viscosity.shock_indicator;
  p->viscosity.alpha = new_alpha;
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
    struct part *restrict p) {
  /* Reset the acceleration. */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* Reset the time derivatives. */
  p->u_dt = 0.0f;
  p->force.h_dt = 0.0f;
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 * @param cosmo The cosmological model.
 * @param pressure_floor The properties of the pressure floor.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values(
    struct part *restrict p, const struct xpart *restrict xp,
    const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor) {
  /* Re-set the predicted velocities */
  p->v[0] = xp->v_full[0];
  p->v[1] = xp->v_full[1];
  p->v[2] = xp->v_full[2];

  /* Re-set the entropy */
  p->u = xp->u_full;

  /* Compute the sound speed */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;

  /* Update the signal velocity, if we need to. */
  p->viscosity.v_sig = max(p->viscosity.v_sig, soundspeed);
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
 * @param pressure_floor The properties of the pressure floor.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part *restrict p, const struct xpart *restrict xp, float dt_drift,
    float dt_therm, float dt_kick_grav, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Predict the internal energy */
  p->u += p->u_dt * dt_therm;

  const float h_inv = 1.f / p->h;

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt_drift;
  if (fabsf(w1) < 0.2f)
    p->h *= approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    p->h *= expf(w1);

  /* Predict density and weighted pressure */
  const float w2 = -hydro_dimension * w1;
  if (fabsf(w2) < 0.2f) {
    const float expf_approx =
        approx_expf(w2); /* 4th order expansion of exp(w) */
    p->rho *= expf_approx;
  } else {
    const float expf_exact = expf(w2);
    p->rho *= expf_exact;
  }

  /* Check against entropy floor - explicitly do this after drifting the
   * density as this has a density dependence. */
  const float floor_A = entropy_floor(p, cosmo, floor_props);
  const float floor_u = gas_internal_energy_from_entropy(p->rho, floor_A);

  /* Check against absolute minimum */
  const float min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  p->u = max(p->u, floor_u);
  p->u = max(p->u, min_u);

  /* Compute the new sound speed */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;

  /* Update signal velocity if we need to */
  p->viscosity.v_sig = max(p->viscosity.v_sig, soundspeed);
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
    struct part *restrict p, const struct cosmology *cosmo) {
  p->force.h_dt *= p->h * hydro_dimension_inv;
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
    struct part *restrict p, struct xpart *restrict xp, float dt_therm,
    float dt_grav, float dt_grav_mesh, float dt_hydro, float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {
  /* Integrate the internal energy forward in time */
  const float delta_u = p->u_dt * dt_therm;

  /* Do not decrease the energy by more than a factor of 2*/
  xp->u_full = max(xp->u_full + delta_u, 0.5f * xp->u_full);

  /* Check against entropy floor */
  const float floor_A = entropy_floor(p, cosmo, floor_props);
  const float floor_u = gas_internal_energy_from_entropy(p->rho, floor_A);

  /* Check against absolute minimum */
  const float min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  /* Take highest of both limits */
  const float energy_min = max(min_u, floor_u);

  if (xp->u_full < energy_min) {
    xp->u_full = energy_min;
    p->u_dt = 0.f;
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
 * @param pressure_floor The properties of the pressure floor.
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {
  /* Convert the physcial internal energy to the comoving one. */
  /* u' = a^(3(g-1)) u */
  const float factor = 1.f / cosmo->a_factor_internal_energy;
  p->u *= factor;
  xp->u_full = p->u;

  /* Apply the minimal energy limit */
  const float min_comoving_energy =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;
  if (xp->u_full < min_comoving_energy) {
    xp->u_full = min_comoving_energy;
    p->u = min_comoving_energy;
    p->u_dt = 0.f;
  }

  /* Set the initial value of the artificial viscosity based on the non-variable
     schemes for safety */

  p->viscosity.alpha = hydro_props->viscosity.alpha;

  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;
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
    struct part *restrict p, struct xpart *restrict xp) {
  p->time_bin = 0;

  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
  xp->u_full = p->u;

  p->viscosity.shock_indicator_previous_step = 0.f;
  p->viscosity.tensor_norm = 0.f;

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
  p->u = u_init;
}

#endif /* SWIFT_GASOLINE_HYDRO_H */
