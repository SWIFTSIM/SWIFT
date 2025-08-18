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
#ifndef SWIFT_SPHENIX_HYDRO_H
#define SWIFT_SPHENIX_HYDRO_H

/**
 * @file SPHENIX/hydro.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added SPHENIX physics (Borrow 2020)  (Non-neighbour loop
 *        equations)
 */

#include "adiabatic_index.h"
#include "approx_math.h"
#include "cosmology.h"
#include "dimension.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "fvpm_geometry.h"
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
__attribute__((always_inline)) INLINE static double
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
__attribute__((always_inline)) INLINE static double
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
__attribute__((always_inline)) INLINE static double
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
__attribute__((always_inline)) INLINE static double
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
__attribute__((always_inline)) INLINE static double hydro_get_comoving_pressure(
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
__attribute__((always_inline)) INLINE static double hydro_get_physical_pressure(
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
__attribute__((always_inline)) INLINE static double hydro_get_comoving_entropy(
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
__attribute__((always_inline)) INLINE static double hydro_get_physical_entropy(
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
__attribute__((always_inline)) INLINE static double
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
__attribute__((always_inline)) INLINE static double
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
__attribute__((always_inline)) INLINE static double
hydro_get_comoving_soundspeed(const struct part *restrict p) {

  return gas_soundspeed_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the physical sound speed of a particle
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static double
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
__attribute__((always_inline)) INLINE static double hydro_get_mass(
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
    struct part *restrict p, double m) {

  p->mass = m;
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static double
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
__attribute__((always_inline)) INLINE static double
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
hydro_set_comoving_internal_energy_dt(struct part *restrict p, double du_dt) {

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
                                      double du_dt) {

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
    const double entropy) {

  /* Note there is no conversion from physical to comoving entropy */
  const double comoving_entropy = entropy;
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
                                   const double u) {

  xp->u_full = u / cosmo->a_factor_internal_energy;
}

__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy_TESTING_SPH_RT(struct part *p,
                                                  const struct cosmology *cosmo,
                                                  const double u) {

  // TODO: This might be a problem. Hacky version to get code to compile.
  // TODO: Cosmology might need attention
  p->u = u / cosmo->a_factor_internal_energy;
}

/**
 * @brief Sets the drifted physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param pressure_floor The #pressure_floor_props used.
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(
    struct part *p, const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor, const double u) {

  /* There is no need to use the floor here as this function is called in the
   * feedback, so the new value of the internal energy should be strictly
   * higher than the old value. */

  p->u = u / cosmo->a_factor_internal_energy;

  /* Now recompute the extra quantities */

  /* Compute the sound speed */
  const double pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const double pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const double soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  /* Update variables. */
  p->force.soundspeed = soundspeed;
  p->force.pressure = pressure_including_floor;

  p->viscosity.v_sig = fmax(p->viscosity.v_sig, 2.0 * soundspeed);
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
                                       const double dv_phys) {

  /* Compute the velocity kick in comoving coordinates */
  const double dv = dv_phys / cosmo->a_factor_sound_speed;

  /* Sound speed */
  const double soundspeed = hydro_get_comoving_soundspeed(p);

  /* Update the signal velocity */
  p->viscosity.v_sig =
      fmax(2.0 * soundspeed, p->viscosity.v_sig + const_viscosity_beta * dv);
}

/**
 * @brief Update the value of the viscosity alpha for the scheme.
 *
 * @param p the particle of interest
 * @param alpha the new value for the viscosity coefficient.
 */
__attribute__((always_inline)) INLINE static void hydro_set_viscosity_alpha(
    struct part *restrict p, double alpha) {
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

  p->diffusion.alpha = hydro_props_default_diffusion_alpha_feedback_reset;
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
__attribute__((always_inline)) INLINE static double hydro_compute_timestep(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct hydro_props *restrict hydro_properties,
    const struct cosmology *restrict cosmo) {

  const double CFL_condition = hydro_properties->CFL_condition;

  /* CFL condition */
  const double dt_cfl = 2.0 * kernel_gamma * CFL_condition * cosmo->a * p->h /
                       (cosmo->a_factor_sound_speed * p->viscosity.v_sig);

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
__attribute__((always_inline)) INLINE static double hydro_signal_velocity(
    const double dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const double mu_ij, const double beta) {

  const double ci = pi->force.soundspeed;
  const double cj = pj->force.soundspeed;

  return ci + cj - beta * mu_ij;
}

/**
 * @brief returns the signal velocity
 *
 * @brief p  the particle
 */
__attribute__((always_inline)) INLINE static double hydro_get_signal_velocity(
    const struct part *restrict p) {

  return p->viscosity.v_sig;
}

/**
 * @brief returns the div_v
 *
 * @brief p  the particle
 */
__attribute__((always_inline)) INLINE static double hydro_get_div_v(
    const struct part *restrict p) {

  return p->viscosity.div_v;
}

/**
 * @brief Does some extra hydro operations once the actual physical time step
 * for the particle is known.
 *
 * @param p The particle to act upon.
 * @param dt Physical time step of the particle during the next step.
 */
__attribute__((always_inline)) INLINE static void hydro_timestep_extra(
    struct part *p, double dt) {}

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

  p->density.wcount = 0.0;
  p->density.wcount_dh = 0.0;
  p->rho = 0.0;
  p->density.rho_dh = 0.0;

  p->diffusion.laplace_u = 0.0;

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  p->N_density = 1; /* Self contribution */
  p->N_force = 0;
  p->N_gradient = 1;
  p->N_density_exact = 0;
  p->N_force_exact = 0;
  p->rho_exact = 0.f;
  p->n_gradient = 0.f;
  p->n_gradient_exact = 0.f;
  p->n_density = 0.f;
  p->n_density_exact = 0.f;
  p->n_force = 0.f;
  p->n_force_exact = 0.f;
  p->inhibited_exact = 0;
  p->limited_part = 0;
#endif

  /* Init geometry for FVPM Radiative Transfer */
  fvpm_geometry_init(p);
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
  const double h = p->h;
  const double h_inv = 1.0 / h;                       /* 1/h */
  const double h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const double h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

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

  /* Finish matrix and volume computations for FVPM Radiative Transfer */
  fvpm_compute_volume_and_matrix(p, h_inv_dim);

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  p->n_density += kernel_root;
  p->n_density *= h_inv_dim;
#endif
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
 * @param pressure_floor The #pressure_floor_props used.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_gradient(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Compute the sound speed  */
  const double pressure = hydro_get_comoving_pressure(p);
  const double pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const double soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  /* Compute the "grad h" term  - Note here that we have \tilde{x}
   * as 1 as we use the local number density to find neighbours. This
   * introduces a j-component that is considered in the force loop,
   * meaning that this cached grad_h_term gives:
   *
   * f_ij = 1.f - grad_h_term_i / m_j */
  const double common_factor = p->h * hydro_dimension_inv / p->density.wcount;
  double grad_h_term;

  /* Ignore changing-kernel effects when h ~= h_max */
  if (p->h > 0.9999 * hydro_props->h_max) {
    grad_h_term = 0.0;
    warning("h ~ h_max for particle with ID %lld (h: %g)", p->id, p->h);
  } else {
    const double grad_W_term = common_factor * p->density.wcount_dh;
    if (grad_W_term < -0.9999) {
      /* if we get here, we either had very small neighbour contributions
         (which should be treated as a no neighbour case in the ghost) or
         a very weird particle distribution (e.g. particles sitting on
         top of each other). Either way, we cannot use the normal
         expression, since that would lead to overflow or excessive round
         off and cause excessively high accelerations in the force loop */
      grad_h_term = 0.0;
      warning(
          "grad_W_term very small for particle with ID %lld (h: %g, wcount: "
          "%g, wcount_dh: %g)",
          p->id, p->h, p->density.wcount, p->density.wcount_dh);
    } else {
      grad_h_term = common_factor * p->density.rho_dh / (1.0 + grad_W_term);
    }
  }

  /* Update variables. */
  p->force.f = grad_h_term;
  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;
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

  p->viscosity.rot_v[0] = 0.0;
  p->viscosity.rot_v[1] = 0.0;
  p->viscosity.rot_v[2] = 0.0;

  p->viscosity.div_v = 0.0;

  p->viscosity.v_sig = 2.0 * p->force.soundspeed;
  p->force.alpha_visc_max_ngb = p->viscosity.alpha;
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
    struct part *p, const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor) {

  /* Some smoothing length multiples. */
  const double h = p->h;
  const double h_inv = 1.0 / h;                       /* 1/h */
  const double h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const double h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  const double rho_inv = 1.0 / p->rho;
  const double a_inv2 = cosmo->a2_inv;

  /* Finish calculation of the velocity curl components */
  p->viscosity.rot_v[0] *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  p->viscosity.rot_v[1] *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  p->viscosity.rot_v[2] *= h_inv_dim_plus_one * a_inv2 * rho_inv;

  /* Finish calculation of the velocity divergence */
  p->viscosity.div_v *= h_inv_dim_plus_one * rho_inv * a_inv2;
  p->viscosity.div_v += cosmo->H * hydro_dimension;

  /* Include the extra factors in the del^2 u */

  p->diffusion.laplace_u *= 2.0 * h_inv_dim_plus_one;

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  p->n_gradient += kernel_root;
#endif
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
  const double h = p->h;
  const double h_inv = 1.0 / h;                 /* 1/h */
  const double h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  warning(
      "Gas particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      p->id, h, p->density.wcount);

  /* Re-set problematic values */
  p->rho = p->mass * kernel_root * h_inv_dim;
  p->viscosity.v_sig = 0.0;
  p->density.wcount = kernel_root * h_inv_dim;
  p->density.rho_dh = 0.0;
  p->density.wcount_dh = 0.0;

  p->viscosity.rot_v[0] = 0.0;
  p->viscosity.rot_v[1] = 0.0;
  p->viscosity.rot_v[2] = 0.0;

  /* Probably not shocking, so this is safe to do */
  p->viscosity.div_v = 0.0;
  p->diffusion.laplace_u = 0.0;
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
 * @param pressure_floor The #pressure_floor_props used.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor, const double dt_alpha,
    const double dt_therm) {

  /* Here we need to update the artificial viscosity */

  /* We use in this function that h is the radius of support */
  const double kernel_support_physical = p->h * cosmo->a * kernel_gamma;
  const double kernel_support_physical_inv = 1.0 / kernel_support_physical;
  const double v_sig_physical = p->viscosity.v_sig * cosmo->a_factor_sound_speed;
  const double pressure = hydro_get_comoving_pressure(p);
  const double pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const double soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);
  const double soundspeed_physical = soundspeed * cosmo->a_factor_sound_speed;
  const double fac_B = cosmo->a_factor_Balsara_eps;

  const double sound_crossing_time_inverse =
      soundspeed_physical * kernel_support_physical_inv;

  /* Compute the norm of the curl */
  const double curl_v = sqrt(p->viscosity.rot_v[0] * p->viscosity.rot_v[0] +
                             p->viscosity.rot_v[1] * p->viscosity.rot_v[1] +
                             p->viscosity.rot_v[2] * p->viscosity.rot_v[2]);

  /* Compute the norm of div v */
  const double abs_div_v = fabs(p->viscosity.div_v);

  /* Get the squares of the quantities necessary for the Balsara-like switch */
  const double fac_B_2 = fac_B * fac_B;
  const double curl_v_2 = curl_v * curl_v;
  const double abs_div_v_2 = abs_div_v * abs_div_v;
  const double soundspeed_2 = soundspeed * soundspeed;
  const double h_2 = p->h * p->h;

  /* Compute the Balsara-like switch (Price et a. (2018), eq. 47; a simplified
   * version of eq. 18 in Cullen & Dehnen (2012)) */
  const double balsara = abs_div_v_2 / (abs_div_v_2 + curl_v_2 +
                                       0.0001f * soundspeed_2 * fac_B_2 / h_2);

  /* Construct time differential of div.v implicitly following the ANARCHY spec
   */
  const double div_v_dt =
      dt_alpha == 0.0
          ? 0.0
          : (p->viscosity.div_v - p->viscosity.div_v_previous_step) / dt_alpha;

  /* Construct the source term for the AV; if shock detected this is _positive_
   * as div_v_dt should be _negative_ before the shock hits */
  /* Source term is only activated if flow is converging (i.e. in the pre-
   * shock region) */
  const double S = p->viscosity.div_v < 0.0
                      ? balsara * kernel_support_physical *
                            kernel_support_physical * fmax(0.0, -1.0 * div_v_dt)
                      : 0.0;

  /* We want the decay to occur based on the thermodynamic properties
   * of the gas - hence use the soundspeed instead of the signal velocity */
  const double soundspeed_square = soundspeed_physical * soundspeed_physical;

  /* Calculate the current appropriate value of the AV based on the above */
  const double alpha_loc =
      hydro_props->viscosity.alpha_max * S / (soundspeed_square + S);

  if (alpha_loc > p->viscosity.alpha) {
    /* Reset the value of alpha to the appropriate value */
    p->viscosity.alpha = alpha_loc;
  } else {
    /* Integrate the alpha forward in time to decay back to alpha = alpha_loc */
    /* This type of integration is stable w.r.t. different time-step lengths
     * (Price 2018) */
    const double timescale_ratio =
        dt_alpha * sound_crossing_time_inverse * hydro_props->viscosity.length;

    p->viscosity.alpha += alpha_loc * timescale_ratio;
    p->viscosity.alpha /= (1.0 + timescale_ratio);
  }

  /* Check that we did not hit the minimum */
  p->viscosity.alpha =
      fmax(p->viscosity.alpha, hydro_props->viscosity.alpha_min);

  /* Set our old div_v to the one for the next loop */
  p->viscosity.div_v_previous_step = p->viscosity.div_v;
  p->viscosity.div_v_dt = div_v_dt;

  /* Now for the diffusive alpha */

  const double diffusion_timescale_physical_inverse =
      v_sig_physical * kernel_support_physical_inv;

  const double sqrt_u_inv = 1.0 / sqrt(p->u);

  /* Calculate initial value of alpha dt before bounding */
  /* Evolution term: following Schaller+ 2015. This is made up of several
     cosmology factors: physical smoothing length, sound speed from laplace(u) /
     sqrt(u), and the 1 / a^2 coming from the laplace operator. */
  double alpha_diff_dt = hydro_props->diffusion.beta * kernel_support_physical *
                        p->diffusion.laplace_u * cosmo->a_factor_sound_speed *
                        sqrt_u_inv * cosmo->a2_inv;

  /* Decay term: not documented in Schaller+ 2015 but was present
   * in the original EAGLE code and in appendix of Schaye+ 2015 */
  alpha_diff_dt -= (p->diffusion.alpha - hydro_props->diffusion.alpha_min) *
                   diffusion_timescale_physical_inverse;

  double new_diffusion_alpha = p->diffusion.alpha;
  new_diffusion_alpha += alpha_diff_dt * dt_alpha;

  /* Consistency checks to ensure min < alpha < max */
  new_diffusion_alpha =
      fmax(new_diffusion_alpha, hydro_props->diffusion.alpha_min);

  /* Now we limit in viscous flows; remove diffusion there. If we
   * don't do that, then we end up diffusing energy away in supernovae.
   * This is an EAGLE-specific fix. We limit based on the maximal
   * viscous alpha over our neighbours in an attempt to keep diffusion
   * low near to supernovae sites. */

  /* This also enforces alpha_diff < alpha_diff_max */

  const double viscous_diffusion_limit =
      hydro_props->diffusion.alpha_max *
      (1.0 - p->force.alpha_visc_max_ngb / hydro_props->viscosity.alpha_max);

  new_diffusion_alpha = fmin(new_diffusion_alpha, viscous_diffusion_limit);

  p->diffusion.alpha = new_diffusion_alpha;
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
  p->a_hydro[0] = 0.0;
  p->a_hydro[1] = 0.0;
  p->a_hydro[2] = 0.0;

  /* Reset the time derivatives. */
  p->u_dt = 0.0;
  p->force.h_dt = 0.0;
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 * @param cosmo The cosmological model.
 * @param pressure_floor The #pressure_floor_props used.
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
  const double pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const double pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const double soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;

  /* Update the signal velocity, if we need to. */
  p->viscosity.v_sig = fmax(p->viscosity.v_sig, 2.0 * soundspeed);
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
    struct part *restrict p, const struct xpart *restrict xp, double dt_drift,
    double dt_therm, double dt_kick_grav, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Predict the internal energy */
  p->u += p->u_dt * dt_therm;

  const double h_inv = 1.0 / p->h;

  /* Predict smoothing length */
  const double w1 = p->force.h_dt * h_inv * dt_drift;
  if (fabs(w1) < 0.2f)
    p->h *= exp(w1); /* 4th order expansion of exp(w) */
  else
    p->h *= exp(w1);

  /* Predict density and weighted pressure */
  const double w2 = -hydro_dimension * w1;
  if (fabs(w2) < 0.2f) {
    const double expf_approx =
        exp(w2); /* 4th order expansion of exp(w) */
    p->rho *= expf_approx;
  } else {
    const double expf_exact = exp(w2);
    p->rho *= expf_exact;
  }

  /* Check against entropy floor - explicitly do this after drifting the
   * density as this has a density dependence. */
  const double floor_A = entropy_floor(p, cosmo, floor_props);
  const double floor_u = gas_internal_energy_from_entropy(p->rho, floor_A);

  /* Check against absolute minimum */
  const double min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  p->u = fmax(p->u, floor_u);
  p->u = fmax(p->u, min_u);

  /* Compute the new sound speed */
  const double pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const double pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const double soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;

  /* Update signal velocity if we need to */
  p->viscosity.v_sig = fmax(p->viscosity.v_sig, 2.0* soundspeed);
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
    struct part *restrict p, struct xpart *restrict xp, double dt_therm,
    double dt_grav, double dt_grav_mesh, double dt_hydro, double dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Integrate the internal energy forward in time */
  const double delta_u = p->u_dt * dt_therm;

  /* Do not decrease the energy by more than a factor of 2*/
  xp->u_full = fmax(xp->u_full + delta_u, 0.5 * xp->u_full);

  /* Check against entropy floor */
  const double floor_A = entropy_floor(p, cosmo, floor_props);
  const double floor_u = gas_internal_energy_from_entropy(p->rho, floor_A);

  /* Check against absolute minimum */
  const double min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  /* Take highest of both limits */
  const double energy_min = fmax(min_u, floor_u);

  if (xp->u_full < energy_min) {
    xp->u_full = energy_min;
    p->u_dt = 0.0;
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
  const double factor = 1.0 / cosmo->a_factor_internal_energy;
  p->u *= factor;
  xp->u_full = p->u;

  /* Apply the minimal energy limit */
  const double min_comoving_energy =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;
  if (xp->u_full < min_comoving_energy) {
    xp->u_full = min_comoving_energy;
    p->u = min_comoving_energy;
    p->u_dt = 0.0;
  }

  /* Set the initial value of the artificial viscosity based on the non-variable
     schemes for safety */

  p->viscosity.alpha = hydro_props->viscosity.alpha;
  /* Initialise this here to keep all the AV variables together */
  p->viscosity.div_v_previous_step = 0.0;
  p->viscosity.div_v_dt = 0.0;

  /* Set the initial values for the thermal diffusion */
  p->diffusion.alpha = hydro_props->diffusion.alpha;

  const double pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const double pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const double soundspeed =
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
hydro_set_init_internal_energy(struct part *p, double u_init) {

  p->u = u_init;
}

#endif /* SWIFT_SPHENIX_HYDRO_H */
