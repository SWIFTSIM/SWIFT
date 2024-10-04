/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_HYDRO_H
#define SWIFT_PLANETARY_HYDRO_H

/**
 * @file Planetary/hydro.h
 * @brief REMIX implementation of SPH (Sandnes et al. 2024)
 */

#include "adiabatic_index.h"
#include "approx_math.h"
#include "cosmology.h"
#include "debug.h"
#include "dimension.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "hydro_kernels.h"
#include "hydro_parameters.h"
#include "hydro_properties.h"
#include "hydro_space.h"
#include "hydro_visc_difn.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "pressure_floor.h"
#include "timestep_sync_part.h"



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

  return gas_pressure_from_internal_energy(p->rho_evol, p->u, p->mat_id);
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

  return cosmo->a_factor_pressure *
         gas_pressure_from_internal_energy(p->rho_evol, p->u, p->mat_id);
}

/**
 * @brief Returns the comoving entropy of a particle
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

  return gas_entropy_from_internal_energy(p->rho_evol, xp->u_full, p->mat_id);
}

/**
 * @brief Returns the physical entropy of a particle
 *
 * For implementations where the main thermodynamic variable
 * is not entropy, this function computes the entropy from
 * the thermodynamic variable and converts it to
 * physical coordinates.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_entropy(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct cosmology *cosmo) {

  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return gas_entropy_from_internal_energy(p->rho_evol, xp->u_full, p->mat_id);
}

/**
 * @brief Returns the comoving entropy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_entropy(const struct part *restrict p) {

  return gas_entropy_from_internal_energy(p->rho_evol, p->u, p->mat_id);
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
  return gas_entropy_from_internal_energy(p->rho_evol, p->u, p->mat_id);
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
 * @brief Returns the time derivative of internal energy of a particle
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
  xp->u_full = gas_internal_energy_from_entropy(p->rho_evol, comoving_entropy,
                                                p->mat_id);
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
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(struct part *p,
                                           const struct cosmology *cosmo,
                                           const float u) {

  p->u = u / cosmo->a_factor_internal_energy;

  /* Now recompute the extra quantities */

  /* Compute the sound speed */
  const float pressure =
      gas_pressure_from_internal_energy(p->rho_evol, p->u, p->mat_id);
  const float soundspeed = hydro_get_comoving_soundspeed(p);

  /* Update variables. */
  p->force.pressure = pressure;
  p->force.soundspeed = soundspeed;

  p->force.v_sig = max(p->force.v_sig, 2.f * soundspeed);
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

__attribute__((always_inline)) INLINE static void
hydro_set_v_sig_based_on_velocity_kick(struct part *p,
                                       const struct cosmology *cosmo,
                                       const float dv_phys) {

  /* Compute the velocity kick in comoving coordinates */
  const float dv = dv_phys / cosmo->a_factor_sound_speed;

  /* Sound speed */
  const float soundspeed = hydro_get_comoving_soundspeed(p);

  // const float viscosity_parameter_factor = (alpha == 0.f) ? 0.f : beta / alpha;
  /* Update the signal velocity */
  p->force.v_sig =
      max(2.f * soundspeed, p->force.v_sig + 2 * dv);
}

__attribute__((always_inline)) INLINE static void hydro_set_velocity1(
    struct part *restrict p, const struct cosmology *cosmo, float v_new1) {

  p->v[0] = v_new1 * cosmo->a_inv;
}

__attribute__((always_inline)) INLINE static void hydro_set_velocity2(
    struct part *restrict p, const struct cosmology *cosmo, float v_new2) {

  p->v[1] = v_new2 * cosmo->a_inv;
}

__attribute__((always_inline)) INLINE static void hydro_set_velocity3(
    struct part *restrict p, const struct cosmology *cosmo, float v_new3) {

  p->v[2] = v_new3 * cosmo->a_inv;
}

__attribute__((always_inline)) INLINE static void
hydro_set_physical_velocity1(struct part *p, struct xpart *xp,
                                   const struct cosmology *cosmo,
                                   const float v_new1) {

  xp->v_full[0] = v_new1 * cosmo->a_inv;
  p->gpart->v_full[0] = xp->v_full[0];
}

__attribute__((always_inline)) INLINE static void
hydro_set_physical_velocity2(struct part *p, struct xpart *xp,
                                   const struct cosmology *cosmo,
                                   const float v_new2) {

  xp->v_full[1] = v_new2 * cosmo->a_inv;
  p->gpart->v_full[1] = xp->v_full[1];
}

__attribute__((always_inline)) INLINE static void
hydro_set_physical_velocity3(struct part *p, struct xpart *xp,
                                   const struct cosmology *cosmo,
                                   const float v_new3) {

  xp->v_full[2] = v_new3 * cosmo->a_inv;
  p->gpart->v_full[2] = xp->v_full[2];
}


/**
 * @brief Update the value of the diffusive coefficients to the
 *        feedback reset value for the scheme.
 *
 * @param p the particle of interest
 */
__attribute__((always_inline)) INLINE static void
hydro_diffusive_feedback_reset(struct part *restrict p) {
  /* This scheme has fixed alpha */
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
                       (cosmo->a_factor_sound_speed * p->force.v_sig);

  return dt_cfl;
}

/**
 * @brief returns the signal velocity
 *
 * @brief p  the particle
 */
__attribute__((always_inline)) INLINE static float hydro_get_signal_velocity(
    const struct part *restrict p) {

  return p->force.v_sig;
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
  p->density.rho_dh = 0.f;
  p->num_unkicked_ngbs = 0;

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  p->N_density = 1; /* Self contribution */
  p->N_force = 0;
  p->N_density_exact = 0;
  p->N_force_exact = 0;
  p->rho_exact = 0.f;
  p->n_density = 0.f;
  p->n_density_exact = 0.f;
  p->n_force = 0.f;
  p->n_force_exact = 0.f;
  p->inhibited_exact = 0;
  p->limited_part = 0;
#endif

  hydro_init_part_extra_kernel(p);
  hydro_init_part_extra_viscosity(p);
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

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  p->n_density += kernel_root;
  p->n_density *= h_inv_dim;
#endif

  hydro_end_density_extra_viscosity(p);
  hydro_end_density_extra_kernel(p);

  p->rho_sph = p->rho; 

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

  /* Re-set problematic values */
  p->rho = p->mass * kernel_root * h_inv_dim;
  p->density.wcount = kernel_root * h_inv_dim;
  p->density.rho_dh = 0.f;
  p->density.wcount_dh = 0.f;

  p->is_h_max = 1;
  p->m0 = p->mass * kernel_root * h_inv_dim / p->rho_evol;
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
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

    p->max_id = hydro_props->max_id;
    
  if (p->h > 0.999f * hydro_props->h_max) {
    p->is_h_max = 1;
  } else {
    p->is_h_max = 0;
  }

  if (p->id < hydro_props->max_id && p->hit_by_jet_feedback > 0 && p->num_unkicked_ngbs == 0) {
    p->id = p->id + hydro_props->max_id;

    p->rho = p->rho_sph; 
    p->rho_evol = p->rho_sph; 
    xp->rho_evol_full = p->rho_sph; 
    p->grad_m0[0] = 0.f;
    p->grad_m0[1] = 0.f;
    p->grad_m0[2] = 0.f;
    p->m0 = 1.f;

  }
  
if (p->id > p->max_id) {
  hydro_prepare_gradient_extra_kernel(p);
  hydro_prepare_gradient_extra_viscosity(p);
}
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
    struct part *restrict p) {}

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
if (p->id > p->max_id) {
  hydro_end_gradient_extra_kernel(p);
  hydro_end_gradient_extra_viscosity(p);
}

  // Set the density to be used in the force loop to be the evolved density
  p->rho = p->rho_evol;
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
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor, const float dt_alpha,
    const float dt_therm) {
    
if (p->id > p->max_id) {
  hydro_prepare_force_extra_kernel(p);
}

#ifdef PLANETARY_FIXED_ENTROPY
  /* Override the internal energy to satisfy the fixed entropy */
  p->u = gas_internal_energy_from_entropy(p->rho_evol, p->s_fixed, p->mat_id);
  xp->u_full = p->u;
#endif

  /* Compute the sound speed */
  const float soundspeed =
      gas_soundspeed_from_internal_energy(p->rho_evol, p->u, p->mat_id);

  float div_v = p->dv_norm_kernel[0][0] + p->dv_norm_kernel[1][1] +
                p->dv_norm_kernel[2][2];
  float curl_v[3];
  curl_v[0] = p->dv_norm_kernel[1][2] - p->dv_norm_kernel[2][1];
  curl_v[1] = p->dv_norm_kernel[2][0] - p->dv_norm_kernel[0][2];
  curl_v[2] = p->dv_norm_kernel[0][1] - p->dv_norm_kernel[1][0];
  float mod_curl_v = sqrtf(curl_v[0] * curl_v[0] + curl_v[1] * curl_v[1] +
                           curl_v[2] * curl_v[2]);

  // Balsara switch using normalised kernel gradients
  float balsara;
  if (div_v == 0.f) {
    balsara = 0.f;
  } else {
    balsara = fabsf(div_v) /
              (fabsf(div_v) + mod_curl_v + 0.0001f * soundspeed / p->h);
  }

  /* Compute the pressure */
  const float pressure =
      gas_pressure_from_internal_energy(p->rho_evol, p->u, p->mat_id);
  p->force.pressure = pressure;
  p->force.soundspeed = soundspeed;
  p->force.balsara = balsara;
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
  p->drho_dt = 0.0f;
  p->force.h_dt = 0.0f;
  p->force.v_sig = p->force.soundspeed;
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
    struct part *restrict p, const struct xpart *restrict xp,
    const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor) {

  /* Re-set the predicted velocities */
  p->v[0] = xp->v_full[0];
  p->v[1] = xp->v_full[1];
  p->v[2] = xp->v_full[2];

  /* Re-set the internal energy */
  p->u = xp->u_full;

  p->rho = xp->rho_evol_full;
  p->rho_evol = xp->rho_evol_full;

  /* Compute the pressure */
  const float pressure =
      gas_pressure_from_internal_energy(p->rho_evol, p->u, p->mat_id);

  /* Compute the sound speed */
  const float soundspeed =
      gas_soundspeed_from_internal_energy(p->rho_evol, p->u, p->mat_id);

  p->force.pressure = pressure;
  p->force.soundspeed = soundspeed;

  p->force.v_sig = max(p->force.v_sig, 2.f * soundspeed);
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
 * @param hydro_props The constants used in the scheme
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part *restrict p, const struct xpart *restrict xp, float dt_drift,
    float dt_therm, float dt_kick_grav, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Predict the internal energy and density */
  p->u += p->u_dt * dt_therm;
  p->rho_evol += p->drho_dt * dt_therm;

  /* compute minimum SPH quantities */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  const float floor_rho = p->mass * kernel_root * h_inv_dim;
  p->rho_evol = max(p->rho_evol, floor_rho);

  /* Check against absolute minimum */
  const float min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  p->u = max(p->u, min_u);

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt_drift;
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

  p->rho = p->rho_evol;

  const float floor_u = FLT_MIN;
  p->u = max(p->u, floor_u);

  /* Compute the new pressure */
  const float pressure =
      gas_pressure_from_internal_energy(p->rho_evol, p->u, p->mat_id);

  /* Compute the new sound speed */
  const float soundspeed =
      gas_soundspeed_from_internal_energy(p->rho_evol, p->u, p->mat_id);

  p->force.pressure = pressure;
  p->force.soundspeed = soundspeed;

  p->force.v_sig = max(p->force.v_sig, 2.f * soundspeed);
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
    const struct entropy_floor_properties *floor_props,const double time) {

  /* Integrate the internal energy forward in time */
  const float delta_u = p->u_dt * dt_therm;

  /* Do not decrease the energy by more than a factor of 2*/
  xp->u_full = max(xp->u_full + delta_u, 0.5f * xp->u_full);

  /* Check against absolute minimum */
  const float min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;
  const float floor_u = FLT_MIN;

  /* Take highest of both limits */
  const float energy_min = max(min_u, floor_u);

  if (xp->u_full < energy_min) {
    xp->u_full = energy_min;
    p->u_dt = 0.f;
  }

  const float delta_rho = p->drho_dt * dt_therm;

  xp->rho_evol_full =
      max(xp->rho_evol_full + delta_rho, 0.5f * xp->rho_evol_full);

  /* Minimum SPH quantities */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  const float floor_rho = p->mass * kernel_root * h_inv_dim;
  if (xp->rho_evol_full < floor_rho) {
    xp->rho_evol_full = floor_rho;
    p->drho_dt = 0.f;
  }

  if (hydro_props->use_2d==1) {
        float aspect_ratio = 1.;
    if (hydro_props->constant_density==1)
      aspect_ratio = 0.3;

    float direction = 0.;
    double delta_x = (p->x[0]-aspect_ratio*hydro_props->box_centre);
    double delta_y = (p->x[1]-hydro_props->box_centre);
    double r = sqrt(delta_x*delta_x + delta_y*delta_y);

    double cos_theta = delta_y/r; /* 1.-(1.-fabs(delta_z/R))*(1.-cos(hydro_props->opening_angle/180.*M_PI)) */
    double sin_theta = sqrt(1.-cos_theta*cos_theta);

    if (hydro_props->use_jets==1 && time<hydro_props->jet_duration && p->id<(hydro_props->jet_power)*time/(0.5*p->mass*hydro_props->v_jet*hydro_props->v_jet) && p->hit_by_jet_feedback<1) {

      if (delta_y > 0.) {
        direction = 1;
      } else {
        direction = -1;
      }

      if (hydro_props->launch_spread==1) {
        float tan_theta = tan(hydro_props->opening_angle/180.*M_PI)*sqrt(1.-fabs(cos_theta));
        cos_theta = direction / sqrt(1.+tan_theta*tan_theta);
        sin_theta = tan_theta / sqrt(1.+tan_theta*tan_theta); 
        /* cos_theta = direction * ((1.-cos(hydro_props->opening_angle/180.*M_PI))*fabs(cos_theta)+1.);
        sin_theta = sqrt(1.-cos_theta*cos_theta); */
      }

      float vel_kick = 0.;
      if (hydro_props->assume_equipartition==1) {
        vel_kick = 5./sqrt(26) * hydro_props->v_jet;
      } else {
        vel_kick = hydro_props->v_jet;
      }

      
      if (delta_x>0) {
        sin_theta = 1.f * sin_theta;
      }
      if (delta_x<0) {
        sin_theta = -1.*sin_theta;
      }
      float vel_kick_vec[2];
      if (hydro_props->launch_parallel==1) {
        vel_kick_vec[0] = 0.;
        vel_kick_vec[1] = 0.;
      } else {
        vel_kick_vec[0] = vel_kick*sin_theta;
        vel_kick_vec[1] = vel_kick*cos_theta;
      }

      float v_new[2];

      v_new[0] = vel_kick_vec[0];
      v_new[1] = vel_kick_vec[1];

      printf("New velocity in z direction: %f\n", v_new[1]);

      p->v[0]=v_new[0];
      p->v[1]=v_new[1];
      xp->v_full[0]=v_new[0];
      xp->v_full[1]=v_new[1];
      
      /* hydro_set_physical_velocity1(p, xp, cosmo, v_new[0]);
      hydro_set_physical_velocity2(p, xp, cosmo, v_new[1]); */
      p->hit_by_jet_feedback = 1;

      if (hydro_props->assume_equipartition==1) {
        const double u_init_jet = hydro_get_physical_internal_energy(p, xp, cosmo);
        const float delta_u_jet = 1./26. * 0.5 * vel_kick * vel_kick;
        const double u_new_jet = u_init_jet + delta_u_jet;

        hydro_set_physical_internal_energy(p, xp, cosmo, u_new_jet);
        hydro_set_drifted_physical_internal_energy(p, cosmo, u_new_jet);
      }

      float v_norm = sqrt(v_new[0]*v_new[0]+v_new[1]*v_new[1]);
      hydro_diffusive_feedback_reset(p); 
      hydro_set_v_sig_based_on_velocity_kick(p, cosmo, v_norm);
      timestep_sync_part(p);
      /* p->timestep_counter = 1; */

      if (hydro_props->launch_spread==0) {
        p->id = p->id + hydro_props->max_id; 
      }
    }
  } 
  if (hydro_props->use_2d==0) {
    float aspect_ratio = 1.;
    float use_full_box = 1.;
    if (hydro_props->constant_density==1) {
      aspect_ratio = 0.3;
      use_full_box = 1.;
    } 
      
    float direction = 0.;
    double delta_x = (p->x[0]-aspect_ratio*hydro_props->box_centre);
    double delta_y = (p->x[1]-aspect_ratio*hydro_props->box_centre);
    double delta_z = (p->x[2]-use_full_box * hydro_props->box_centre);
    double r = sqrt(delta_x*delta_x + delta_y*delta_y);
    double R = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
    
    double phi = 0.;
    if (delta_y>0.) {
      phi = acos(delta_x/r);
    } else if (delta_y<0.) {
      phi = -1.*acos(delta_x/r);
    } else if (delta_y==0.) {
      if (delta_x>0.) {
        phi = 0.;
      } else {
        phi = 180.;
      }
    }

    double cos_theta = delta_z/R; /* 1.-(1.-fabs(delta_z/R))*(1.-cos(hydro_props->opening_angle/180.*M_PI)) */
    double sin_theta = sqrt(1.-cos_theta*cos_theta);
    
     /*long long id_to_check = p->id;
    if (id_to_check%2!=0) {
      id_to_check = id_to_check + 1;
    } */
    if (hydro_props->use_jets==1 && time<hydro_props->jet_duration && p->id<(hydro_props->jet_power*time)/(0.5*p->mass*hydro_props->v_jet*hydro_props->v_jet) && p->hit_by_jet_feedback<1) {
      if (delta_z > 0.) {
        direction = 1;
      } else {
        direction = -1;
      }
      if (hydro_props->launch_spread==1) {
        float tan_theta = tan(hydro_props->opening_angle/180.*M_PI)*sqrt(1.-fabs(cos_theta));
        cos_theta = direction / sqrt(1.+tan_theta*tan_theta);
        sin_theta = tan_theta / sqrt(1.+tan_theta*tan_theta); 
        /* cos_theta = direction * ((1.-cos(hydro_props->opening_angle/180.*M_PI))*fabs(cos_theta)+1.);
        sin_theta = sqrt(1.-cos_theta*cos_theta); */
      }
      float vel_kick = 0.;
      if (hydro_props->assume_equipartition==1) {
        vel_kick = 5./sqrt(26) * hydro_props->v_jet;
      } else {
        vel_kick = hydro_props->v_jet;
      }
      float vel_kick_vec[3];
      if (hydro_props->launch_parallel==1) {
        vel_kick_vec[0] = 0.;
        vel_kick_vec[1] = 0.;
        vel_kick_vec[2] = direction*vel_kick;
      } else {
        vel_kick_vec[0] = vel_kick*sin_theta*cos(phi);
        vel_kick_vec[1] = vel_kick*sin_theta*sin(phi);
        vel_kick_vec[2] = vel_kick*cos_theta;
      }
      float v_new[3];

      v_new[0] = vel_kick_vec[0];
      v_new[1] = vel_kick_vec[1];
      v_new[2] = vel_kick_vec[2];

      printf("Height: %f\n",p->x[2]);
      printf("New velocity in z direction: %f\n", v_new[2]);

      hydro_set_velocity1(p, cosmo, v_new[0]);
      hydro_set_velocity2(p, cosmo, v_new[1]);
      hydro_set_velocity3(p, cosmo, v_new[2]);

      hydro_set_physical_velocity1(p, xp, cosmo, v_new[0]);
      hydro_set_physical_velocity2(p, xp, cosmo, v_new[1]);
      hydro_set_physical_velocity3(p, xp, cosmo, v_new[2]);
      p->hit_by_jet_feedback = 1;

      if (hydro_props->assume_equipartition==1) {
        const double u_init_jet = hydro_get_physical_internal_energy(p, xp, cosmo);
        const float delta_u_jet = 1./26. * 0.5 * vel_kick * vel_kick;
        const double u_new_jet = u_init_jet + delta_u_jet;

        hydro_set_physical_internal_energy(p, xp, cosmo, u_new_jet);
        hydro_set_drifted_physical_internal_energy(p, cosmo, u_new_jet);
      }

      double v_norm = sqrt(v_new[0]*v_new[0] + v_new[1]*v_new[1] + v_new[2]*v_new[2]);
      hydro_diffusive_feedback_reset(p); 
      hydro_set_v_sig_based_on_velocity_kick(p, cosmo, v_norm);
      timestep_sync_part(p);
      p->timestep_counter = 1;
      /* p->id = p->id + hydro_props->max_id; */

    }
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
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Compute the pressure */
  const float pressure =
      gas_pressure_from_internal_energy(p->rho_evol, p->u, p->mat_id);

  /* Compute the sound speed */
  const float soundspeed =
      gas_soundspeed_from_internal_energy(p->rho_evol, p->u, p->mat_id);

  p->force.pressure = pressure;
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

  p->num_unkicked_ngbs = 0;
    
  p->timestep_counter = 0;
  p->hit_by_jet_feedback = 0;

  p->rho_evol = p->rho;
  xp->rho_evol_full = p->rho_evol;

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

/**
 * @brief Operations performed when a particle gets removed from the
 * simulation volume.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @param time The simulation time.
 */
__attribute__((always_inline)) INLINE static void hydro_remove_part(
    const struct part *p, const struct xpart *xp, const double time) {

  /* Print the particle info as csv to facilitate later analysis, e.g. with
   * grep '## Removed' -A 1 --no-group-separator output.txt > removed.txt
   */
  printf(
      "## Removed particle: "
      "id, x, y, z, vx, vy, vz, m, u, P, rho, h, mat_id, time \n"
      "%lld, %.7g, %.7g, %.7g, %.7g, %.7g, %.7g, "
      "%.7g, %.7g, %.7g, %.7g, %.7g, %d, %.7g \n",
      p->id, p->x[0], p->x[1], p->x[2], p->v[0], p->v[1], p->v[2], p->mass,
      p->u, p->force.pressure, p->rho, p->h, p->mat_id, time);
}

#endif /* SWIFT_PLANETARY_HYDRO_H */
