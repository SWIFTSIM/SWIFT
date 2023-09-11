/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2015 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

  return gas_internal_energy_from_entropy(p->rho, xp->entropy_full);
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

  return gas_internal_energy_from_entropy(p->rho * cosmo->a3_inv,
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

  return gas_internal_energy_from_entropy(p->rho, p->entropy);
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

  return gas_internal_energy_from_entropy(p->rho * cosmo->a3_inv, p->entropy);
}

/**
 * @brief Returns the comoving pressure of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_pressure(
    const struct part *restrict p) {

  return gas_pressure_from_entropy(p->rho, p->entropy);
}

/**
 * @brief Returns the physical pressure of a particle
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_pressure(
    const struct part *restrict p, const struct cosmology *cosmo) {

  return gas_pressure_from_entropy(p->rho * cosmo->a3_inv, p->entropy);
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
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
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
 * @brief Returns the time derivative of co-moving internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy_dt(const struct part *restrict p) {

  return gas_internal_energy_from_entropy(p->rho, p->entropy_dt);
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

  return gas_internal_energy_from_entropy(p->rho * cosmo->a3_inv,
                                          p->entropy_dt);
}

/**
 * @brief Sets the time derivative of the co-moving internal energy of a
 * particle
 *
 * We assume a constant density for the conversion to entropy.
 *
 * @param p The particle of interest.
 * @param du_dt The new time derivative of the comoving internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_comoving_internal_energy_dt(struct part *restrict p,
                                      const float du_dt) {

  p->entropy_dt = gas_entropy_from_internal_energy(p->rho, du_dt);
}

/**
 * @brief Sets the time derivative of the physical internal energy of a particle
 *
 * We assume a constant density for the conversion to entropy.
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param du_dt The time derivative of the physical internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy_dt(struct part *restrict p,
                                      const struct cosmology *restrict cosmo,
                                      const float du_dt) {
  p->entropy_dt =
      gas_entropy_from_internal_energy(p->rho * cosmo->a3_inv, du_dt);
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

  xp->entropy_full =
      gas_entropy_from_internal_energy(p->rho * cosmo->a3_inv, u);
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

  p->entropy = gas_entropy_from_internal_energy(p->rho * cosmo->a3_inv, u);

  /* Now recompute the extra quantities */

  /* Inverse of the co-moving density */
  const float rho_inv = 1.f / p->rho;

  /* Compute the pressure */
  float comoving_pressure = gas_pressure_from_entropy(p->rho, p->entropy);
  comoving_pressure = pressure_floor_get_comoving_pressure(
      p, pressure_floor, comoving_pressure, cosmo);

  /* Compute the sound speed */
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, comoving_pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float P_over_rho2 = comoving_pressure * rho_inv * rho_inv;

  /* Update variables. */
  p->force.P_over_rho2 = P_over_rho2;
  p->force.soundspeed = soundspeed;

  p->force.v_sig = max(p->force.v_sig, 2.f * soundspeed);
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
  p->force.v_sig =
      max(2.f * soundspeed, p->force.v_sig + const_viscosity_beta * dv);
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
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 * @param hydro_properties The constants used in the scheme
 * @param cosmo The cosmological model.
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
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const float mu_ij, const float beta) {

  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

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
    struct part *p, float dt) {}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various density tasks
 *
 * @param p The particle to act upon
 * @param hs #hydro_space containing hydro specific space information.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part *restrict p, const struct hydro_space *hs) {

#ifdef DEBUG_INTERACTIONS_SPH
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS; ++i) p->ids_ngbs_density[i] = -1;
  p->num_ngb_density = 0;
#endif

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
  p->density.rho_dh -= hydro_dimension * p->mass * kernel_root;
  p->density.wcount += kernel_root;
  p->density.wcount_dh -= hydro_dimension * kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  p->rho *= h_inv_dim;
  p->density.rho_dh *= h_inv_dim_plus_one;
  p->density.wcount *= h_inv_dim;
  p->density.wcount_dh *= h_inv_dim_plus_one;

  const float rho_inv = 1.f / p->rho;
  const float a_inv2 = cosmo->a2_inv;

  /* Finish calculation of the (physical) velocity curl components */
  p->density.rot_v[0] *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  p->density.rot_v[1] *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  p->density.rot_v[2] *= h_inv_dim_plus_one * a_inv2 * rho_inv;

  /* Finish calculation of the (physical) velocity divergence */
  p->density.div_v *= h_inv_dim_plus_one * a_inv2 * rho_inv;
}

/**
 * @brief Prepare a particle for the gradient calculation.
 *
 * This function is called after the density loop and before the gradient loop.
 * Nothing to do in this scheme as the gradient loop is not used.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 * @param hydro_props Hydrodynamic properties.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_gradient(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props) {}

/**
 * @brief Resets the variables that are required for a gradient calculation.
 *
 * This function is called after hydro_prepare_gradient.
 * Nothing to do in this scheme as the gradient loop is not used.
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
 * Nothing to do in this scheme as the gradient loop is not used.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_end_gradient(
    struct part *p) {}

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

  warning(
      "Gas particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      p->id, h, p->density.wcount);

  /* Re-set problematic values */
  p->rho = p->mass * kernel_root * h_inv_dim;
  p->density.wcount = kernel_root * h_inv_dim;
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

  const float fac_Balsara_eps = cosmo->a_factor_Balsara_eps;

  /* Inverse of the co-moving density */
  const float rho_inv = 1.f / p->rho;

  /* Inverse of the smoothing length */
  const float h_inv = 1.f / p->h;

  /* Compute the norm of the curl */
  const float curl_v = sqrtf(p->density.rot_v[0] * p->density.rot_v[0] +
                             p->density.rot_v[1] * p->density.rot_v[1] +
                             p->density.rot_v[2] * p->density.rot_v[2]);

  /* Compute the norm of div v including the Hubble flow term */
  const float div_physical_v = p->density.div_v + hydro_dimension * cosmo->H;
  const float abs_div_physical_v = fabsf(div_physical_v);

  /* Compute the pressure */
  float comoving_pressure = gas_pressure_from_entropy(p->rho, p->entropy);
  comoving_pressure = pressure_floor_get_comoving_pressure(
      p, pressure_floor, comoving_pressure, cosmo);

  /* Compute the sound speed */
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, comoving_pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float P_over_rho2 = comoving_pressure * rho_inv * rho_inv;

  /* Compute the Balsara switch */
  /* Pre-multiply in the AV factor; hydro_props are not passed to the iact
   * functions */
  const float balsara = hydro_props->viscosity.alpha * abs_div_physical_v /
                        (abs_div_physical_v + curl_v +
                         0.0001f * fac_Balsara_eps * soundspeed * h_inv);

  /* Compute the "grad h" term */
  float rho_dh = p->density.rho_dh;
  /* Ignore changing-kernel effects when h ~= h_max */
  if (p->h > 0.9999f * hydro_props->h_max) {
    rho_dh = 0.f;
    warning("h ~ h_max for particle with ID %lld (h: %g)", p->id, p->h);
  }
  const float grad_rho_term = hydro_dimension_inv * p->h * rho_dh * rho_inv;
  float omega_inv;
  if (grad_rho_term < -0.9999f) {
    omega_inv = 1.f;
    warning(
        "grad_rho_term very small for particle with ID %lld (h: %g, rho: %g, "
        "rho_dh: %g).",
        p->id, p->h, p->rho, rho_dh);
  } else {
    omega_inv = 1.f / (1.f + grad_rho_term);
  }

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

#ifdef DEBUG_INTERACTIONS_SPH
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS; ++i) p->ids_ngbs_force[i] = -1;
  p->num_ngb_force = 0;
#endif

  /* Reset the acceleration. */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* Reset the time derivatives. */
  p->entropy_dt = 0.0f;
  p->force.h_dt = 0.0f;

  /* Reset maximal signal velocity */
  p->force.v_sig = 2.f * p->force.soundspeed;
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
  p->entropy = xp->entropy_full;

  /* Re-compute the pressure */
  float comoving_pressure = gas_pressure_from_entropy(p->rho, p->entropy);
  comoving_pressure = pressure_floor_get_comoving_pressure(
      p, pressure_floor, comoving_pressure, cosmo);

  /* Compute the new sound speed */
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, comoving_pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_inv = 1.f / p->rho;
  const float P_over_rho2 = comoving_pressure * rho_inv * rho_inv;

  /* Update variables */
  p->force.soundspeed = soundspeed;
  p->force.P_over_rho2 = P_over_rho2;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param p The particle
 * @param xp The extended data of the particle
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

  /* Predict the entropy */
  p->entropy += p->entropy_dt * dt_therm;

  const float h_inv = 1.f / p->h;

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt_drift;
  if (fabsf(w1) < 0.2f) {
    p->h *= approx_expf(w1); /* 4th order expansion of exp(w) */
  } else {
    p->h *= expf(w1);
  }

  /* Predict density */
  const float w2 = -hydro_dimension * w1;
  if (fabsf(w2) < 0.2f) {
    p->rho *= approx_expf(w2); /* 4th order expansion of exp(w) */
  } else {
    p->rho *= expf(w2);
  }

  /* Check against entropy floor */
  const float floor_A = entropy_floor(p, cosmo, floor_props);

  /* Check against absolute minimum; recall that A_physical == A_comoving
   * by definition so no conversion is necessary */
  const float min_u_physical = hydro_props->minimal_internal_energy;
  /* Conversion done in physical space; multiplication by a3_inv is faster
   * than division by a_factor_internal_energy to do in comoving space */
  const float min_A =
      gas_entropy_from_internal_energy(p->rho * cosmo->a3_inv, min_u_physical);

  p->entropy = max(p->entropy, floor_A);
  p->entropy = max(p->entropy, min_A);

  /* Re-compute the pressure */
  float comoving_pressure = gas_pressure_from_entropy(p->rho, p->entropy);
  comoving_pressure = pressure_floor_get_comoving_pressure(
      p, pressure_floor, comoving_pressure, cosmo);

  /* Compute the new sound speed */
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, comoving_pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_inv = 1.f / p->rho;
  const float P_over_rho2 = comoving_pressure * rho_inv * rho_inv;

  /* Update variables */
  p->force.soundspeed = soundspeed;
  p->force.P_over_rho2 = P_over_rho2;

  p->force.v_sig = max(p->force.v_sig, 2.f * soundspeed);
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

  p->entropy_dt =
      0.5f * gas_entropy_from_internal_energy(p->rho, p->entropy_dt);
}

/**
 * @brief Kick the additional variables
 *
 * @param p The particle to act upon
 * @param xp The particle extended data to act upon
 * @param dt_therm The time-step for this kick (for thermodynamic quantities)
 * @param dt_grav The time-step for this kick (for gravity forces)
 * @param dt_grav_mesh The time-step for this kick (mesh gravity).
 * @param dt_hydro The time-step for this kick (for hydro forces)
 * @param dt_kick_corr The time-step for this kick (for correction of the kick)
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part *restrict p, struct xpart *restrict xp, float dt_therm,
    float dt_grav, float dt_grav_mesh, float dt_hydro, float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Integrate the entropy forward in time */
  const float delta_entropy = p->entropy_dt * dt_therm;

  /* Do not decrease the entropy by more than a factor of 2 */
  xp->entropy_full =
      max(xp->entropy_full + delta_entropy, 0.5f * xp->entropy_full);

  /* Check against entropy floor */
  const float floor_A = entropy_floor(p, cosmo, floor_props);

  /* Check against absolute minimum; recall that A_physical == A_comoving
   * by definition so no conversion is necessary */
  const float min_u_physical = hydro_props->minimal_internal_energy;
  /* Conversion done in physical space; multiplication by a3_inv is faster
   * than division by a_factor_internal_energy to do in comoving space */
  const float min_A =
      gas_entropy_from_internal_energy(p->rho * cosmo->a3_inv, min_u_physical);

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
 * @param p The particle to act upon.
 * @param xp The extended data.
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param pressure_floor The properties of the pressure floor.
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

  /* We read u in the entropy field. We now get (comoving) A from (physical) u
   * and (physical) rho. Note that comoving A (A') == physical A */
  xp->entropy_full =
      gas_entropy_from_internal_energy(p->rho * cosmo->a3_inv, p->entropy);
  p->entropy = xp->entropy_full;

  /* Apply the minimal energy limit */
  const float physical_density = p->rho * cosmo->a3_inv;
  const float min_physical_energy = hydro_props->minimal_internal_energy;
  const float min_physical_entropy =
      gas_entropy_from_internal_energy(physical_density, min_physical_energy);
  const float min_comoving_entropy = min_physical_entropy; /* A' = A */
  if (xp->entropy_full < min_comoving_entropy) {
    xp->entropy_full = min_comoving_entropy;
    p->entropy = min_comoving_entropy;
    p->entropy_dt = 0.f;
  }

  /* Compute the pressure */
  float comoving_pressure = gas_pressure_from_entropy(p->rho, p->entropy);
  comoving_pressure = pressure_floor_get_comoving_pressure(
      p, pressure_floor, comoving_pressure, cosmo);

  /* Compute the sound speed */
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, comoving_pressure);

  /* Divide the pressure by the density squared to get the SPH term */
  const float rho_inv = 1.f / p->rho;
  const float P_over_rho2 = comoving_pressure * rho_inv * rho_inv;

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
 * @param time The simulation time.
 */
__attribute__((always_inline)) INLINE static void hydro_remove_part(
    const struct part *p, const struct xpart *xp, const double time) {}

#endif /* SWIFT_GADGET2_HYDRO_H */
