/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_MAGMA_HYDRO_H
#define SWIFT_MAGMA_HYDRO_H

/**
 * @file MAGMA/hydro.h
 * @brief MAGMA2 implementation of SPH (Non-neighbour loop equations)
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

/* System includes */
#include <string.h>
#include <float.h>

/**
 * @brief Returns the comoving internal energy of a particle at the last
 * time the particle was kicked.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy(const struct part *p,
                                   const struct xpart *xp) {

  return xp->u_full;
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
hydro_get_physical_internal_energy(const struct part *p, const struct xpart *xp,
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
hydro_get_drifted_comoving_internal_energy(const struct part *p) {

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
hydro_get_drifted_physical_internal_energy(const struct part *p,
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
    const struct part *p) {

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
    const struct part *p, const struct cosmology *cosmo) {

  return cosmo->a_factor_pressure *
         gas_pressure_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the comoving entropy of a particle at the last
 * time the particle was kicked.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_entropy(
    const struct part *p, const struct xpart *xp) {

  return gas_entropy_from_internal_energy(p->rho, xp->u_full);
}

/**
 * @brief Returns the physical entropy of a particle at the last
 * time the particle was kicked.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_entropy(
    const struct part *p, const struct xpart *xp,
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
hydro_get_drifted_comoving_entropy(const struct part *p) {

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
hydro_get_drifted_physical_entropy(const struct part *p,
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
hydro_get_comoving_soundspeed(const struct part *p) {

  return p->force.soundspeed;
}

/**
 * @brief Returns the physical sound speed of a particle
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_soundspeed(const struct part *p,
                              const struct cosmology *cosmo) {

  return cosmo->a_factor_sound_speed * p->force.soundspeed;
}

/**
 * @brief Returns the comoving density of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_density(
    const struct part *p) {

  return p->rho;
}

/**
 * @brief Returns the comoving density of a particle.
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_density(
    const struct part *p, const struct cosmology *cosmo) {

  return cosmo->a3_inv * p->rho;
}

/**
 * @brief Returns the mass of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_mass(
    const struct part *p) {

  return p->mass;
}

/**
 * @brief Sets the mass of a particle
 *
 * @param p The particle of interest
 * @param m The mass to set.
 */
__attribute__((always_inline)) INLINE static void hydro_set_mass(struct part *p,
                                                                 float m) {

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
hydro_get_comoving_internal_energy_dt(const struct part *p) {

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
hydro_get_physical_internal_energy_dt(const struct part *p,
                                      const struct cosmology *cosmo) {

  return p->u_dt * cosmo->a_factor_internal_energy;
}

/**
 * @brief Sets the time derivative of the co-moving internal energy of a
 * particle
 *
 * We assume a constant density for the conversion to entropy.
 *
 * @param p The particle of interest.
 * @param du_dt The new time derivative of the internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_comoving_internal_energy_dt(struct part *p, float const du_dt) {

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
hydro_set_physical_internal_energy_dt(struct part *p,
                                      const struct cosmology *cosmo,
                                      const float du_dt) {

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
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(
    struct part *p, const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor, const float u) {

  p->u = u / cosmo->a_factor_internal_energy;

  /* Now recompute the extra quantities */

  /* Compute the sound speed */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Update variables. */
  p->force.pressure = pressure;
  p->force.soundspeed = soundspeed;
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
                                       const float dv_phys) {}

/**
 * @brief Update the value of the viscosity alpha for the scheme.
 *
 * @param p the particle of interest
 * @param alpha the new value for the viscosity coefficient.
 */
__attribute__((always_inline)) INLINE static void hydro_set_viscosity_alpha(
    struct part *p, const float alpha) {
  /* This scheme has fixed alpha */
}

/**
 * @brief Update the value of the diffusive coefficients to the
 *        feedback reset value for the scheme.
 *
 * @param p the particle of interest
 */
__attribute__((always_inline)) INLINE static void
hydro_diffusive_feedback_reset(struct part *p) {
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
    const struct part *p, const struct xpart *xp,
    const struct hydro_props *hydro_properties, const struct cosmology *cosmo) {

  const float CFL_condition = hydro_properties->CFL_condition;

  /* Criterion based on acceleration (eq. 35) */
  const float norm_a = p->a_hydro[0] * p->a_hydro[0] +
                       p->a_hydro[1] * p->a_hydro[1] +
                       p->a_hydro[2] * p->a_hydro[2];
  const float dt_acc = sqrtf(p->h / sqrtf(norm_a));

  /* Criterion based on acceleration (eq. 35) */
  const float c = p->force.soundspeed;
  float dt_Courant =
      p->h / (c + 0.6f * const_viscosity_alpha * (c + 2.f * p->force.mu_tilde));
  
  /* MHD variables ----------------------------------------------*/
  
  /* Retrieve and compute relevant cosmological factors */
  const float a = cosmo->a;
  const float a2 = a * a;

  /* Retrieve relevant particle attributes */
  const float rho = p->rho;

  float B_over_rho[3];
  float B_over_rho_dt[3];
  for (int k = 0; k < 3; k++) {
    B_over_rho[k] = p->mhd.B_over_rho[k];
    B_over_rho_dt[k] = p->mhd.B_over_rho_dt[k];
  }

  const float psi_over_ch = p->mhd.psi_over_ch;
  const float psi_over_ch_dt = p->mhd.psi_over_ch_dt;

  /* Compute the norm squared of vectors of interest */
  float B_over_rho2 = 0.0f;
  for (int k = 0; k < 3; k++) {
    B_over_rho2 += B_over_rho[k] * B_over_rho[k];
  }

  /* Compute metric to evaluate dynamical significance of Dedner scalar field */
  const float vpsi_tp_vB =
      B_over_rho2 ? fabsf(psi_over_ch) / sqrtf(B_over_rho2 * rho * rho) : 0.0f;

  /* Condition to limit the per time-step change in the magnitude of the
   * magnetic field */
  const float maxRelChangeBoverRho = hydro_properties->mhd.maxRelChangeBoverRho;

  float denum_dt_deltaB2 = 0.f;
  for (int k = 0; k < 3; k++) {
    denum_dt_deltaB2 += B_over_rho_dt[k] * B_over_rho_dt[k];
  }

  const float dt_deltaB =
      B_over_rho2 && denum_dt_deltaB2
          ? maxRelChangeBoverRho * a2 * sqrtf(B_over_rho2 / denum_dt_deltaB2)
          : FLT_MAX;

  /* Condition to limit the per time-step change in the magnitude of the Dedner
   * scalar field */
  const float maxRelChangePsiOverCh =
      hydro_properties->mhd.maxRelChangePsiOverCh;
  const float R_ePsi_to_eB = hydro_properties->mhd.R_ePsi_to_eB;

  const float denum_dt_deltaPsi = fabsf(psi_over_ch_dt);

  const float dt_deltaPsi =
      (vpsi_tp_vB > R_ePsi_to_eB) && (denum_dt_deltaPsi != 0.0f)
          ? maxRelChangePsiOverCh * a2 * fabsf(psi_over_ch) / denum_dt_deltaPsi
          : FLT_MAX;

  /* Keep the minimum of two preious conditions */
  const float dt_deltaField = fminf(dt_deltaB, dt_deltaPsi);

  /* Compute new time-step because of physical diffusion */
  const float dt_eta = p->mhd.resistive_eta != 0.f
                           ? a * a * p->h * p->h /
                                 p->mhd.resistive_eta
                           : FLT_MAX;

  /* Keep the minimum of all MHD time-steps */
  const float dt_mhd = fminf(dt_deltaField, dt_eta);
  
  dt_Courant = fminf(dt_Courant, dt_mhd);

  /* END MHD variables -------------------------------------------*/

  return CFL_condition * fminf(dt_acc, dt_Courant);
}

/**
 * @brief Compute Alfven speed
 */
__attribute__((always_inline)) INLINE static float
mhd_get_comoving_Alfven_speed(const struct part *p, const float mu_0) {

  /* Recover some data */
  const float rho = p->rho;
  const float B[3] = {p->mhd.B_over_rho[0] * rho, p->mhd.B_over_rho[1] * rho,
                      p->mhd.B_over_rho[2] * rho};

  /* B squared */
  const float B2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];

  /* Square of Alfven speed */
  const float vA2 = B2 / (mu_0 * rho);

  return sqrtf(vA2);
}

/**
 * @brief Compute magnetosonic speed
 */
__attribute__((always_inline)) INLINE static float
mhd_get_comoving_magnetosonic_speed(const struct part *p, const float mu_0) {

  /* Compute fast magnetosonic speed, the Pythagorean addition of the sound
   * speed and Alfven speed */
  const float cs = hydro_get_comoving_soundspeed(p);
  const float cs2 = cs * cs;

  const float vA = mhd_get_comoving_Alfven_speed(p, mu_0);
  const float vA2 = vA * vA;

  const float cms2 = cs2 + vA2;

  return sqrtf(cms2);
}

/**
 * @brief Compute fast magnetosonic wave phase veolcity
 */
__attribute__((always_inline)) INLINE static float
mhd_get_comoving_fast_magnetosonic_wave_phase_velocity(const float dx[3],
                                                       const struct part *p,
                                                       const float a,
                                                       const float mu_0) {

  /* Get r and 1/r. */
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float rho = p->rho;
  const float B[3] = {
      p->mhd.B_over_rho[0] * rho,
      p->mhd.B_over_rho[1] * rho,
      p->mhd.B_over_rho[2] * rho,
  };

  /* B dot r. */
  const float Br = B[0] * dx[0] + B[1] * dx[1] + B[2] * dx[2];
  const float permeability_inv = 1.0f / mu_0;

  /* Compute effective sound speeds */
  const float cs = p->force.soundspeed;
  const float cs2 = cs * cs;
  const float c_ms = mhd_get_comoving_magnetosonic_speed(p, mu_0);
  const float c_ms2 = c_ms * c_ms;
  const float projection_correction = c_ms2 * c_ms2 - 4.0f * permeability_inv *
                                                          cs2 * Br * r_inv *
                                                          Br * r_inv / rho;

  const float v_fmsw2 = 0.5f * (c_ms2 + sqrtf(projection_correction));

  return sqrtf(v_fmsw2);
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
    const float dx[3], const struct part *pi, const struct part *pj,
    const float mu_ij, const float beta, const float mu_0) {

  /* Compute pairwise MHD signal velocity,
   * the sum of two particles' fast magnetosonic speeds
   * (i.e. the Pythagorean addition of their sound speed and Alfven speed),
   * and a von Neumann-type correction. */

  const float v_sigi = mhd_get_comoving_magnetosonic_speed(pi, mu_0);
  const float v_sigj = mhd_get_comoving_magnetosonic_speed(pj, mu_0);

  return v_sigi + v_sigj - beta * mu_ij;
}

/**
 * @brief returns the signal velocity
 *
 * @brief p  the particle
 */
__attribute__((always_inline)) INLINE static float hydro_get_signal_velocity(
    const struct part *p) {

  return 0.;
}
/**
 * @brief returns the div_v
 *
 * @brief p  the particle
 */
__attribute__((always_inline)) INLINE static float hydro_get_div_v(
    const struct part *p) {

  return 0.;
}

/**
 * @brief Does some extra hydro operations once the actual physical time step
 * for the particle is known.
 *
 * @param p The particle to act upon.
 * @param dt Physical time step of the particle during the next step.
 */
__attribute__((always_inline)) INLINE static void hydro_timestep_extra(
    struct part *p, const float dt) {}

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
mhd_set_v_sig_based_on_velocity_kick(struct part *p,
                                     const struct cosmology *cosmo,
                                     const float dv_phys, 
				     const float mu_0) {

  /* Compute the velocity kick in comoving coordinates */
  const float dv = dv_phys / cosmo->a_factor_sound_speed;

  /* Fast magnetosonic speed */
  const float cms = mhd_get_comoving_magnetosonic_speed(p,mu_0);

  /* Update the signal velocity */
  p->mhd.v_sig =
      fmaxf(2.f * cms, p->mhd.v_sig + const_viscosity_beta * dv);
}

/**
 * @brief Calculate time derivative of dedner scalar
 *
 * @param p The particle to act upon
 * @param a The current value of the cosmological scale factor
 */
__attribute__((always_inline)) INLINE static float mhd_get_psi_over_ch_dt(
    struct part *p, const float a, const float a_factor_sound_speed,
    const float H, const struct hydro_props *hydro_props, const float mu_0) {

  /* Retrieve inverse of smoothing length. */
  const float h = p->h;
  const float h_inv = 1.0f / h;

  /* Compute Dedner cleaning speed. */
  const float ch = 0.5f * p->mhd.v_sig;

  /* Compute Dedner cleaning scalar time derivative. */
  const float hyp = hydro_props->mhd.hyp_dedner;
  const float hyp_divv = hydro_props->mhd.hyp_dedner_divv;
  const float par = hydro_props->mhd.par_dedner;

  const float divV = p->force.gradient_vx[0] + p->force.gradient_vy[1] + p->force.gradient_vz[2];
  const float div_B = p->mhd.divB;
  const float div_v = a * a * divV - 3.0f * a * a * H;
  const float psi_over_ch = p->mhd.psi_over_ch;

  const float cp = ch;
  const float tau_inv = par * cp * h_inv;

  const float hyperbolic_term =
      -hyp * a * a * a_factor_sound_speed * a_factor_sound_speed * ch * div_B;
  const float hyperbolic_divv_term = -hyp_divv * psi_over_ch * div_v;
  const float parabolic_term =
      -a * a_factor_sound_speed * psi_over_ch * tau_inv;
  const float Hubble_term = a * a * H * psi_over_ch;

  return hyperbolic_term + hyperbolic_divv_term + parabolic_term + Hubble_term;
}


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
    struct part *p, const struct hydro_space *hs) {

  p->density.wcount = 0.f;
  p->density.wcount_dh = 0.f;
  p->rho = 0.f;
  p->density.rho_dh = 0.f;
  p->force.f = 0.f;
}

/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 *
 * Also adds/multiplies the cosmological terms if need be.
 *
 * @param p The particle to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part *p, const struct cosmology *cosmo) {

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
}

/**
 * @brief Prepare a particle for the gradient calculation.
 *
 * This function is called after the density loop and before the gradient loop.
 * Nothing to do in this scheme as there are no terms needing a preparation
 * after the density loop.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 * @param hydro_props Hydrodynamic properties.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_gradient(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor, const float mu_0) {
  /* MHD variables ----------------------------------------------*/
  p->mhd.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);
  /* END MHD variables ------------------------------------------*/
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
    struct part *p) {

  zero_sym_matrix(&p->gradient.c_matrix_inv);
  for (int i = 0; i < 3; ++i) p->gradient.gradient_vx[i] = 0.f;
  for (int i = 0; i < 3; ++i) p->gradient.gradient_vy[i] = 0.f;
  for (int i = 0; i < 3; ++i) p->gradient.gradient_vz[i] = 0.f;
  for (int i = 0; i < 3; ++i) p->gradient.gradient_u[i] = 0.f;

  /* MHD variables ----------------------------------------------*/
  /* Zero the fields updated by the mhd gradient loop */
  p->mhd.curl_B[0] = 0.0f;
  p->mhd.curl_B[1] = 0.0f;
  p->mhd.curl_B[2] = 0.0f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      p->mhd.grad_B_tensor[i][j] = 0.0f;
    }
  }

  p->mhd.plasma_beta_rms = 0.0f;
  p->mhd.neighbour_number = 0.0f;

  /* Initialise MHD signal velocity */
  //const float cms = mhd_get_comoving_magnetosonic_speed(p, mu_0);
  //p->mhd.v_sig = 2.0f * cms;

  /* END MHD variables -------------------------------------------*/
}

/**
 * @brief Finishes the gradient calculation.
 *
 * Multiplies the C-matrix by the appropiate constants.
 *
 * Also adds/multiplies the cosmological terms if need be.
 * Nothing to do in this scheme as the gradient loop is not used.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_end_gradient(
    struct part *p, const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor, const float mu_0) {

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Finish the construction of the inverse of the c-matrix by
   * multiplying in the factors of h coming from W */
  sym_matrix_multiply_by_scalar(&p->gradient.c_matrix_inv, h_inv_dim);

  /* Finish the construction of the inverse of the velocity gradient
   * multiplying in the factors of h coming from W */
  for (int i = 0; i < 3; ++i) p->gradient.gradient_vx[i] *= h_inv_dim;
  for (int i = 0; i < 3; ++i) p->gradient.gradient_vy[i] *= h_inv_dim;
  for (int i = 0; i < 3; ++i) p->gradient.gradient_vz[i] *= h_inv_dim;

  /* Finish the construction of the inverse of the internal energy gradient
   * multiplying in the factors of h coming from W */
  for (int i = 0; i < 3; ++i) p->gradient.gradient_u[i] *= h_inv_dim;
  
  /* MHD variables ----------------------------------------------*/
  /* Recover some data */
  const float rho = p->rho;
  const float P = p->force.pressure;

  float B[3];
  for (int k = 0; k < 3; k++) {
    B[k] = p->mhd.B_over_rho[k] * rho;
  }

  const float B2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];

  /* Finalise local plasma beta mean square calculation */
  const float Pmag_inv = B2 ? 2.0f * mu_0 / B2 : FLT_MAX;
  const float plasma_beta = P * Pmag_inv;

  p->mhd.neighbour_number += 1.0f;
  p->mhd.plasma_beta_rms += plasma_beta * plasma_beta;

  p->mhd.plasma_beta_rms = sqrtf(
      p->mhd.plasma_beta_rms /
      p->mhd.neighbour_number); /* Divisor guaranteed to be strictly positive */
  
  for (int i = 0; i < 3; ++i) 
    for (int j = 0; j < 3; ++j) 
      p->mhd.grad_B_tensor[i][j] *= h_inv_dim;

  /* END MHD variables -------------------------------------------*/
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
    struct part *p, struct xpart *xp, const struct cosmology *cosmo) {

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
  
  /* MHD variables ----------------------------------------------*/
  
  p->mhd.divB = 0.0f;
  p->mhd.curl_B[0] = 0.0f;
  p->mhd.curl_B[1] = 0.0f;
  p->mhd.curl_B[2] = 0.0f;
  
  /* END MHD variables -------------------------------------------*/
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
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor, const float dt_alpha,
    const float dt_therm) {

  /* Compute the pressure */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);

  /* Compute the sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Invert the c-matrix */
  sym_matrix_invert(&p->force.c_matrix, &p->gradient.c_matrix_inv);

  /* Finish computation of velocity gradient (eq. 18) */
  sym_matrix_multiply_by_vector(p->force.gradient_vx, &p->force.c_matrix,
                                p->gradient.gradient_vx);
  sym_matrix_multiply_by_vector(p->force.gradient_vy, &p->force.c_matrix,
                                p->gradient.gradient_vy);
  sym_matrix_multiply_by_vector(p->force.gradient_vz, &p->force.c_matrix,
                                p->gradient.gradient_vz);

  /* Finish computation of u gradient (same as eq. 18) */
  sym_matrix_multiply_by_vector(p->force.gradient_u, &p->force.c_matrix,
                                p->gradient.gradient_u);

  /* Update other variables. */
  p->force.pressure = pressure;
  p->force.soundspeed = soundspeed;
  
  /* MHD variables ----------------------------------------------*/

  float B[3];
  B[0] = p->mhd.B_over_rho[0] * p->rho;
  B[1] = p->mhd.B_over_rho[1] * p->rho;
  B[2] = p->mhd.B_over_rho[2] * p->rho;

  const float B2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
  const float normB = sqrtf(B2);

  float grad_B_mean_square = 0.0f;
  
  /* Finish computation of Grad_B gradient (eq. 18) */
  float tmp_Grad_in[3], tmp_Grad_out[3];
  for (int i = 0; i < 3; i++){
     for (int j = 0; j < 3; j++) 
        tmp_Grad_in[j] = p->mhd.grad_B_tensor[i][j];
     sym_matrix_multiply_by_vector(tmp_Grad_out, &p->force.c_matrix,
                                tmp_Grad_in);
     for (int j = 0; j < 3; j++) 
	p->mhd.grad_B_tensor[i][j] = tmp_Grad_out[j] ;
  }
  p->mhd.divB = p->mhd.grad_B_tensor[0][0] + p->mhd.grad_B_tensor[1][1] + p->mhd.grad_B_tensor[2][2];
  p->mhd.curl_B[0] = p->mhd.grad_B_tensor[2][1] - p->mhd.grad_B_tensor[1][2];
  p->mhd.curl_B[1] = p->mhd.grad_B_tensor[0][2] - p->mhd.grad_B_tensor[2][0];
  p->mhd.curl_B[2] = p->mhd.grad_B_tensor[1][0] - p->mhd.grad_B_tensor[0][1];
  
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      grad_B_mean_square +=
          p->mhd.grad_B_tensor[i][j] * p->mhd.grad_B_tensor[i][j];
    }
  }

  const float alpha_AR_max = p->mhd.art_diff_beta;

  p->mhd.alpha_AR =
      normB ? fminf(alpha_AR_max, p->h * sqrtf(grad_B_mean_square) / normB) : 0.0f;

  /* END MHD variables -------------------------------------------*/
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
    struct part *p) {

  /* Reset the acceleration. */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* Reset the time derivatives. */
  p->u_dt = 0.0f;
  p->force.h_dt = 0.0f;
  p->force.mu_tilde = 0.0f;
  
  /* MHD variables ----------------------------------------------*/
  
  /* Zero the fields updated by the mhd force loop */

  p->mhd.B_over_rho_dt[0] = 0.0f;
  p->mhd.B_over_rho_dt[1] = 0.0f;
  p->mhd.B_over_rho_dt[2] = 0.0f;

  p->mhd.B_over_rho_dt_AR[0] = 0.0f;
  p->mhd.B_over_rho_dt_AR[1] = 0.0f;
  p->mhd.B_over_rho_dt_AR[2] = 0.0f;

  p->mhd.u_dt_AR = 0.0f;

  /* END MHD variables -------------------------------------------*/
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 * @param cosmo The cosmological model
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values(
    struct part *p, const struct xpart *xp, const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor, const float mu_0) {

  /* Re-set the predicted velocities */
  p->v[0] = xp->v_full[0];
  p->v[1] = xp->v_full[1];
  p->v[2] = xp->v_full[2];

  /* Re-set the entropy */
  p->u = xp->u_full;

  /* Re-compute the pressure */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Update variables */
  p->force.pressure = pressure;
  p->force.soundspeed = soundspeed;

  /* MHD variables ----------------------------------------------*/

  /* Re-set the predicted magnetic flux densities */
  p->mhd.B_over_rho[0] = xp->mhd.B_over_rho_full[0];
  p->mhd.B_over_rho[1] = xp->mhd.B_over_rho_full[1];
  p->mhd.B_over_rho[2] = xp->mhd.B_over_rho_full[2];

  p->mhd.psi_over_ch = xp->mhd.psi_over_ch_full;
  
  p->mhd.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);
  /* Re-set MHD signal velocity */
  const float cms = mhd_get_comoving_magnetosonic_speed(p, mu_0);
  p->mhd.v_sig = fmaxf(p->mhd.v_sig, 2.0f * cms);
  
  /* END MHD variables -------------------------------------------*/
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
    struct part *p, const struct xpart *xp, float dt_drift, float dt_therm,
    float dt_kick_grav, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct pressure_floor_props *pressure_floor, const float mu_0) {

  /* Predict the internal energy */
  p->u += p->u_dt * dt_therm;

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
  const float floor_u = gas_internal_energy_from_entropy(p->rho, floor_A);

  /* Check against absolute minimum */
  const float min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  p->u = max(p->u, floor_u);
  p->u = max(p->u, min_u);

  /* Compute the new pressure */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  p->force.pressure = pressure;
  p->force.soundspeed = soundspeed;

  /* MHD variables ----------------------------------------------*/

  /* Predict the magnetic flux density */
  p->mhd.B_over_rho[0] += p->mhd.B_over_rho_dt[0] * dt_therm;
  p->mhd.B_over_rho[1] += p->mhd.B_over_rho_dt[1] * dt_therm;
  p->mhd.B_over_rho[2] += p->mhd.B_over_rho_dt[2] * dt_therm;

  p->mhd.psi_over_ch += p->mhd.psi_over_ch_dt * dt_therm;
  
  p->mhd.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);
  /* Initialise MHD signal velocity */
  const float cms = mhd_get_comoving_magnetosonic_speed(p, mu_0);
  p->mhd.v_sig = fmaxf(p->mhd.v_sig, 2.0f * cms);
  
  /* END MHD variables -------------------------------------------*/
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
    struct part *p, const struct cosmology *cosmo, const struct hydro_props * hydro_props, 
    const float mu_0) {

  p->force.h_dt *= p->h * hydro_dimension_inv;
  
  /* MHD variables ----------------------------------------------*/
  
  /* Get time derivative of Dedner scalar */
  p->mhd.psi_over_ch_dt = mhd_get_psi_over_ch_dt(
      p, cosmo->a, cosmo->a_factor_sound_speed, cosmo->H, hydro_props, mu_0);


  /* Hubble expansion contribution to induction equation */
  const float Hubble_induction_pref =
      cosmo->a * cosmo->a * cosmo->H * (1.5f * hydro_gamma - 2.f);
  p->mhd.B_over_rho_dt[0] +=
      Hubble_induction_pref * p->mhd.B_over_rho[0];
  p->mhd.B_over_rho_dt[1] +=
      Hubble_induction_pref * p->mhd.B_over_rho[1];
  p->mhd.B_over_rho_dt[2] +=
      Hubble_induction_pref * p->mhd.B_over_rho[2];

  /*const float divV = 
  p->force.gradient_vx[0]+
  p->force.gradient_vy[1]+
  p->force.gradient_vz[2];

  p->mhd.B_over_rho_dt[0] -=
      divV * p->mhd.B_over_rho[0];
  p->mhd.B_over_rho_dt[1] -=
      divV * p->mhd.B_over_rho[1];
  p->mhd.B_over_rho_dt[2] -=
      divV * p->mhd.B_over_rho[2];*/
 /*
  const float plasma_beta_i = p->mhd.plasma_beta_rms;
  const float scale_i = 0.125f * (10.0f - plasma_beta_i);
  const float tensile_correction_scale_i = fmaxf(0.0f, fminf(scale_i, 1.0f));

  for (int i = 0; i < 3; i++){
     p->a_hydro[i] += 1.f / (mu_0) * ( 
        p->mhd.B_over_rho[0] * p->mhd.grad_B_tensor[i][0] +
        p->mhd.B_over_rho[1] * p->mhd.grad_B_tensor[i][1] +
        p->mhd.B_over_rho[2] * p->mhd.grad_B_tensor[i][2]
        - tensile_correction_scale_i * p->mhd.divB * p->mhd.B_over_rho[i]);
 }*/ 
  /* END MHD variables -------------------------------------------*/
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
 * @param hydro_props The constants used in the scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part *p, struct xpart *xp, float dt_therm, float dt_grav,
    float dt_grav_mesh, float dt_hydro, float dt_kick_corr,
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
  
  /* MHD variables ----------------------------------------------*/
  
  /* Integrate the magnetic flux density forward in time */
  const float delta_Bx = p->mhd.B_over_rho_dt[0] * dt_therm;
  const float delta_By = p->mhd.B_over_rho_dt[1] * dt_therm;
  const float delta_Bz = p->mhd.B_over_rho_dt[2] * dt_therm;

  /* Do not decrease the magnetic flux density by more than a factor of 2*/
  xp->mhd.B_over_rho_full[0] = xp->mhd.B_over_rho_full[0] + delta_Bx;
  xp->mhd.B_over_rho_full[1] = xp->mhd.B_over_rho_full[1] + delta_By;
  xp->mhd.B_over_rho_full[2] = xp->mhd.B_over_rho_full[2] + delta_Bz;

  xp->mhd.psi_over_ch_full += p->mhd.psi_over_ch_dt * dt_therm;
  /* END MHD variables -------------------------------------------*/
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
 * @param mu_0 Vacuum permeability constant
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor, const float mu_0) {

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

  /* Compute the pressure */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);

  /* Compute the sound speed */
  const float soundspeed = gas_soundspeed_from_internal_energy(p->rho, p->u);

  p->force.pressure = pressure;
  p->force.soundspeed = soundspeed;

  /* MHD variables ----------------------------------------------*/
  
  /* Set Restitivity Eta */
  p->mhd.resistive_eta = hydro_props->mhd.mhd_eta;
  /* Set Monopole subtraction factor */
  p->mhd.monopole_beta = hydro_props->mhd.monopole_subtraction;
  /* Set Artificial Difussion */
  p->mhd.art_diff_beta = hydro_props->mhd.art_diffusion;

  /* Convert B into B/rho */
  p->mhd.B_over_rho[0] /= p->rho;
  p->mhd.B_over_rho[1] /= p->rho;
  p->mhd.B_over_rho[2] /= p->rho;

  /* Convert to co-moving B/rho */
  p->mhd.B_over_rho[0] *= powf(cosmo->a, 1.5f * hydro_gamma);
  p->mhd.B_over_rho[1] *= powf(cosmo->a, 1.5f * hydro_gamma);
  p->mhd.B_over_rho[2] *= powf(cosmo->a, 1.5f * hydro_gamma);

  /* Instantiate full step magnetic field */
  xp->mhd.B_over_rho_full[0] = p->mhd.B_over_rho[0];
  xp->mhd.B_over_rho_full[1] = p->mhd.B_over_rho[1];
  xp->mhd.B_over_rho_full[2] = p->mhd.B_over_rho[2];

  /* Instantiate full step magnetic Dedner scalar */
  xp->mhd.psi_over_ch_full = p->mhd.psi_over_ch;

  /* Instantiate Alfven speed */
  p->mhd.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);
  
  /* END MHD variables -------------------------------------------*/
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
    struct part *p, struct xpart *xp) {

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
    const struct part *p, const struct xpart *xp, const double time) {}

#endif /* SWIFT_MAGMA_HYDRO_H */
