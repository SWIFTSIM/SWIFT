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
#ifndef SWIFT_MINIMAL_ENGINEERING_H
#define SWIFT_MINIMAL_ENGINEERING_H

/**
 * @file Minimal_Engineering/hydro.h
 * @brief Minimal weakly compressible implementation of SPH 
 *
 * corresponds to Morris, Fox and Zhu,  J. Comp. Physics, 136, 214 (1997)
 */

#include "adiabatic_index.h" /* RGB: this function needs to be updated for gamma=1, as used by Morris */
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
 * *** not used in weakly compressible flow ***
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy(const struct part *restrict p,
                                   const struct xpart *restrict xp) {
  error("internal energy is not defined in weakly copressible EOS");
  return 0.f;
}

/**
 * @brief Returns the physical internal energy of a particle at the last
 * time the particle was kicked.
 * *** not used in weakly compressible flow ***
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy(const struct part *restrict p,
                                   const struct xpart *restrict xp,
                                   const struct cosmology *cosmo) {

  error("internal energy is not defined in weakly copressible EOS");
  return 0.f;
}

/**
 * @brief Returns the comoving internal energy of a particle drifted to the
 * current time.  *** not used ***
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_internal_energy(const struct part *restrict p) {
  error("internal energy is not defined in weakly copressible EOS");
  return 0.f;
}

/**
 * @brief Returns the physical internal energy of a particle drifted to the
 * current time. *** not used ***
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_physical_internal_energy(const struct part *restrict p,
                                           const struct cosmology *cosmo) {
  error("internal energy is not defined in weakly copressible EOS");
  return 0.f;
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

  return cosmo->a_factor_pressure *
         gas_pressure_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the comoving entropy of a particle at the last
 * time the particle was kicked.  ***  not used ***
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_entropy(
    const struct part *restrict p, const struct xpart *restrict xp) {
  error("entropy is not defined in weakly copressible EOS");
  return 0.f;
}

/**
 * @brief Returns the physical entropy of a particle at the last
 * time the particle was kicked.  *** not used ***
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
  error("entropy is not defined in weakly copressible EOS");
  return 0.f;
}

/**
 * @brief Returns the comoving entropy of a particle drifted to the
 * current time. *** not used ***
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_entropy(const struct part *restrict p) {
  error("entropy is not defined in weakly copressible EOS");
  return 0.f;
}

/**
 * @brief Returns the physical entropy of a particle drifted to the
 * current time. *** not used ***
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_physical_entropy(const struct part *restrict p,
                                   const struct cosmology *cosmo) {

  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  error("entropy is not defined in weakly copressible EOS");
  return 0.f;
}

/**
 * @brief Returns the comoving sound speed of a particle. RGB no need to store sound speed of each particle.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_soundspeed(const struct part *restrict p) {
  const float dummy = 0.f;  /* does not depend on particle properties */
  return gas_soundspeed_from_internal_energy( dummy, dummy );
}

/**
 * @brief Returns the physical sound speed of a particle. RGB no need to store sound speed of each particle.
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_soundspeed(const struct part *restrict p,
                              const struct cosmology *cosmo) {
  const float dummy = 0.f;  /* does not depend on particle properties */
  return cosmo->a_factor_sound_speed * gas_soundspeed_from_internal_energy( dummy, dummy );
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
 * @brief Returns the time derivative of co-moving internal energy of a particle
 *
 * We assume a constant density.  *** not used ***
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy_dt(const struct part *restrict p) {
  error("internal energy is not defined in weakly copressible EOS");
  return 0.f;
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density. *** not used ***
 *
 * @param p The particle of interest
 * @param cosmo Cosmology data structure
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy_dt(const struct part *restrict p,
                                      const struct cosmology *cosmo) {
  error("internal energy is not defined in weakly copressible EOS");
  return 0.f;
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
hydro_set_comoving_internal_energy_dt(struct part *restrict p, float du_dt) {
  error("internal energy is not defined in weakly copressible EOS");
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.  *** not used ***
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param du_dt The new time derivative of the internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy_dt(struct part *restrict p,
                                      const struct cosmology *cosmo,
                                      float du_dt) {
  error("internal energy is not defined in weakly copressible EOS");
}

/**
 * @brief Sets the physical entropy of a particle  *** not used ***
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param entropy The physical entropy
 */
__attribute__((always_inline)) INLINE static void hydro_set_physical_entropy(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const float entropy) {
  error("internal energy is not defined in weakly copressible EOS");
}

/**
 * @brief Sets the physical internal energy of a particle *** not used ***
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
  error("internal energy is not defined in weakly copressible EOS");
}

/**
 * @brief Sets the drifted physical internal energy of a particle.
 *  RGB: Although internal energy is not used, the other quantities need to be set here.
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(struct part *p,
                                           const struct cosmology *cosmo,
                                           const float u) {
  /* internal energy is not used */
  // p->u = u / cosmo->a_factor_internal_energy; 
  
  /* Now recompute the extra quantities *** this is still required *** */

  /* Compute the sound speed */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Update variables. */
  p->force.pressure = pressure;
  // p->force.soundspeed = soundspeed; RGB: this is a EOS constant.

  p->force.v_sig = max(p->force.v_sig, 2.f * soundspeed);  //RGB v_sig may require modifaction later.
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
 * Many quantities are not used in Minimal_Engineering
 *
 * @param p The particle to act upon
 * @param hs #hydro_space containing hydro specific space information. 
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part *restrict p, const struct hydro_space *hs) {

  p->density.wcount = 0.f;       // not required, but neightbour cound may be useful.
  //p->density.wcount_dh = 0.f;  // h is fixed and not adaptive.
  p->rho = 0.f;
  //p->density.rho_dh = 0.f;     // h is fixed and not adaptive.
  p->density.div_v = 0.f;
  //p->density.rot_v[0] = 0.f;   // curl v is not needed.
  //p->density.rot_v[1] = 0.f;
  //p->density.rot_v[2] = 0.f;
}

/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 * Additional quantities such as velocity gradients will also get the final
 * terms added to them here.
 * Many quantities are not used in Minimal_Engineering
 *
 * Also adds/multiplies the cosmological terms if need be.
 *
 * @param p The particle to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part *restrict p, const struct cosmology *cosmo) {

  /*RGB for weakly compressible, just compute density by continuity equation.
    h should not vary, but left as a particle property fro now. */

  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish calculation of the (physical) velocity divergence */
  const float rho_inv = 1.f / p->rho;
  const float a_inv2 = cosmo->a2_inv;
  p->density.div_v *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  p->rho = - p->rho * p->density.div_v ;   //RGB: density updated based on div V [** check this! **]

  // RGB these parts are not needed (keep wcount for debugging for now).
  /* Final operation on the density (add self-contribution). */
  //p->rho += p->mass * kernel_root;
  //p->density.rho_dh -= hydro_dimension * p->mass * kernel_root;
  p->density.wcount += kernel_root;
  //p->density.wcount_dh -= hydro_dimension * kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  //p->rho *= h_inv_dim;
  //p->density.rho_dh *= h_inv_dim_plus_one;
  p->density.wcount *= h_inv_dim;
  //p->density.wcount_dh *= h_inv_dim_plus_one;

  //RGB not required as no Balsara switch is applied.
  /* Finish calculation of the (physical) velocity curl components */
  //p->density.rot_v[0] *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  //p->density.rot_v[1] *= h_inv_dim_plus_one * a_inv2 * rho_inv;
  //p->density.rot_v[2] *= h_inv_dim_plus_one * a_inv2 * rho_inv;

}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * In the desperate case where a particle has no neighbours (likely because
 * of the h_max ceiling), set the particle fields to something sensible to avoid
 * NaNs in the next calculations.
 *  *** this should not happend in weakly compressible fluid. throw an error. ***
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_part_has_no_neighbours(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo) {

  error("a particle has no neighbours! This cannot happen in weakly compressible fluid.");
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
 * *** RGB: many quantitites not need to weakly compressible flow ***
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

  // RGB: no balsara switvch.
  // const float fac_Balsara_eps = cosmo->a_factor_Balsara_eps;

  /* Inverse of the smoothing length */
  const float h_inv = 1.f / p->h;

  /* Compute the norm of the curl */
  // const float curl_v = sqrtf(p->density.rot_v[0] * p->density.rot_v[0] +
  //                           p->density.rot_v[1] * p->density.rot_v[1] +
  //                           p->density.rot_v[2] * p->density.rot_v[2]);

  /* Compute the norm of div v including the Hubble flow term */
  // const float div_physical_v = p->density.div_v + hydro_dimension * cosmo->H;
  // const float abs_div_physical_v = fabsf(div_physical_v);

  /* Compute the pressure */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);

  /* Compute the sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  //RGB not required in Engineering version since we compute density from div v
  /* Compute the "grad h" term */
  //const float rho_inv = 1.f / p->rho;
  //float rho_dh = p->density.rho_dh;
  ///* Ignore changing-kernel effects when h ~= h_max */
  //if (p->h > 0.9999f * hydro_props->h_max) {
  //  rho_dh = 0.f;
  //}
  //const float grad_h_term =
  //    1.f / (1.f + hydro_dimension_inv * p->h * rho_dh * rho_inv);

  /* Compute the Balsara switch */
  /* Pre-multiply in the AV factor; hydro_props are not passed to the iact
   * functions */
  //const float balsara = hydro_props->viscosity.alpha * abs_div_physical_v /
  //                      (abs_div_physical_v + curl_v +
  //                       0.0001f * fac_Balsara_eps * soundspeed * h_inv);

  /* Update variables. */
  //RGB not needed... p->force.f = grad_h_term;
  p->force.f = 0.;
  p->force.pressure = pressure;
  //p->force.soundspeed = soundspeed;
  //p->force.balsara = balsara;
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking  place in the various force tasks.
 * RGB: in weakly compressible limit, some quantities are not required.
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
  //p->u_dt = 0.0f;
  //p->force.h_dt = 0.0f;
  const float dummy = 0.f;
  p->force.v_sig = 2.f * gas_soundspeed_from_internal_energy( dummy, dummy);  // does not depend on properties of particle. Remove particle property later.
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time.  *** RGB: many quantities are not used in weakly compressible flow ***
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 * @param cosmo The cosmological model
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values(
    struct part *restrict p, const struct xpart *restrict xp,
    const struct cosmology *cosmo) {

  /* Re-set the predicted velocities */
  p->v[0] = xp->v_full[0];
  p->v[1] = xp->v_full[1];
  p->v[2] = xp->v_full[2];

  /* Re-set the entropy */
  // p->u = xp->u_full;

  /* Re-compute the pressure */
  // RGB for WC-SPH compute the pressure from the density via EOS. Routine has same name, but internal energy is not used.
  const float dummy = 0.f;
  const float pressure = gas_pressure_from_internal_energy(p->rho, dummy);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  /* Update variables */
  p->force.pressure = pressure;
  //p->force.soundspeed = soundspeed;

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
 * RGB: mnay properties are not required for weakly compressible flow.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
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

  /* Predict the internal energy */
  //p->u += p->u_dt * dt_therm;

  const float h_inv = 1.f / p->h;

  /* Predict smoothing length */
  //const float w1 = p->force.h_dt * h_inv * dt_drift;
  //if (fabsf(w1) < 0.2f) {
  //  p->h *= approx_expf(w1); /* 4th order expansion of exp(w) */
  //} else {
  //  p->h *= expf(w1);
  //}

  /* Predict density */
  // RGB: assume we do not predict the density... Is this appropriate?.
  //const float w2 = -hydro_dimension * w1;
  //if (fabsf(w2) < 0.2f) {
  //  p->rho *= approx_expf(w2); /* 4th order expansion of exp(w) */
  //} else {
  //  p->rho *= expf(w2);
  //}

  /* Check against entropy floor */
  //const float floor_A = entropy_floor(p, cosmo, floor_props);
  //const float floor_u = gas_internal_energy_from_entropy(p->rho, floor_A);

  /* Check against absolute minimum */
  //const float min_u =
  //    hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  //p->u = max(p->u, floor_u);
  //p->u = max(p->u, min_u);

  /* Compute the new pressure */
  const float dummy = 0.f ;
  const float pressure = gas_pressure_from_internal_energy(p->rho, dummy);

  /* Compute the new sound speed */
  const float soundspeed = gas_soundspeed_from_pressure(p->rho, pressure);

  p->force.pressure = pressure;
  //p->force.soundspeed = soundspeed;

  p->force.v_sig = max(p->force.v_sig, 2.f * soundspeed);
}

/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the force and accelerations by the appropiate constants
 * and add the self-contribution term. In most cases, there is little
 * to do here. RGB: nothing for Minimal_Engineering.
 *
 * Cosmological terms are also added/multiplied here.
 *
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_end_force(
    struct part *restrict p, const struct cosmology *cosmo) {

  //p->force.h_dt *= p->h * hydro_dimension_inv;
}

/**
 * @brief Kick the additional variables
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 * RGB: many variables are not used in Minimal_Engineering
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 * @param dt_grav The time-step for this kick (for gravity quantities).
 * @param dt_hydro The time-step for this kick (for hydro quantities).
 * @param dt_kick_corr The time-step for this kick (for gravity corrections).
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part *restrict p, struct xpart *restrict xp, float dt_therm,
    float dt_grav, float dt_hydro, float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Integrate the internal energy forward in time */
  //const float delta_u = p->u_dt * dt_therm;

  /* Do not decrease the energy by more than a factor of 2*/
  //xp->u_full = max(xp->u_full + delta_u, 0.5f * xp->u_full);

  /* Check against entropy floor */
  //const float floor_A = entropy_floor(p, cosmo, floor_props);
  //const float floor_u = gas_internal_energy_from_entropy(p->rho, floor_A);

  /* Check against absolute minimum */
  //const float min_u =
  //    hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  /* Take highest of both limits */
  //const float energy_min = max(min_u, floor_u);

  //if (xp->u_full < energy_min) {
  //  xp->u_full = energy_min;
  //  p->u_dt = 0.f;
  //}
}

/**
 * @brief Converts hydro quantity of a particle at the start of a run
 *
 * This function is called once at the end of the engine_init_particle()
 * routine (at the start of a calculation) after the densities of
 * particles have been computed.
 * This can be used to convert internal energy into entropy for instance.
 * RGB: many properties do not need to be set in Minimal_Engineering.
 *
 * @param p The particle to act upon
 * @param xp The extended particle to act upon
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props) {

  /* Convert the physcial internal energy to the comoving one. */
  /* u' = a^(3(g-1)) u */
  //const float factor = 1.f / cosmo->a_factor_internal_energy;
  //p->u *= factor;
  //xp->u_full = p->u;

  /* Apply the minimal energy limit */
  //const float min_comoving_energy =
  //    hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;
  //if (xp->u_full < min_comoving_energy) {
  //  xp->u_full = min_comoving_energy;
  //  p->u = min_comoving_energy;
  //  p->u_dt = 0.f;
  //}

  /* Compute the pressure */
  const float dummy = 0.f ;
  const float pressure = gas_pressure_from_internal_energy(p->rho, dummy);

  /* Compute the sound speed */
  //const float soundspeed = gas_soundspeed_from_internal_energy(p->rho, p->u);

  p->force.pressure = pressure;
  //p->force.soundspeed = soundspeed;
}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 * RGB: for Minimal_engineering, many properties do not need to be set.
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
  xp->a_grav[0] = 0.f;
  xp->a_grav[1] = 0.f;
  xp->a_grav[2] = 0.f;
  //xp->u_full = p->u;

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
 * ***not used*** for Minimal_Engineering
 *
 * @param p The #part to write to.
 * @param u_init The new initial internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_init_internal_energy(struct part *p, float u_init) {
  error("internal energy is not used in Minimal Engineering hydro");
  //p->u = u_init;
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

#endif /* SWIFT_MINIMAL_ENGINEERING_H */
