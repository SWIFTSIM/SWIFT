/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_STRENGTH_DEFAULT_H
#define SWIFT_PLANETARY_STRENGTH_DEFAULT_H

/**
 * @file Planetary/strength/default/hydro_strength.h
 * @brief Planetary implementation of SPH with default material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength.h"

__attribute__((always_inline)) INLINE static float hydro_compute_timestep_strength(
    const struct part *restrict p, const struct hydro_props *restrict hydro_properties, 
    const float dt_cfl) {
    
  const float elastic_timestep_factor = hydro_properties->CFL_condition; // ### Set as same as CFL factor for now. Treat this similarly to CFL
  const float norm_dS_dt = norm_sym_matrix(&p->strength_data.dS_dt);
  const float shear_mod = material_shear_mod(p->mat_id);

  float dt_elastic;
    
  if (norm_dS_dt * dt_cfl > elastic_timestep_factor * shear_mod) {
    dt_elastic = elastic_timestep_factor * shear_mod / norm_dS_dt;
  } else {
    dt_elastic = FLT_MAX;
  }

  float dt_strength = dt_elastic;
    
#if defined(STRENGTH_DAMAGE)
  // Update dt_strength with damage contribution
  float dt_damage;
  damage_timestep(p, dt_cfl, &dt_damage);
  dt_strength = fminf(dt_elastic, dt_damage);
#endif /* STRENGTH_DAMAGE */
    
  return dt_strength;
}

/**
 * @brief Prepares extra strength parameters for a particle for the density
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_init_part_extra_strength(struct part *restrict p) {}

/**
 * @brief Finishes extra strength parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_density_extra_strength(struct part *restrict p) {}

/**
 * @brief Prepares extra strength parameters for a particle for the force
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_force_extra_strength(struct part *restrict p, struct xpart *restrict xp,
                                   const float density, const float u) {

  // Set the density to be used in the force loop to be the evolved density
  p->rho = p->strength_data.rho_evol;

#ifdef PLANETARY_FIXED_ENTROPY
  /* Override the internal energy to satisfy the fixed entropy. 
    * Needs to happen here as well because of the updated rho. */
  p->u = gas_internal_energy_from_entropy(p->rho, p->s_fixed, p->mat_id);
  xp->u_full = p->u;
#endif
    
  const float pressure =
      gas_pressure_from_internal_energy(p->rho, p->u, p->mat_id);

  hydro_set_stress_tensor(p, pressure);
}

/**
 * @brief Resets strength time derivative fields in preparation
  * for the sums taking  place in the various force tasks.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_reset_acceleration_strength(struct part *restrict p) {

  p->strength_data.drho_dt = 0.0f;
    
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      p->strength_data.dv_force_loop[i][j] = 0.f;
    }
  }

  zero_sym_matrix(&p->strength_data.dS_dt);
}

/**
 * @brief Finishes extra strength parts of the force calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_force_extra_strength(struct part *restrict p) {

 calculate_dS_dt(p);

  // Calculate dD/dt for timestep
#if defined(STRENGTH_DAMAGE)
  const int phase_state = p->phase_state;
  const float density = p->strength_data.rho_evol;
  const float u = p->u;
  const float damage = p->strength_data.damage;
  const float yield_stress = compute_yield_stress_damaged(p, phase_state, density, u, damage);

  p->strength_data.dD_dt = 0.f;
    
  float tensile_cbrtD_dt = 0.f;
  float number_of_activated_flaws = 0;
  calculate_tensile_cbrtD_dt(p, &tensile_cbrtD_dt, &number_of_activated_flaws, 
                                        p->strength_data.deviatoric_stress_tensor, damage, density, u);
  if (p->strength_data.tensile_damage < number_of_activated_flaws / (float)p->strength_data.number_of_flaws) {
    // Chain rule d(D^(1/3))/dt  = d(D^(1/3))/dD * dD/dt
    p->strength_data.dD_dt += 3.f * powf(p->strength_data.tensile_damage, 2.f / 3.f) * tensile_cbrtD_dt;
  }
    
  float shear_dD_dt = 0.f;
  calculate_shear_dD_dt(p, &shear_dD_dt, p->strength_data.deviatoric_stress_tensor, yield_stress, density, u);

  if (p->strength_data.shear_damage < 1.f) {
    p->strength_data.dD_dt += shear_dD_dt;
  } 

#endif /* STRENGTH_DAMAGE */
}

/**
 * @brief Sets the values of additional particle strength properties at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values_extra_strength(
    struct part *restrict p, const struct xpart *restrict xp) {

  p->rho = xp->strength_data.rho_evol_full;
  p->strength_data.rho_evol = xp->strength_data.rho_evol_full;
  p->strength_data.deviatoric_stress_tensor = xp->strength_data.deviatoric_stress_tensor_full;
  #if defined(STRENGTH_DAMAGE)
    p->strength_data.damage = xp->strength_data.damage_full;
    p->strength_data.tensile_damage = xp->strength_data.tensile_damage_full;
    p->strength_data.shear_damage = xp->strength_data.shear_damage_full;
  #endif /* STRENGTH_DAMAGE */
}

/**
 * @brief Predict additional particle strength properties forward in time when
 * drifting. At beginning of hydro function, before hydro quantities have been drifted.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra_strength_beginning(
    struct part *restrict p, const float dt_therm) {

  // ### FOR LEAPFROG: dS_dt is calculated similarly to e.g. du_dt and S is updated similarly to u
  // ### dDamage_dt does not depend on positions or velocities of particle neighbours, only on properties
  // ### of the particle itself. Therefore dDamage_dt is recalculated each time damage is updated.
  // ## Damage depends on S and S can depend on damage through Y makes.

  const int phase_state = p->phase_state;
  const float density = p->strength_data.rho_evol;
  const float u = p->u;

  #if defined(STRENGTH_DAMAGE)
    const float damage = p->strength_data.damage;
    const float yield_stress = compute_yield_stress_damaged(p, phase_state, density, u, damage);

    evolve_damage(p, &p->strength_data.tensile_damage, &p->strength_data.shear_damage, &p->strength_data.damage, p->strength_data.deviatoric_stress_tensor, yield_stress, density, u, dt_therm);
  #else
    const float yield_stress = compute_yield_stress(p, phase_state, density, u);
  #endif /* STRENGTH_DAMAGE */

  evolve_deviatoric_stress(p, &p->strength_data.deviatoric_stress_tensor, phase_state, dt_therm);

  adjust_deviatoric_stress_tensor_by_yield_stress(
        p, &p->strength_data.deviatoric_stress_tensor, yield_stress, density, u);
}

/**
 * @brief Predict additional particle strength properties forward in time when
 * drifting. At end of hydro function, after hydro quantities have been drifted.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra_strength_end(
    struct part *restrict p, const float dt_therm) {

  p->strength_data.rho_evol += p->strength_data.drho_dt * dt_therm;

  /* compute minimum density */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  const float min_rho = p->mass * kernel_root * h_inv_dim;

  // Overwrite stored hydro qunatities with those calculated based on evolved density.
  p->strength_data.rho_evol = max(p->strength_data.rho_evol, min_rho);
  p->rho = p->strength_data.rho_evol;

  const float density = p->strength_data.rho_evol;
  const float u = p->u;
  const float pressure =
      gas_pressure_from_internal_energy(density, u, p->mat_id);  
  const float soundspeed =
      gas_soundspeed_from_internal_energy(p->rho, p->u, p->mat_id);

  p->force.pressure = pressure;
  p->force.soundspeed = soundspeed;

  p->force.v_sig = max(p->force.v_sig, 2.f * soundspeed);

  p->phase_state =
    (enum mat_phase_state)material_phase_state_from_internal_energy(
     p->rho, p->u, p->mat_id);
    
  hydro_set_stress_tensor(p, pressure);
}

/**
 * @brief Kick the additional particle strength properties
 * At beginning of hydro function, before hydro quantities have been kicked.
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra_strength_beginning(
    struct part *restrict p, struct xpart *restrict xp, float dt_therm) {

  const int phase_state = xp->phase_state_full;
  const float density = xp->strength_data.rho_evol_full;
  const float u = xp->u_full;

  #if defined(STRENGTH_DAMAGE)
    const float damage = xp->strength_data.damage_full;
    const float yield_stress = compute_yield_stress_damaged(p, phase_state, density, u, damage);

    evolve_damage(p, &xp->strength_data.tensile_damage_full, &xp->strength_data.shear_damage_full, &xp->strength_data.damage_full, xp->strength_data.deviatoric_stress_tensor_full, yield_stress, density, u,  dt_therm);
  #else
      const float yield_stress = compute_yield_stress(p, phase_state, density, u);
  #endif /* STRENGTH_DAMAGE */

  evolve_deviatoric_stress(p, &xp->strength_data.deviatoric_stress_tensor_full, phase_state, dt_therm);

  adjust_deviatoric_stress_tensor_by_yield_stress(
        p, &xp->strength_data.deviatoric_stress_tensor_full, yield_stress, density, u);

  const float delta_rho = p->strength_data.drho_dt * dt_therm;

  xp->strength_data.rho_evol_full =
      max(xp->strength_data.rho_evol_full + delta_rho, 0.5f * xp->strength_data.rho_evol_full);

  /* Minimum SPH quantities */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  const float floor_rho = p->mass * kernel_root * h_inv_dim;
  if (xp->strength_data.rho_evol_full < floor_rho) {
    xp->strength_data.rho_evol_full = floor_rho;
    p->strength_data.drho_dt = 0.f;
  }

  // Overwrite based on evolved density
  xp->phase_state_full =
    (enum mat_phase_state)material_phase_state_from_internal_energy(
     xp->strength_data.rho_evol_full, xp->u_full, p->mat_id);
}

/**
 * @brief Kick the additional particle strength properties
 * At end of hydro function, after hydro quantities have been kicked.
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra_strength_end(
    struct part *restrict p, struct xpart *restrict xp, float dt_therm) {}

/**
 * @brief Initialises the strength properties for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_first_init_part_strength(
    struct part *restrict p, struct xpart *restrict xp) {
    
  p->strength_data.rho_evol = p->rho;
  xp->strength_data.rho_evol_full = p->strength_data.rho_evol;

  for (int i = 0; i < 6; i++) {
    p->strength_data.deviatoric_stress_tensor.elements[i] = 0.f;
    xp->strength_data.deviatoric_stress_tensor_full.elements[i] = 0.f;
  }
  #ifdef STRENGTH_DAMAGE
    p->damage = 0.f;
    p->tensile_damage = 0.f;
    p->shear_damage = 0.f;

    xp->damage_full = 0.f;
    xp->tensile_damage_full = 0.f;
    xp->shear_damage_full = 0.f;
 #endif /* STRENGTH_DAMAGE */
}

#endif /* SWIFT_PLANETARY_STRENGTH_DEFAULT_H */
