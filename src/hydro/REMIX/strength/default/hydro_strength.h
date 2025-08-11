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
#ifndef SWIFT_REMIX_STRENGTH_DEFAULT_H
#define SWIFT_REMIX_STRENGTH_DEFAULT_H

/**
 * @file REMIX/strength/default/hydro_strength.h
 * @brief REMIX implementation of SPH with default material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength.h"

__attribute__((always_inline)) INLINE static void hydro_compute_timestep_strength(
    const struct part *restrict p, const struct hydro_props *restrict hydro_properties, 
    float dt_cfl) {
    
  const float elastic_timestep_factor = hydro_properties->CFL_condition; // ### Set as same as CFL factor for now. Treat this similarly to CFL
  const float norm_dS_dt = norm_sym_matrix(&p->strength_data.dS_dt);
  const float shear_mod = material_shear_mod(p->mat_id);
    
  if (norm_dS_dt * dt_cfl > elastic_timestep_factor * shear_mod) {
    dt_cfl = elastic_timestep_factor * shear_mod / norm_dS_dt;
  }
    
  // Update dt_strength with damage contribution
  damage_timestep(p, dt_cfl);
}

__attribute__((always_inline)) INLINE static void
hydro_compute_max_wave_speed_strength(const struct part *restrict p, const float soundspeed, const float density, float *wave_speed) {
  if (p->phase_state == mat_phase_state_solid) {
    const float shear_mod = material_shear_mod(p->mat_id);
    
    // Speed of longitudinal elastic wave
    *wave_speed = sqrtf(soundspeed * soundspeed + (4.f / 3.f) * shear_mod / density);
  }
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
 * @brief Prepares extra strength parameters for a particle for the gradient
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_extra_strength(struct part *restrict p) {}

/**
 * @brief Finishes extra strength parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_extra_strength(struct part *restrict p) {}

/**
 * @brief Prepares extra strength parameters for a particle for the force
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_force_extra_strength(struct part *restrict p,
                                   const float density, const float u) {

  const float pressure =
      gas_pressure_from_internal_energy(density, u, p->mat_id);

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

  // Update p->strength_data.dD/dt for timestep
  const int phase_state = p->phase_state;
  const float density = p->rho_evol;
  const float u = p->u;
  update_dD_dt(p, phase_state, density, u);
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

  p->strength_data.deviatoric_stress_tensor = xp->strength_data.deviatoric_stress_tensor_full;

  hydro_reset_predicted_values_extra_damage(p, xp);
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
  const float density = p->rho_evol;
  const float u = p->u;

  const float yield_stress = compute_yield_stress_damage(p, phase_state, density, u);
  hydro_predict_evolve_damage(p, yield_stress, density, u, dt_therm);

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

  const float density = p->rho_evol;
  const float u = p->u;
  const float pressure =
      gas_pressure_from_internal_energy(density, u, p->mat_id);  
    
  hydro_set_stress_tensor(p, pressure);
}

/**
 * @brief Kick the additional particle strength properties.
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
  const float density = xp->rho_evol_full;
  const float u = xp->u_full;


  const float yield_stress = compute_yield_stress_damage_full(p, xp, phase_state, density, u);
  hydro_kick_evolve_damage(p, xp, yield_stress, density, u,  dt_therm);

  evolve_deviatoric_stress(p, &xp->strength_data.deviatoric_stress_tensor_full, phase_state, dt_therm);

  adjust_deviatoric_stress_tensor_by_yield_stress(
        p, &xp->strength_data.deviatoric_stress_tensor_full, yield_stress, density, u);
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

  for (int i = 0; i < 6; i++) {
    p->strength_data.deviatoric_stress_tensor.elements[i] = 0.f;
    xp->strength_data.deviatoric_stress_tensor_full.elements[i] = 0.f;
  }
  
  hydro_first_init_part_damage(p, xp);
}

#endif /* SWIFT_REMIX_STRENGTH_DEFAULT_H */
