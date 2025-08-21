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

/**
 * @brief Updates the hydro (+ strength) time-step based on strength methods.
 *
 *
 * @param dt_cfl The hydro (+ strength) time-step.
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void hydro_compute_timestep_strength(
    float *dt_cfl, const struct part *restrict p, const struct hydro_props *restrict hydro_properties) {

  /* Update dt_cfl with elastic contribution. */
  strength_compute_timestep_stress_tensor(dt_cfl, p, hydro_properties);

  /* Update dt_cfl with damage contribution. */
  strength_compute_timestep_damage(dt_cfl, p);
}

/**
 * @brief Updates the max wave speed based on strength methods.
 *
 * @param wave_speed The wave speed to be updated.
 * @param p The particle of interest.
 * @param soundspeed The sound speed.
 * @param density The sound density.
 */
__attribute__((always_inline)) INLINE static void
hydro_compute_max_wave_speed_strength(float *wave_speed, const struct part *restrict p, const float soundspeed, const float density) {

  strength_compute_max_wave_speed_stress_tensor(wave_speed, p, soundspeed, density);
}

/**
 * @brief Prepares extra strength parameters for a particle for the density
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_init_part_strength(struct part *restrict p) {}

/**
 * @brief Finishes extra strength parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_density_strength(struct part *restrict p) {}

/**
 * @brief Prepares extra strength parameters for a particle for the force
 * calculation.
 *
 * @param p The particle to act upon
 * @param density The density
 * @param u The specific internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_force_strength(struct part *restrict p, struct xpart *restrict xp,
                                   const float density, const float u) {

  /* Set the density to be used in the force loop to be the evolved density. */
  p->rho = p->strength_data.rho_evol;

#ifdef PLANETARY_FIXED_ENTROPY
  /* Override the internal energy to satisfy the fixed entropy.
    * Needs to happen here as well because of the updated rho. */
  p->u = gas_internal_energy_from_entropy(p->rho, p->s_fixed, p->mat_id);
  xp->u_full = p->u;
#endif

  const float pressure =
      gas_pressure_from_internal_energy(p->rho, p->u, p->mat_id);

  /* Compute stress tensor. */
  strength_compute_stress_tensor(p, pressure);

  /* Compute principal stresses. */
  sym_matrix_compute_eigenvalues(p->strength_data.principal_stress_eigen, p->strength_data.stress_tensor);
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

  memset(p->strength_data.dv_force_loop, 0.f, 3 * 3 * sizeof(float));
  zero_sym_matrix(&p->strength_data.dS_dt);
}

/**
 * @brief Finishes extra strength parts of the force calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_force_strength(struct part *restrict p) {

 /* Update dS/dt for timestep. */
 stress_tensor_compute_dS_dt(p, p->strength_data.dv_force_loop);

  /* Get quntities needed for dD/dt calculation. */
  const int mat_id = p->mat_id;
  const int phase_state = p->phase_state;
  const float mass = p->mass;
  const float density = p->rho_evol;
  const float u = p->u;
  const float pressure = gas_pressure_from_internal_energy(density, u, mat_id);
  const float damage = strength_get_damage(p);
  const float yield_stress = yield_compute_yield_stress(mat_id, phase_state, density, u, damage);
  const struct sym_matrix deviatoric_stress_tensor = p->strength_data.deviatoric_stress_tensor;

  struct sym_matrix stress_tensor;
  const struct sym_matrix damaged_deviatoric_stress_tensor = yield_compute_damaged_deviatoric_stress_tensor(deviatoric_stress_tensor, damage);
  damage_compute_stress_tensor(&stress_tensor, damaged_deviatoric_stress_tensor, pressure, damage);

  /* Update dD/dt for timestep. */
  damage_compute_dD_dt(p, stress_tensor, deviatoric_stress_tensor, mat_id, mass, density, u, yield_stress);
}

/**
 * @brief Sets the values of additional particle strength properties at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values_strength(
    struct part *restrict p, const struct xpart *restrict xp) {

  p->rho = xp->strength_data.rho_evol_full;
  p->strength_data.rho_evol = xp->strength_data.rho_evol_full;

  strength_reset_predicted_values_stress_tensor(p, xp);
  strength_reset_predicted_values_damage(p, xp);
  strength_reset_predicted_values_extra(p, xp);
}

/**
 * @brief Predict additional particle strength properties forward in time when
 * drifting. At beginning of hydro function, before hydro quantities have been drifted.
 *
 * @param p The particle to act upon
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_strength_beginning(
    struct part *restrict p, const float dt_therm) {

  // ### FOR LEAPFROG: dS_dt is calculated similarly to e.g. du_dt and S is updated similarly to u
  // ### dDamage_dt does not depend on positions or velocities of particle neighbours, only on properties
  // ### of the particle itself. Therefore dDamage_dt is recalculated each time damage is updated.
  // ## Damage depends on S and S can depend on damage through Y makes.

  /* Get quantities needed for strength evolution. */
  const int mat_id = p->mat_id;
  const int phase_state = p->phase_state;
  const float mass = p->mass;
  const float density = p->strength_data.rho_evol;
  const float u = p->u;
  const float pressure = gas_pressure_from_internal_energy(density, u, mat_id);
  const float damage = strength_get_damage(p);
  const float yield_stress = yield_compute_yield_stress(mat_id, phase_state, density, u, damage);
  const struct sym_matrix deviatoric_stress_tensor = p->strength_data.deviatoric_stress_tensor;

  struct sym_matrix stress_tensor;
  const struct sym_matrix damaged_deviatoric_stress_tensor = yield_compute_damaged_deviatoric_stress_tensor(deviatoric_stress_tensor, damage);
  damage_compute_stress_tensor(stress_tensor, damaged_deviatoric_stress_tensor, pressure, damage);

  // ### Since I have to calc dD/d now for timesteps, could this just use p->dD/dt?
  /* Evolve damage. */
  damage_predict_evolve(&stress_tensor, deviatoric_stress_tensor, mat_id, mass,
                              density, u, yield_stress, dt_therm);

  /* Evolve deviatoric stress tensor. */
  stress_tensor_evolve_deviatoric_stress_tensor(&p->strength_data.deviatoric_stress_tensor, p, phase_state, dt_therm);

  /* Apply yield stress to deviatoric stress tensor. */
  yield_apply_yield_stress_to_deviatoric_stress_tensor(
        &p->strength_data.deviatoric_stress_tensor, yield_stress, density, u);

  strength_predict_extra_beginning(p, dt_therm);
}

/**
 * @brief Predict additional particle strength properties forward in time when
 * drifting. At end of hydro function, after hydro quantities have been drifted.
 *
 * @param p The particle to act upon
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_strength_end(
    struct part *restrict p, const float dt_therm) {

  /* Evolve density. */
  p->strength_data.rho_evol += p->strength_data.drho_dt * dt_therm;

  /* Compute minimum density */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  const float min_rho = p->mass * kernel_root * h_inv_dim;

  /* Overwrite stored hydro qunatities with those calculated based on evolved density. */
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

  /* Compute updated stress tensor. */
  strength_compute_stress_tensor(p, pressure);
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
__attribute__((always_inline)) INLINE static void hydro_kick_strength_beginning(
    struct part *restrict p, struct xpart *restrict xp, const float dt_therm) {

  /* Get quantities needed for strength evolution. */
  const int mat_id = p->mat_id;
  const int phase_state = xp->phase_state_full;
  const float mass = p->mass;
  const float density = xp->strength_data.rho_evol_full;
  const float u = xp->u_full;
  const float pressure = gas_pressure_from_internal_energy(density, u, mat_id);
  const float damage = strength_get_damage_full(xp);
  const float yield_stress = yield_compute_yield_stress(mat_id, phase_state, density, u, damage);
  const struct sym_matrix deviatoric_stress_tensor = xp->strength_data.deviatoric_stress_tensor_full;

  struct sym_matrix stress_tensor;
  const struct sym_matrix damaged_deviatoric_stress_tensor = yield_compute_damaged_deviatoric_stress_tensor(deviatoric_stress_tensor, damage);
  damage_compute_stress_tensor(&stress_tensor, damaged_deviatoric_stress_tensor, pressure, damage);

  /* Evolve damage. */
  damage_kick_evolve(stress_tensor, deviatoric_stress_tensor, mat_id, mass,
                           density, u, yield_stress, dt_therm);

  /* Evolve deviatoric stress tensor. */
  stress_tensor_evolve_deviatoric_stress_tensor(&xp->strength_data.deviatoric_stress_tensor_full, p, phase_state, dt_therm);

  /* Apply yield stress to deviatoric stress tensor. */
  yield_apply_yield_stress_to_deviatoric_stress_tensor(
        &xp->strength_data.deviatoric_stress_tensor_full, yield_stress, density, u);

  strength_kick_extra_beginning(p, xp, dt_therm);
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
__attribute__((always_inline)) INLINE static void hydro_kick_strength_end(
    struct part *restrict p, struct xpart *restrict xp, const float dt_therm) {

  /* Evolve density. Note this comes after e.g. calculation of stress tensor for
   * strength evolution, since the stress tensor used in the evolution myst be
   * at the current time. */
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

  /* Overwrite based on evolved density. */
  xp->phase_state_full =
    (enum mat_phase_state)material_phase_state_from_internal_energy(
     xp->strength_data.rho_evol_full, xp->u_full, p->mat_id);
}

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

  strength_first_init_part_stress_tensor(p, xp);
  strength_first_init_part_damage(p, xp);
}

#endif /* SWIFT_PLANETARY_STRENGTH_DEFAULT_H */
