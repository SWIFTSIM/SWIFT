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

/**
 * @brief Initialises the strength state of a particle from its drifted values.
 *
 * @param s The strength state to be initialised.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void hydro_strength_state_init_drift(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {

  strength_state_init(s, p->mat_id, p->mass,
                      p->rho_evol, p->u,
                      strength_get_damage(p),
                      p->strength_data.deviatoric_stress_tensor,
                      p->strength_data.strain_rate_tensor, mp);
  strength_state_set_method_fields(s, p, mp);
}

/**
 * @brief Initialises the strength state of a particle from its kicked values.
 *
 * @param s The strength state to be initialised.
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void hydro_strength_state_init_kick(
    struct strength_state *s, const struct part *restrict p,
    const struct xpart *restrict xp, const struct mat_params *mp) {

  strength_state_init(s, p->mat_id, p->mass,
                      xp->rho_evol_full, xp->u_full,
                      strength_get_damage_full(xp),
                      xp->strength_data.deviatoric_stress_tensor_full,
                      p->strength_data.strain_rate_tensor, mp);
  strength_state_set_method_fields(s, p, mp);
}
/**
 * @brief Stores the quantities derived from the strength state on the particle.
 *
 * Every quantity on the particle that is derived from the strength state is set
 * here, so that they always come from the same state.
 *
 * @param p The particle of interest.
 * @param s The strength state of the particle, from its drifted values.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void hydro_strength_state_store_drift(
    struct part *restrict p, struct strength_state *s,
    const struct mat_params *mp) {

  p->phase = strength_state_get_phase(s, mp);
  p->strength_data.P_min = strength_state_get_P_min(s, mp);
  p->strength_data.stress_tensor = strength_state_get_stress_tensor(s, mp);

  const float *principal_stresses = strength_state_get_principal_stresses(s, mp);
  p->strength_data.principal_stress_eigen[0] = principal_stresses[0];
  p->strength_data.principal_stress_eigen[1] = principal_stresses[1];
  p->strength_data.principal_stress_eigen[2] = principal_stresses[2];
}

/**
 * @brief Stores the quantities derived from the strength state on the extended
 * particle.
 *
 * @param xp The extended data of the particle of interest.
 * @param s The strength state of the particle, from its kicked values.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void hydro_strength_state_store_kick(
    struct xpart *restrict xp, struct strength_state *s,
    const struct mat_params *mp) {

  xp->phase_full = strength_state_get_phase(s, mp);
}


/**
 * @brief Limits the tension of a pressure from the EoS, by the minimum pressure
 * from the strength state and by damage, depending on the yield stress method.
 *
 * @param p The particle of interest.
 * @param pressure The pressure from the EoS.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_pressure_strength(const struct part *restrict p,
                                     const float pressure) {

  float limited_pressure = fmaxf(pressure, p->strength_data.P_min);
  yield_apply_damage_to_pressure(&limited_pressure,
                                 strength_get_damage(p).total);

  return limited_pressure;
}

/**
 * @brief Updates the hydro (+ strength) time-step based on strength methods.
 *
 *
 * @param dt_cfl The hydro (+ strength) time-step.
 * @param p The particle of interest.
 * @param hydro_properties The hydro properties.
 */
__attribute__((always_inline)) INLINE static void hydro_compute_timestep_strength(
    float *dt_cfl, const struct part *restrict p, const struct hydro_props *restrict hydro_properties) {

  /* Strength state from the drifted values. */
  const struct mat_params *mp = material_params(p->mat_id);
  struct strength_state s;
  hydro_strength_state_init_drift(&s, p, mp);

  /* Update dt_cfl with elastic contribution. */
  strength_compute_timestep_stress_tensor(dt_cfl, &s, mp, strength_get_props());

  /* Update dt_cfl with damage contribution. */
  strength_compute_timestep_damage(dt_cfl, p, strength_get_props());
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

  const struct mat_params *mp = material_params(p->mat_id);

  strength_compute_max_wave_speed_stress_tensor(wave_speed, p, soundspeed, density, mp);
}

/**
 * @brief Whether material strength adds a condition for artificial viscosity
 * velocity reconstruction between a pair of particles.
 *
 * @param pi First particle.
 * @param pj Second particle.
 */
__attribute__((always_inline)) INLINE static int
hydro_visc_no_reconstruction_strength(const struct part *restrict pi,
                                      const struct part *restrict pj) {

  return strength_is_solid_interface(pi, pj);
}

/**
 * @brief Whether material strength excludes a pair of particles from the
 * artificial diffusion of u and rho.
 *
 * @param pi First particle.
 * @param pj Second particle.
 */
__attribute__((always_inline)) INLINE static int
hydro_difn_exclude_pair_strength(const struct part *restrict pi,
                                 const struct part *restrict pj) {

  return strength_is_solid_interface(pi, pj);
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
 * @brief Prepares extra strength parameters for a particle for the gradient
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_strength(struct part *restrict p) {}

/**
 * @brief Finishes extra strength parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_strength(struct part *restrict p) {}

/**
 * @brief Prepares extra strength parameters for a particle for the force
 * calculation.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param density The density
 * @param u The specific internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_force_strength(struct part *restrict p, struct xpart *restrict xp,
                                   const float density, const float u) {

  /* Strength state from the drifted values. */
  const struct mat_params *mp = material_params(p->mat_id);
  struct strength_state s;
  strength_state_init(&s, p->mat_id, p->mass,
                      density, u,
                      strength_get_damage(p),
                      p->strength_data.deviatoric_stress_tensor,
                      p->strength_data.strain_rate_tensor, mp);
  strength_state_set_method_fields(&s, p, mp);

  /* Store the quantities derived from the state. */
  hydro_strength_state_store_drift(p, &s, mp);
}

/**
 * @brief Resets strength time derivative fields in preparation
  * for the sums taking  place in the various force tasks.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_reset_acceleration_strength(struct part *restrict p) {

  memset(p->strength_data.dv_force_loop, 0, 3 * 3 * sizeof(float));
  zero_sym_matrix(&p->strength_data.dS_dt);
}


/**
 * @brief Finishes extra strength parts of the force calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_force_strength(struct part *restrict p) {

  /* Compute strain rate and rotation rate. */
  float strain_rate_tensor[3][3], rotation_rate_tensor[3][3];
  strength_compute_strain_rate_tensor(strain_rate_tensor, p->strength_data.dv_force_loop);
  strength_compute_rotation_rate_tensor(rotation_rate_tensor, p->strength_data.dv_force_loop);
  get_sym_matrix_from_matrix(&p->strength_data.strain_rate_tensor, strain_rate_tensor);

  /* Strength state from the drifted values. */
  const struct mat_params *mp = material_params(p->mat_id);
  struct strength_state s;
  hydro_strength_state_init_drift(&s, p, mp);

  /* Update dS/dt for timestep. */
  stress_tensor_compute_dS_dt(&p->strength_data.dS_dt, rotation_rate_tensor, &s, mp);

  /* Rates of damage accumulation, for the drift, kicks and time-step. */
  damage_compute_rates(p, &s, mp);
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

  strength_reset_predicted_values_stress_tensor(p, xp);
  strength_reset_predicted_values_damage(p, xp);

  /* Store the quantities derived from the state with the reset values. */
  const struct mat_params *mp = material_params(p->mat_id);
  struct strength_state s;
  hydro_strength_state_init_drift(&s, p, mp);
  hydro_strength_state_store_drift(p, &s, mp);
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

  // ### FOR LEAPFROG: dS_dt and the damage rates are calculated similarly to e.g. du_dt, and S and damage are updated similarly to u
  // ## Damage depends on S and S can depend on damage through Y makes.

  /* Strength state from the drifted values, before they are evolved. */
  const struct mat_params *mp = material_params(p->mat_id);
  struct strength_state s;
  hydro_strength_state_init_drift(&s, p, mp);

  /* Evolve damage. */
  damage_predict_evolve(p, dt_therm);

  /* Evolve deviatoric stress tensor. */
  stress_tensor_evolve_deviatoric_stress_tensor(&p->strength_data.deviatoric_stress_tensor,
                                                p->strength_data.dS_dt, &s, mp, dt_therm);

  /* Apply yield stress to the stepped deviatoric stress tensor. The yield
   * stress comes from the state, at the start of the step like the rates, but
   * J_2 must be of the stepped tensor, since that is what gets limited. */
  p->strength_data.is_yielding =
      yield_apply_yield_stress_to_sym_matrix(&p->strength_data.deviatoric_stress_tensor,
                                             p->strength_data.deviatoric_stress_tensor, &s, mp);
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

  /* Strength state from the drifted values. */
  const struct mat_params *mp = material_params(p->mat_id);
  struct strength_state s;
  hydro_strength_state_init_drift(&s, p, mp);

  /* Store the quantities derived from the state. */
  hydro_strength_state_store_drift(p, &s, mp);
}

/**
 * @brief Kick the additional particle strength properties.
 * At beginning of hydro function, before hydro quantities have been kicked.
 *
 * Additional hydrodynamic quantities are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 */
__attribute__((always_inline)) INLINE static void hydro_kick_strength_beginning(
    struct part *restrict p, struct xpart *restrict xp, const float dt_therm) {

  /* Strength state from the kicked values, before they are evolved. */
  const struct mat_params *mp = material_params(p->mat_id);
  struct strength_state s;
  hydro_strength_state_init_kick(&s, p, xp, mp);

  /* Evolve damage. */
  damage_kick_evolve(p, xp, dt_therm);

  /* Evolve deviatoric stress tensor. */
  stress_tensor_evolve_deviatoric_stress_tensor(&xp->strength_data.deviatoric_stress_tensor_full,
                                                p->strength_data.dS_dt, &s, mp, dt_therm);

  /* Apply yield stress to the stepped deviatoric stress tensor. The yield
   * stress comes from the state, at the start of the step like the rates, but
   * J_2 must be of the stepped tensor, since that is what gets limited. */
  yield_apply_yield_stress_to_sym_matrix(&xp->strength_data.deviatoric_stress_tensor_full,
                                         xp->strength_data.deviatoric_stress_tensor_full, &s, mp);
}

/**
 * @brief Kick the additional particle strength properties
 * At end of hydro function, after hydro quantities have been kicked.
 *
 * Additional hydrodynamic quantities are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 */
__attribute__((always_inline)) INLINE static void hydro_kick_strength_end(
    struct part *restrict p, struct xpart *restrict xp, const float dt_therm) {

  /* Strength state from the kicked values. */
  const struct mat_params *mp = material_params(p->mat_id);
  struct strength_state s;
  hydro_strength_state_init_kick(&s, p, xp, mp);

  /* Store the quantities derived from the state. */
  hydro_strength_state_store_kick(xp, &s, mp);
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

  strength_first_init_part_stress_tensor(p, xp);
  strength_first_init_part_damage(p, xp);

  /* Store the quantities derived from the initial state. */
  const struct mat_params *mp = material_params(p->mat_id);
  struct strength_state s;
  hydro_strength_state_init_drift(&s, p, mp);
  hydro_strength_state_store_drift(p, &s, mp);
  hydro_strength_state_init_kick(&s, p, xp, mp);
  hydro_strength_state_store_kick(xp, &s, mp);
}

#endif /* SWIFT_REMIX_STRENGTH_DEFAULT_H */
