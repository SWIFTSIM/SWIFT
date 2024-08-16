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
#ifndef SWIFT_PLANETARY_STRENGTH_H
#define SWIFT_PLANETARY_STRENGTH_H
#ifdef MATERIAL_STRENGTH

/**
 * @file Planetary/hydro_strength.h
 * @brief REMIX implementation of SPH with material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_kernels.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength_utilities.h"
#include "strength_damage.h"
#include "strength_stress.h"
#include "strength_yield.h"

/**
 * @brief Prepares extra strength parameters for a particle for the density
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_init_part_extra_strength(struct part *restrict p) {}

/**
 * @brief Extra strength density interaction between two particles
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_density_extra_strength(struct part *restrict pi,
                                         struct part *restrict pj,
                                         const float dx[3], const float wi,
                                         const float wj, const float wi_dx,
                                         const float wj_dx) {}

/**
 * @brief Extra strength density interaction between two particles
 * (non-symmetric)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_density_extra_strength(struct part *restrict pi,
                                                const struct part *restrict pj,
                                                const float dx[3],
                                                const float wi,
                                                const float wi_dx) {}

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
 * @brief Extra strength gradient interaction between two particles
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_gradient_extra_strength(struct part *restrict pi,
                                          struct part *restrict pj,
                                          const float dx[3], const float wi,
                                          const float wj, const float wi_dx,
                                          const float wj_dx) {}

/**
 * @brief Extra strength gradient interaction between two particles
 * (non-symmetric)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_gradient_extra_strength(struct part *restrict pi,
                                                 const struct part *restrict pj,
                                                 const float dx[3],
                                                 const float wi,
                                                 const float wi_dx) {}

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
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      p->dv_force_loop[i][j] = 0.f;
    }
  }
    
  const float pressure =
      gas_pressure_from_internal_energy(density, u, p->mat_id);
    
  hydro_set_stress_tensor(p, pressure);
}

/**
 * @brief Extra strength force interaction between two particles
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_force_extra_strength(struct part *restrict pi,
                                       struct part *restrict pj,
                                       const float dx[3], const float Gi[3],
                                       const float Gj[3]) {

  // Compute velocity gradient if both particles are solid
  if ((pi->phase_state == mat_phase_state_solid) &&
      (pj->phase_state == mat_phase_state_solid)) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        pi->dv_force_loop[i][j] +=
            (pj->v[j] - pi->v[j]) * Gi[i] * (pj->mass / pj->rho_evol);
        pj->dv_force_loop[i][j] +=
            (pi->v[j] - pj->v[j]) * Gj[i] * (pi->mass / pi->rho_evol);
      }
    }
  }
}

/**
 * @brief Extra strength force interaction between two particles (non-symmetric)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_force_extra_strength(struct part *restrict pi,
                                              const struct part *restrict pj,
                                              const float dx[3],
                                              const float Gi[3]) {

  // Compute velocity gradient if both particles are solid
  if ((pi->phase_state == mat_phase_state_solid) &&
      (pj->phase_state == mat_phase_state_solid)) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        pi->dv_force_loop[i][j] +=
            (pj->v[j] - pi->v[j]) * Gi[i] * (pj->mass / pj->rho_evol);
      }
    }
  }
}

/**
 * @brief Finishes extra strength parts of the force calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_force_extra_strength(struct part *restrict p) {

 calculate_dS_dt(p); 
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

  p->deviatoric_stress_tensor = xp->deviatoric_stress_tensor_full;
  #if defined(STRENGTH_DAMAGE)  
    p->damage = xp->damage_full;
    p->tensile_damage = xp->tensile_damage_full;
    p->shear_damage = xp->shear_damage_full;
  #endif /* STRENGTH_DAMAGE */   
}
/**
 * @brief Predict additional particle strength properties forward in time when
 * drifting
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra_strength(
    struct part *restrict p, const float dt_therm) {

  // ### FOR LEAPFROG: dS_dt is calculated similarly to e.g. du_dt and S is updated similarly to u
  // ### dDamage_dt does not depend on positions or velocities of particle neighbours, only on properties 
  // ### of the particle itself. Therefore dDamage_dt is recalculated each time damage is updated.
  // ## Damage depends on S and S can depend on damage through Y makes.

  const int phase_state = p->phase_state; 
  const float density = p->rho_evol;
  const float u = p->u;  
     
  #if defined(STRENGTH_DAMAGE)
    const float damage = p->damage;
    const float yield_stress = compute_yield_stress_damaged(p, phase_state, density, u, damage);
    
    evolve_damage(p, &p->tensile_damage, &p->shear_damage, &p->damage, p->deviatoric_stress_tensor, yield_stress, density, u, dt_therm);  
  #else
    const float yield_stress = compute_yield_stress(p, phase_state, density, u);
  #endif /* STRENGTH_DAMAGE */  

  evolve_deviatoric_stress(p, &p->deviatoric_stress_tensor, phase_state, dt_therm);      

  adjust_deviatoric_stress_tensor_by_yield_stress(
        p, &p->deviatoric_stress_tensor, yield_stress, density, u);
}

/**
 * @brief Kick the additional particle strength properties
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra_strength(
    struct part *restrict p, struct xpart *restrict xp, float dt_therm) {

  const int phase_state = xp->phase_state_full;  
  const float density = xp->rho_evol_full;
  const float u = xp->u_full;
    
  #if defined(STRENGTH_DAMAGE)
    const float damage = xp->damage_full;
    const float yield_stress = compute_yield_stress_damaged(p, phase_state, density, u, damage);
    
    evolve_damage(p, &xp->tensile_damage_full, &xp->shear_damage_full, &xp->damage_full, xp->deviatoric_stress_tensor_full, yield_stress, density, u,  dt_therm);  
  #else
      const float yield_stress = compute_yield_stress(p, phase_state, density, u); 
  #endif /* STRENGTH_DAMAGE */  
   
  evolve_deviatoric_stress(p, &xp->deviatoric_stress_tensor_full, phase_state, dt_therm);  

  adjust_deviatoric_stress_tensor_by_yield_stress(
        p, &xp->deviatoric_stress_tensor_full, yield_stress, density, u);
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

  for (int i = 0; i < 6; i++) {
    p->deviatoric_stress_tensor.elements[i] = 0.f;
    xp->deviatoric_stress_tensor_full.elements[i] = 0.f; 
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
#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_STRENGTH_H */
