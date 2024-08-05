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

/**
 * @brief Set the (symmetric) stress tensor by combining the deviatoric with the pressure.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_set_stress_tensor(
    struct part *restrict p, const float pressure) {

    for (int i = 0; i < 6; i++) {
        p->stress_tensor.elements[i] = p->deviatoric_stress_tensor.elements[i];
    }

    float pressure_hat = pressure;

    /*
    // If strength, if damage
    // See schafer 2016 for this...
    if (pressure < 0.f){
        pressure_hat *= (1 - p->damage_D);
    }
    */

    p->stress_tensor.xx -= pressure_hat;
    p->stress_tensor.yy -= pressure_hat;
    p->stress_tensor.zz -= pressure_hat;

    // Compute principal stresses
    sym_matrix_compute_eigenvalues(p->principal_stresses, p->stress_tensor);    
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
 * @brief Calculates the stress tensor with strength for the force interaction.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_set_pairwise_stress_tensors_strength(
    float pairwise_stress_tensor_i[3][3], float pairwise_stress_tensor_j[3][3],
    const struct part *restrict pi, const struct part *restrict pj, const float r) {
    
    // Use the full stress tensor for solid particles
    if ((pi->phase_state == eos_phase_state_solid) &&
        (pj->phase_state == eos_phase_state_solid)) {

      get_matrix_from_sym_matrix(pairwise_stress_tensor_i, &pi->stress_tensor);
      get_matrix_from_sym_matrix(pairwise_stress_tensor_j, &pj->stress_tensor);

      // Artificial stress (Monaghan (2000))
      /*  
      const float mean_h = 0.5f * (pi->h + pj->h);
      const float delta_p = 0.5f * (powf(pi->mass / pi->rho, 1/hydro_dimension) + 
                                    powf(pj->mass / pj->rho, 1/hydro_dimension)); 
        
      float wij_delta_p;
      kernel_eval(delta_p / mean_h, &wij_delta_p);
    
      float wij_r;
      kernel_eval(r / mean_h, &wij_r);

      // This factor should be set in extra parameter file  
      const float n = 4.f;  
      const float f_factor = powf(wij_r / wij_delta_p, n);
    
      // This factor should be set in extra parameter file
      const float stress_epsilon = 0.2f;  

      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          if (pairwise_stress_tensor_i[i][j] > 0) 
            pairwise_stress_tensor_i[i][j] -= f_factor * stress_epsilon * 
                                                pairwise_stress_tensor_i[i][j];
          if (pairwise_stress_tensor_j[i][j] > 0) 
            pairwise_stress_tensor_j[i][j] -= f_factor * stress_epsilon * 
                                                pairwise_stress_tensor_j[i][j];
        }
      } 
      */

      // New version that is independent of basis:
      // Principal stresses are basis-independent  
      // Adding the same constant to all diagonal elemenets is equivalent in any basis:
      // M + c*I = P (M + c*I) P^-1 = P M P^-1 + P (c*I) P^-1 = P M P^-1 + c*I
      // This is a lot easier than e.g. Gray, Monaghan, and Swift (2001)  


      const float mean_h = 0.5f * (pi->h + pj->h);
      const float delta_p = 0.5f * (powf(pi->mass / pi->rho, 1/hydro_dimension) + 
                                    powf(pj->mass / pj->rho, 1/hydro_dimension)); 
        
      float wij_delta_p;
      kernel_eval(delta_p / mean_h, &wij_delta_p);
    
      float wij_r;
      kernel_eval(r / mean_h, &wij_r);

      // This factor should be set in extra parameter file  
      const float n = 4.f;  
      const float f_factor = powf(wij_r / wij_delta_p, n);
    
      // This factor should be set in extra parameter file
      const float stress_epsilon = 0.2f; 
        
      float max_principal_stress_i = pi->principal_stresses[0];
      if (pi->principal_stresses[1] > max_principal_stress_i)
        max_principal_stress_i = pi->principal_stresses[1];
      if (pi->principal_stresses[2] > max_principal_stress_i)
        max_principal_stress_i = pi->principal_stresses[2];  

      float max_principal_stress_j = pj->principal_stresses[0];
      if (pj->principal_stresses[1] > max_principal_stress_j)
        max_principal_stress_j = pj->principal_stresses[1];
      if (pj->principal_stresses[2] > max_principal_stress_j)
        max_principal_stress_j = pj->principal_stresses[2];    

      for (int i = 0; i < 3; ++i) {
        if (max_principal_stress_i > 0) 
            pairwise_stress_tensor_i[i][i] -= f_factor * stress_epsilon * 
                                                max_principal_stress_i;
        if (max_principal_stress_j > 0) 
            pairwise_stress_tensor_j[i][i] -= f_factor * stress_epsilon * 
                                                max_principal_stress_j;
      }
    }
}

/**
 * @brief Prepares extra strength parameters for a particle for the force
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_force_extra_strength(struct part *restrict p, const float pressure) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      p->dv_lin_repr_kernel[i][j] = 0.f;
    }
  }

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
  if ((pi->phase_state == eos_phase_state_solid) &&
      (pj->phase_state == eos_phase_state_solid)) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        pi->dv_lin_repr_kernel[i][j] +=
            (pj->v[j] - pi->v[j]) * Gi[i] * (pj->mass / pj->rho_evol);
        pj->dv_lin_repr_kernel[i][j] +=
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
  if ((pi->phase_state == eos_phase_state_solid) &&
      (pj->phase_state == eos_phase_state_solid)) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        pi->dv_lin_repr_kernel[i][j] +=
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

  float strain_rate_tensor[3][3], rotation_rate_tensor[3][3], rotation_term[3][3];
  float deviatoric_stress_tensor[3][3];
  get_matrix_from_sym_matrix(deviatoric_stress_tensor, &p->deviatoric_stress_tensor);  

  // Set the strain and rotation rates
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_rate_tensor[i][j] =
          0.5f * (p->dv_lin_repr_kernel[i][j] + p->dv_lin_repr_kernel[j][i]);
      rotation_rate_tensor[i][j] =
          0.5f * (p->dv_lin_repr_kernel[j][i] - p->dv_lin_repr_kernel[i][j]);
    }
  }

  // Set rotation to transform into the solid body rotation frame
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rotation_term[i][j] = 0.f;

      for (int k = 0; k < 3; k++) {
        // See Dienes 1978 (eqn 4.8)
        rotation_term[i][j] += deviatoric_stress_tensor[i][k] *
                                   rotation_rate_tensor[k][j] -
                               rotation_rate_tensor[i][k] *
                                   deviatoric_stress_tensor[k][j];
      }
    }
  }

  // Compute time derivative of the deviatoric stress tensor (Hooke's law)
  float shear_modulus = material_shear_mod(p->mat_id);
  float dS_dt[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dS_dt[i][j] = 2.0f * shear_modulus * strain_rate_tensor[i][j] + rotation_term[i][j];
    }

    dS_dt[i][i] -= 2.0f * shear_modulus *
        (strain_rate_tensor[0][0] + strain_rate_tensor[1][1] + strain_rate_tensor[2][2]) / 3.f;
  }

  get_sym_matrix_from_matrix(&p->dS_dt, dS_dt);
}

/**
 * @brief Evolve the deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_deviatoric_stress(
    struct part *restrict p, const float dt_therm) {

  if (p->phase_state == eos_phase_state_fluid) {
    // No stress for fluids
    zero_sym_matrix(&p->deviatoric_stress_tensor);
      
  } else {
    // Update solid stress
    for (int i = 0; i < 6; i++) {
      p->deviatoric_stress_tensor.elements[i] += p->dS_dt.elements[i] * dt_therm;
    }
  }
}

/**
 * @brief Predict additional particle strength properties forward in time when
 * drifting
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra_strength(
    struct part *restrict p, const float dt_therm) {

    evolve_deviatoric_stress(p, dt_therm);
}

#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_STRENGTH_H */
