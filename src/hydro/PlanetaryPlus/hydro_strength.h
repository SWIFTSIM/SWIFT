/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_PLANETARY_HYDRO_STRENGTH_H
#define SWIFT_PLANETARY_HYDRO_STRENGTH_H

/**
 * @file Planetary/hydro_strength.h
 * @brief Material strength calculations, including fracturing, etc.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Prepares extra strength parameters for a particle for the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_init_part_extra_strength(
    struct part *restrict p) {}


/**
 * @brief Extra strength density interaction between two particles
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_runner_iact_density_extra_strength(
    struct part *restrict pi, struct part *restrict pj, const float dx[3], const float wi, const float wj, const float wi_dx, const float wj_dx) {}




/**
 * @brief Extra strength density interaction between two particles (non-symmetric)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_runner_iact_nonsym_density_extra_strength(
    struct part *restrict pi, const struct part *restrict pj, const float dx[3], const float wi, const float wi_dx) {}




/**
 * @brief Finishes extra strength parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_end_density_extra_strength(
    struct part *restrict p) {}



/**
 * @brief Prepares extra strength parameters for a particle for the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_gradient_extra_strength(
    struct part *restrict p) {

// ### If strength  
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
        p->grad_v[i][j] = 0.f;        
    }
  }    
}

/**
 * @brief Extra strength gradient interaction between two particles
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_runner_iact_gradient_extra_strength(
    struct part *restrict pi, struct part *restrict pj, const float dx[3], const float wi, const float wj, const float wi_dx, const float wj_dx) {
    
    
// ### If strength       
#if defined PLANETARY_MATRIX_INVERSION    
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {     
        // Note here using vi - vj in grad
        pi->grad_v[i][j] += (pj->v[i] - pi->v[i]) * dx[j] * wi * (pj->mass / pj->rho);
        pj->grad_v[i][j] += (pj->v[i] - pi->v[i]) * dx[j] * wj * (pi->mass / pi->rho);
    }
  }
#else
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {     
            // Note here using vi - vj in grad
            pi->grad_v[i][j] += (pj->v[i] - pi->v[i]) * (dx[j]*wi_dx*r_inv) * (pj->mass / pj->rho);
            pj->grad_v[i][j] += (pj->v[i] - pi->v[i]) * (dx[j]*wj_dx*r_inv) * (pi->mass / pi->rho);
        }
      }  

#endif    
    
}




/**
 * @brief Extra strength gradient interaction between two particles (non-symmetric)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_runner_iact_nonsym_gradient_extra_strength(
    struct part *restrict pi, const struct part *restrict pj, const float dx[3], const float wi, const float wi_dx) {
    
// ### If strength    
#if defined PLANETARY_MATRIX_INVERSION    
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {     
        // Note here using vi - vj in grad
        pi->grad_v[i][j] += (pj->v[i] - pi->v[i]) * dx[j] * wi * (pj->mass / pj->rho);
    }
  }
#else
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {     
            // Note here using vi - vj in grad
            pi->grad_v[i][j] += (pj->v[i] - pi->v[i]) * (dx[j]*wi_dx*r_inv) * (pj->mass / pj->rho);
        }
      }  
    
#endif    
}




/**
 * @brief Finishes extra strength parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_end_gradient_extra_strength(
    struct part *restrict p) {
    
// ### If strength      
int i, j, k;    
    
#if defined PLANETARY_MATRIX_INVERSION   
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
        p->strain_rate_tensor_epsilon_dot[i][j] = 0.f;
        p->rotation_rate_tensor_R[i][j] = 0.f;
    }
  }
    
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
            // grad_v is missing C matrix multiplication
            p->strain_rate_tensor_epsilon_dot[i][j] += 0.5f * (p->grad_v[i][k] * p->C[k][j] + p->grad_v[j][k] *  p->C[k][i]);
            p->rotation_rate_tensor_R[i][j] += 0.5f * (p->grad_v[i][k] * p->C[k][j] - p->grad_v[j][k] *  p->C[k][i]);
        }
    }
  }
#else
  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */
    
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
        // grad_v is missing h_inv_dim_plus_one
        p->strain_rate_tensor_epsilon_dot[i][j] = 0.5f * (p->grad_v[i][j] + p->grad_v[j][i]) * h_inv_dim_plus_one;
        p->rotation_rate_tensor_R[i][j] = 0.5f * (p->grad_v[i][j] - p->grad_v[j][i]) * h_inv_dim_plus_one;
    }
  }  
    
#endif 
    
    
    
  // Rotation terms
  float rotation_term[3][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
        rotation_term[i][j] = 0.f;
    }
  }
 for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
            //check this one: is S symmetrixc and R antisymmetric
            rotation_term[i][j] += -p->deviatoric_stress_tensor_S[i][k] * p->rotation_rate_tensor_R[k][j] + p->rotation_rate_tensor_R[i][k] * p->deviatoric_stress_tensor_S[k][j] ;
        }
    }
  }      
    
    
 for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
        p->dS_dt[i][j] = 2.0f * p->shear_modulus_mu * p->strain_rate_tensor_epsilon_dot[i][j] + rotation_term[i][j];
        if (i == j){
            p->dS_dt[i][j] -= 2.0f * p->shear_modulus_mu * (p->strain_rate_tensor_epsilon_dot[0][0] + p->strain_rate_tensor_epsilon_dot[1][1] + p->strain_rate_tensor_epsilon_dot[2][2]) / 3.f;   
        }
    }
  }         
}






/**
 * @brief Calculates the yield stress 
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void calculate_yield_stress(
    struct part *restrict p, const float pressure, const float temperature) {
    
    
    // Jutzi 2015 notation
    
    // 0 is placeholder if we ever want to add in Y_0
    float Y_0 = 0.f;
    // From Collins (2004)
    float mu_i = 2.f;
    float mu_d = 0.8f;
    
    float Y_intact = Y_0 + mu_i * pressure / (1.f + (mu_i * pressure) / (p->Y_M - Y_0));
    
    // This is from Emsenhuber+ (2017)
    float xi = 1.2f;
    Y_intact *= tanhf(xi * (p->T_m / temperature - 1.f)); 
    
    float Y_damaged = mu_d * pressure;
    Y_damaged = min(Y_damaged, Y_intact);
    
    p->yield_stress_Y = (1.f - p->damage_D) * Y_intact + p->damage_D * Y_damaged;
    p->yield_stress_Y = min(p->yield_stress_Y, Y_intact);
    
 
}    





/**
 * @brief Predict additional particle strength properties forward in time when drifting
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_deviatoric_stress(
    struct part *restrict p, const float pressure, const float temperature, const float dt_therm) {

    
// ### If strength
    
int i, j;  
if (temperature > p->T_m){
        
   // No deviatoric stress
    for (i = 0; i < 3; i++) {  
        for (j = 0; j < 3; j++) {
            p->deviatoric_stress_tensor_S[i][j] = 0.f;
        }
    }   
    
    // ## Could easily set damage to 0 here if we want to remove damage when melted to reset when solidified
    
    
}else{
    
     for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            p->deviatoric_stress_tensor_S[i][j] += p->dS_dt[i][j]  * dt_therm;
        }
     }
    

    calculate_yield_stress(p, pressure, temperature);
        
    p->J_2 = 0.f;
    for (i = 0; i < 3; i++) {  
        for (j = 0; j < 3; j++) {
            p->J_2 += 0.5f * p->deviatoric_stress_tensor_S[i][j] * p->deviatoric_stress_tensor_S[i][j];
        }
    }
    // B&A
   //float f = min( (p->yield_stress_Y * p->yield_stress_Y) / (3.f * p->J_2), 1);
    
    // Jutzi 2015 / Collins
    float f = min(p->yield_stress_Y / sqrtf(p->J_2), 1);
 
        
    for (i = 0; i < 3; i++) {  
        for (j = 0; j < 3; j++) {
            p->deviatoric_stress_tensor_S[i][j] *= f;
        }      
    }    
}
}    




/**
 * @brief Evolves particle damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage(
    struct part *restrict p, const float soundspeed, const float dt_therm) {
    
    
        // This follows B&A
    
    
    // Find max eignenvalue of stress_tensor_sigma:

    float max_stress_tensor_sigma = 0.f;
    
    // .......
    
    float E = 9 * p->bulk_modulus_K * p->shear_modulus_mu / (3 * p->bulk_modulus_K + p->shear_modulus_mu);
    
    
    float local_scalar_strain = max_stress_tensor_sigma / ((1 - p->damage_D) * E);
    
    
    p->number_of_activated_flaws = 0;
    for (int i = 0; i < p->number_of_flaws; i++) {
        if (local_scalar_strain > p->activation_thresholds_epsilon_act_ij[i]){
            p->number_of_activated_flaws += 1; 
        }
    }
    
    int in_tension = 0;
    // .........
    
    
    p->dD1_3_dt = 0.f;
    if (in_tension){
        float crack_velocity = 0.4f * soundspeed;
        p->dD1_3_dt = crack_velocity * cbrtf(p->rho / p->mass);
    }

 
    float D_max = cbrtf(p->number_of_activated_flaws / p->number_of_flaws);
    float dD = powf(p->dD1_3_dt * dt_therm, 3.f);
    float dD_max = max(0.f, D_max - p->damage_D);
    if (dD > dD_max) 
        dD = dD_max;
    p->damage_D += dD; 

}   
    
    

#endif /* SWIFT_PLANETARY_HYDRO_STRENGTH_H */
