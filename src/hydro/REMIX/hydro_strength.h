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
        struct part *restrict p) {  }

    /**
     * @brief Extra strength gradient interaction between two particles
     *
     * @param p The particle to act upon
     */
    __attribute__((always_inline)) INLINE static void hydro_runner_iact_gradient_extra_strength(
        struct part *restrict pi, struct part *restrict pj, const float dx[3], const float wi, const float wj, const float wi_dx, const float wj_dx) {}




    /**
     * @brief Extra strength gradient interaction between two particles (non-symmetric)
     *
     * @param p The particle to act upon
     */
    __attribute__((always_inline)) INLINE static void hydro_runner_iact_nonsym_gradient_extra_strength(
        struct part *restrict pi, const struct part *restrict pj, const float dx[3], const float wi, const float wi_dx) {}




    /**
     * @brief Finishes extra strength parts of the gradient calculation.
     *
     * @param p The particle to act upon
     */
    __attribute__((always_inline)) INLINE static void hydro_end_gradient_extra_strength(
        struct part *restrict p) {}




/**
 * @brief Prepares extra strength parameters for a particle for the force calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force_extra_strength(
          struct part *restrict p) {

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      p->grad_v[i][j] = 0.f;
    }
    }
}

/**
 * @brief Extra strength force interaction between two particles
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_runner_iact_force_extra_strength(
    struct part *restrict pi, struct part *restrict pj, const float dx[3], const float Gi[3], const float Gj[3]) {

    if(pi->temperature < pi->T_m && pj->temperature < pj->T_m){
  for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
          // Note here using vi - vj in grad
          pi->grad_v[i][j] += (pj->v[j] - pi->v[j]) * Gi[i] * (pj->mass / pj->rho_evolved);
          pj->grad_v[i][j] += (pi->v[j] - pj->v[j]) * (-Gj[i]) * (pi->mass / pi->rho_evolved);
      }
    }
    }
 }



/**
 * @brief Extra strength force interaction between two particles (non-symmetric)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_runner_iact_nonsym_force_extra_strength(
    struct part *restrict pi, const struct part *restrict pj, const float dx[3], const float Gi[3]) {

    if(pi->temperature < pi->T_m && pj->temperature < pj->T_m){
          for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < 3; ++j) {
                  // Note here using vi - vj in grad
                  pi->grad_v[i][j] += (pj->v[j] - pi->v[j]) * Gi[i] * (pj->mass / pj->rho_evolved);
              }
            }
          }
 }



/**
 * @brief Finishes extra strength parts of the force calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_end_force_extra_strength(
    struct part *restrict p) {
/*
      const float h = p->h;
          const float h_inv = 1.0f / h;
          const float h_inv_dim = pow_dimension(h_inv);
          const float h_inv_dim_plus_one = h_inv_dim * h_inv;


      float G_selfcontrib[3];
    for (int i = 0; i < 3; ++i) {
        G_selfcontrib[i] = kernel_root * h_inv_dim * p->grad_A[i];

        G_selfcontrib[i] += kernel_root * h_inv_dim * p->A * p->B[i];

        G_selfcontrib[i] -= 0.5f * hydro_dimension * kernel_root * p->dh_sphgrad[i] * h_inv_dim_plus_one * p->A;


    }


      for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
              // Note here using vi - vj in grad
              p->grad_v[i][j] += p->v[j] * G_selfcontrib[i] * (p->mass / p->rho_evolved);
          }
        }
*/

          // ### If strength
          int i, j, k;

            for (i = 0; i < 3; i++) {
              for (j = 0; j < 3; j++) {
                  p->strain_rate_tensor_epsilon_dot[i][j] = 0.5f * (p->grad_v[i][j] + p->grad_v[j][i]);
                  p->rotation_rate_tensor_R[i][j] = 0.5f * (p->grad_v[j][i] - p->grad_v[i][j]);
              }
            }



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
                  //check this one: is S symmetric and R antisymmetric
                  // switched signs to maych Schafer. See Dienes 1978 (eq 4.8) for derivation
                  rotation_term[i][j] += p->deviatoric_stress_tensor_S[i][k] * p->rotation_rate_tensor_R[k][j] - p->rotation_rate_tensor_R[i][k] * p->deviatoric_stress_tensor_S[k][j] ;
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

    // From Collins (2004)
    float mu_i = p->coefficient_friction_intact_mu_i;//2.f;
    float mu_d = p->coefficient_friction_damaged_mu_d;//0.8f;

    
    // Should be able to decrease if negative pressures until Y=0?
    // miluphcuda does this so maybe not wrong?
    float Y_intact = p->Y_0;
    if(pressure > 0){
        Y_intact = p->Y_0 + mu_i * pressure / (1.f + (mu_i * pressure) / (p->Y_M - p->Y_0));
    }
    
    //aluminium hard-coded for now. See Luther et al. 2022 appendix. this should be 0.85xref density.
    float rho_weak = 0.85f * p->rho_0;
    if (p->rho < rho_weak){
        Y_intact *= powf(p->rho / rho_weak, 4.f);
    }

    float Y_damaged = 0.f;
    if(pressure > 0){
        Y_damaged = mu_d * pressure;
    }
    
    //Maybe Y_damaged also needs have a max value indep of Y_intact? See e.g. Güldemeister et al. 2015; Winkler et al. 2018
    
    //one of these might not be needed
    Y_damaged = min(Y_damaged, Y_intact);

    // Temporarily set to 0 for examples
  //  Y_damaged = 0.f;

    p->yield_stress_Y = (1.f - p->damage_D) * Y_intact + p->damage_D * Y_damaged;
    //one of these might not be needed
    p->yield_stress_Y = min(p->yield_stress_Y, Y_intact);


        // This is from Emsenhuber+ (2017)
    // See SENFT AND STEWART 2007 for why this comes here and not just for intact
    if (temperature > p->T_m){
        p->yield_stress_Y  = 0.f;
    }else{    
        float xi = p->thermal_softening_parameter_xi;//1.2f;
        p->yield_stress_Y *= tanhf(xi * (p->T_m / temperature - 1.f));
    }
    
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
    
    p->yield_stress_Y  = 0.f;

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

   // commented out for cylinders
    calculate_yield_stress(p, pressure, temperature);

    p->J_2 = 0.f;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            p->J_2 += 0.5f * p->deviatoric_stress_tensor_S[i][j] * p->deviatoric_stress_tensor_S[i][j];
        }
    }
    // B&A
  // float f = min( (p->yield_stress_Y * p->yield_stress_Y) / (3.f * p->J_2), 1);

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
    struct part *restrict p, const float pressure, const float soundspeed, const float dt_therm) {


    // This follows B&A

    /*
    // Find max eignenvalue of stress_tensor_sigma:
    float eigen_val1, eigen_val2, eigen_val3;
    compute_eigenvalues_symmetric_3x3(&eigen_val1, &eigen_val2, &eigen_val3, p->stress_tensor_sigma);
    float max_stress_tensor_sigma = eigen_val1;
    if (eigen_val2 > max_stress_tensor_sigma)
        max_stress_tensor_sigma = eigen_val2;
    if (eigen_val3 > max_stress_tensor_sigma)
        max_stress_tensor_sigma = eigen_val3;
    */

         // Temp debug
    /*
    float M[3][3];
     M[0][0] = 1.f;
     M[0][1] = 2.f;
     M[0][2] = 3.f;
     M[1][0] = 2.f;
     M[1][1] = 2.f;
     M[1][2] = 1.f;
     M[2][0] = 3.f;
     M[2][1] = 1.f;
     M[2][2] = 3.f;
     compute_eigenvalues_symmetric_3x3(&eigen_val1, &eigen_val2, &eigen_val3, M);
    printf("%f", eigen_val1);
    printf("\n");
    printf("%f", eigen_val2);
    printf("\n");
    printf("%f", eigen_val3);
    printf("\n");
      exit(0);
   */

    /*
    p->dD1_3_dt = 0.f;
    // Damage can only accumulate if in tension (this will have to change if we add fracture under compression)
    // commented out cof cylinders

    // If there are no flaws, we can not accumulate damage
    if(p->number_of_flaws > 0){
        if (max_stress_tensor_sigma > 0.f){

            
            // tenisle damage
            
            float E = 9.f * p->bulk_modulus_K * p->shear_modulus_mu / (3.f * p->bulk_modulus_K + p->shear_modulus_mu);//5.3e10;//

            p->local_scalar_strain = max_stress_tensor_sigma / ((1.f - p->damage_D) * E);

            int number_of_activated_flaws = 0;
            for (int i = 0; i < p->number_of_flaws; i++) {
                if (p->local_scalar_strain > p->activation_thresholds_epsilon_act_ij[i]){
                    number_of_activated_flaws += 1;
                }
            }
            p->number_of_activated_flaws = max(p->number_of_activated_flaws, number_of_activated_flaws);// / (float)p->number_of_flaws;

            // Schafer
            float longitudinal_wave_speed = sqrtf((p->bulk_modulus_K + (4.f / 3.f) * (1.f - p->damage_D) * p->shear_modulus_mu) / p->rho);
            float crack_velocity = 0.4f * longitudinal_wave_speed;

            p->dD1_3_dt = number_of_activated_flaws * crack_velocity * cbrtf(p->rho / p->mass);
            float D_max_cbrt = cbrtf(number_of_activated_flaws / (float)p->number_of_flaws);
            float dD_cbrt = p->dD1_3_dt * dt_therm;
            float dD_max_cbrt = max(0.f, D_max_cbrt - cbrtf(p->tensile_damage));
            if (dD_cbrt > dD_max_cbrt)
                dD_cbrt = dD_max_cbrt;
            float evolved_D_cbrt = cbrtf(p->tensile_damage) + dD_cbrt;
            p->tensile_damage = powf(evolved_D_cbrt, 3.f);
            
            
        }
    }
    */
    
    // Se Collins for this.
    
    if(sqrtf(p->J_2) > p->yield_stress_Y){
        // do I need this or can it happen when p is negative as well?
        if(pressure > 0){    
            // shear damage
            float plastic_strain_at_the_point_of_failure_epsilon_f;
            if(pressure < p->brittle_to_ductile_transition_pressure){
                
                // Is this meant to be inear? maybe not clear in paper
                plastic_strain_at_the_point_of_failure_epsilon_f = 0.04f * (pressure / p->brittle_to_ductile_transition_pressure) + 0.01f;
                
            }else if(pressure < p->brittle_to_plastic_transition_pressure){
                
                float slope = (0.1f - 0.05f) / (p->brittle_to_plastic_transition_pressure - p->brittle_to_ductile_transition_pressure);
                float intercept = p->brittle_to_ductile_transition_pressure - slope * 0.05f;    
                
                plastic_strain_at_the_point_of_failure_epsilon_f = slope * pressure + intercept;


            }else{
                plastic_strain_at_the_point_of_failure_epsilon_f = 1.f;
            }

           
            
            // Do I need to have rotation terms here?
            float deviatoric_strain_rate_tensor[3][3];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    deviatoric_strain_rate_tensor[i][j] = p->strain_rate_tensor_epsilon_dot[i][j];
                    if (i == j){
                        deviatoric_strain_rate_tensor[i][j] -= (p->strain_rate_tensor_epsilon_dot[0][0] + p->strain_rate_tensor_epsilon_dot[1][1] + p->strain_rate_tensor_epsilon_dot[2][2]) / 3.f;
                    }
                }
            }
            
            float strain_rate_invariant = 0.f;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    strain_rate_invariant += 0.5f * deviatoric_strain_rate_tensor[i][j] * deviatoric_strain_rate_tensor[i][j];
                }
            }
            
            strain_rate_invariant = sqrtf(strain_rate_invariant);

            // can this be negative?
            p->shear_damage += max(0.f, (strain_rate_invariant / plastic_strain_at_the_point_of_failure_epsilon_f) * dt_therm);
            
            p->shear_damage = min(1.f, p->shear_damage);
            

        }

    }





        // total damage 
        p->damage_D = min(1.f, p->tensile_damage + p->shear_damage);

}



#endif /* SWIFT_PLANETARY_HYDRO_STRENGTH_H */
