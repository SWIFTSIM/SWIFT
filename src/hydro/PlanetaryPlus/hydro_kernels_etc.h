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
#ifndef SWIFT_PLANETARY_HYDRO_KERNELS_ETC_H
#define SWIFT_PLANETARY_HYDRO_KERNELS_ETC_H

/**
 * @file Planetary/hydro_kernels_etc.h
 * @brief Utilities for hydro kernels and e.g. gradient estimates.
 */

#include "const.h"
#include "hydro_misc_utils.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Prepares extra kernel parameters for a particle for the density
 * calculation.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part_extra_kernel(
    struct part *restrict p) {     
    
  p->sum_grad_w[0] = 0.f;
  p->sum_grad_w[1] = 0.f; 
  p->sum_grad_w[2] = 0.f;   
    
    p->sum_dw_dh = 0.f;
    
    p->m0_density_loop = 0.f;
    
    
    p->grad_m0_density_loop[0] = 0.f;
  p->grad_m0_density_loop[1] = 0.f; 
  p->grad_m0_density_loop[2] = 0.f;   
    
}

/**
 * @brief Extra kernel density interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_density_extra_kernel(struct part *restrict pi,
                                       struct part *restrict pj,
                                       const float dx[3], const float wi,
                                       const float wj, const float wi_dx,
                                       const float wj_dx) {

      /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;
    
  pi->sum_grad_w[0] += dx[0] * wi_dx * r_inv;
  pi->sum_grad_w[1] += dx[1] * wi_dx * r_inv;
  pi->sum_grad_w[2] += dx[2] * wi_dx * r_inv;   
    
  pj->sum_grad_w[0] += -dx[0] * wj_dx * r_inv;
  pj->sum_grad_w[1] += -dx[1] * wj_dx * r_inv;
  pj->sum_grad_w[2] += -dx[2] * wj_dx * r_inv;   
    
    pi->sum_dw_dh += -(hydro_dimension * wi + r / pi->h * wi_dx); 
  pj->sum_dw_dh += -(hydro_dimension * wj + r / pj->h * wj_dx); 
    
    pi->m0_density_loop += pj->mass * wi / pj->rho_evolved;
    pj->m0_density_loop += pi->mass * wj / pi->rho_evolved;
    
    for (int i = 0; i < 3; i++) {
      pi->grad_m0_density_loop[i] += (pj->mass / pj->rho_evolved) * dx[i] * wi_dx * r_inv;
      pj->grad_m0_density_loop[i] += (pi->mass / pi->rho_evolved) * (-dx[i] * wj_dx * r_inv);
    }
}

/**
 * @brief Extra kernel density interaction between two particles (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_density_extra_kernel(struct part *restrict pi,
                                              const struct part *restrict pj,
                                              const float dx[3], const float wi,
                                              const float wi_dx) {

   /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;
    
  pi->sum_grad_w[0] += dx[0] * wi_dx * r_inv;
  pi->sum_grad_w[1] += dx[1] * wi_dx * r_inv;
  pi->sum_grad_w[2] += dx[2] * wi_dx * r_inv;     
    
    pi->sum_dw_dh += -(hydro_dimension * wi + r / pi->h * wi_dx);  
    
    
    pi->m0_density_loop += pj->mass * wi / pj->rho_evolved;
    
        for (int i = 0; i < 3; i++) {
      pi->grad_m0_density_loop[i] += (pj->mass / pj->rho_evolved) * dx[i] * wi_dx * r_inv;
    }
    
}

/**
 * @brief Finishes extra kernel parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_density_extra_kernel(struct part *restrict p) {

    
    
    const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */
 
  p->sum_grad_w[0] *= h_inv_dim_plus_one;
  p->sum_grad_w[1] *= h_inv_dim_plus_one;
  p->sum_grad_w[2] *= h_inv_dim_plus_one;
  
  p->sum_dw_dh += -(hydro_dimension * kernel_root); 
  p->sum_dw_dh *= h_inv_dim_plus_one;
  
  
  
  float sph_volume = 1.f / p->density.wcount; 
  
  float grad_volume[3];
  grad_volume[0] = -p->sum_grad_w[0] * sph_volume * sph_volume;
  grad_volume[1] = -p->sum_grad_w[1] * sph_volume * sph_volume;
  grad_volume[2] = -p->sum_grad_w[2] * sph_volume * sph_volume;
  
  p->grad_h[0] = p->h * grad_volume[0] / (hydro_dimension * sph_volume);
  p->grad_h[1] = p->h * grad_volume[1] / (hydro_dimension * sph_volume);
  p->grad_h[2] = p->h * grad_volume[2] / (hydro_dimension * sph_volume);
  
  float total_derivative_factor = -p->h * sph_volume * p->sum_dw_dh / hydro_dimension;
  
  p->grad_h[0] /= (1 - total_derivative_factor);
  p->grad_h[1] /= (1 - total_derivative_factor);
  p->grad_h[2] /= (1 - total_derivative_factor);
  
  
  p->m0_density_loop += p->mass * kernel_root / p->rho_evolved;
  p->m0_density_loop *= h_inv_dim;
  
   for (int i = 0; i < 3; i++) {
      p->grad_m0_density_loop[i] *= h_inv_dim_plus_one;
    }
}


/**
 * @brief Prepares extra kernel parameters for a particle for the gradient
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_extra_kernel(struct part *restrict p) {

#if defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      p->Cinv[i][j] = 0.f;
    }
  }
#endif



    p->m0 = 0.f;
    for (int i = 0; i < 3; i++) {
      p->m1[i] = 0.f;
      p->grad_m0[i] = 0.f;
      p->grad_m0_gradhterm[i] = 0.f;
      for (int j = 0; j < 3; j++) {
        p->m2[i][j] = 0.f;
        p->grad_m1_term1[i][j] = 0.f;
        p->grad_m1_term1_gradhterm[i][j] = 0.f;
        p->grad_m1_term2[i][j] = 0.f;
        for (int k = 0; k < 3; k++) {
            p->grad_m2_term1[i][j][k] = 0.f;
            p->grad_m2_term1_gradhterm[i][j][k] = 0.f;
            p->grad_m2_term2[i][j][k] = 0.f;
        }
      }
    }   

}

/**
 * @brief Extra kernel gradient interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_gradient_extra_kernel(struct part *restrict pi,
                                        struct part *restrict pj,
                                        const float dx[3], const float wi,
                                        const float wj, const float wi_dx,
                                        const float wj_dx) {

#if defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC

  float sph_volume_i = pi->mass / pi->sph_rho;
  float sph_volume_j = pj->mass / pj->sph_rho;

  //planetary_smoothing_correction_tweak_volume(&volume_i, pi);
 // planetary_smoothing_correction_tweak_volume(&volume_j, pj);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      /* Inverse of C matrix (eq 6 in Rosswog 2020) */
      pi->Cinv[i][j] += dx[i] * dx[j] * wi * sph_volume_j;
      pj->Cinv[i][j] += dx[i] * dx[j] * wj * sph_volume_i;
    }
  }
#endif /* defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC */
    
    
    
     float volume_i = pi->mass / pi->rho_evolved;
float volume_j = pj->mass / pj->rho_evolved;   
    
      const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;
    
    
      const float hi = pi->h;
  const float hi_inv = 1.0f / hi;                 /* 1/h */
  const float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
    const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */
    
    
      const float hj = pj->h;
  const float hj_inv = 1.0f / hj;                 /* 1/h */
  const float hj_inv_dim = pow_dimension(hj_inv); /* 1/h^d */
    const float hj_inv_dim_plus_one = hj_inv_dim * hj_inv; /* 1/h^(d+1) */
    
    
    
float wi_term = 0.5f * (wi * hi_inv_dim + wj * hj_inv_dim);    
float wj_term = 0.5f * (wi * hi_inv_dim + wj * hj_inv_dim);       
float wi_dx_term[3], wj_dx_term[3];
for (int i = 0; i < 3; i++) {   
    wi_dx_term[i] = dx[i] * r_inv * 0.5f * (wi_dx * hi_inv_dim_plus_one + wj_dx * hj_inv_dim_plus_one);    
    wj_dx_term[i] = -dx[i] * r_inv * 0.5f * (wi_dx * hi_inv_dim_plus_one + wj_dx * hj_inv_dim_plus_one);           
}
    
    // h term
float wi_dx_gradhterm[3], wj_dx_gradhterm[3];
  for (int i = 0; i < 3; i++) {   
    wi_dx_gradhterm[i] = -0.5f *(hydro_dimension * wi + (r / hi) * wi_dx) * hi_inv_dim_plus_one; 
    wj_dx_gradhterm[i] = -0.5f *(hydro_dimension * wj + (r / hj) * wj_dx) * hj_inv_dim_plus_one;  
}      

    
pi->m0 += volume_j * wi_term;
pj->m0 += volume_i * wj_term;
for (int i = 0; i < 3; i++) {
  pi->m1[i] += dx[i] * volume_j * wi_term;
  pj->m1[i] += -dx[i] * volume_i * wj_term;

  pi->grad_m0[i] += volume_j * wi_dx_term[i];
  pj->grad_m0[i] += volume_i * wj_dx_term[i];
    
    pi->grad_m0_gradhterm[i] += volume_j * wi_dx_gradhterm[i];
  pj->grad_m0_gradhterm[i] += volume_i * wj_dx_gradhterm[i];
  for (int j = 0; j < 3; j++) {
    pi->m2[i][j] += dx[i] * dx[j] * volume_j * wi_term;
    pj->m2[i][j] += dx[i] * dx[j] * volume_i * wj_term;

    pi->grad_m1_term1[i][j] += volume_j * dx[j] * wi_dx_term[i];
    pj->grad_m1_term1[i][j] += -volume_i * dx[j] * wj_dx_term[i];
      
      pi->grad_m1_term1_gradhterm[i][j] += volume_j * dx[j] * wi_dx_gradhterm[i];
    pj->grad_m1_term1_gradhterm[i][j] += -volume_i * dx[j] * wj_dx_gradhterm[i];
    if (i == j){
      pi->grad_m1_term2[i][j] += volume_j * wi_term;
      pj->grad_m1_term2[i][j] += volume_i * wj_term;
    }

    for (int k = 0; k < 3; k++) {
        pi->grad_m2_term1[i][j][k] += volume_j * dx[k] * dx[j] * wi_dx_term[i];
        pj->grad_m2_term1[i][j][k] += volume_i * dx[k] * dx[j] * wj_dx_term[i];
        
        pi->grad_m2_term1_gradhterm[i][j][k] += volume_j * dx[k] * dx[j] * wi_dx_gradhterm[i];
        pj->grad_m2_term1_gradhterm[i][j][k] += volume_i * dx[k] * dx[j] * wj_dx_gradhterm[i];

        if (i == j){
          pi->grad_m2_term2[i][j][k] += volume_j * dx[k] * wi_term;
          pj->grad_m2_term2[i][j][k] += -volume_i * dx[k] * wj_term;
        }
        if (i == k){
          pi->grad_m2_term2[i][j][k] += volume_j * dx[j] * wi_term;
          pj->grad_m2_term2[i][j][k] += -volume_i * dx[j] * wj_term;
        }
    }
  }
}
}

/**
 * @brief Extra kernel gradient interaction between two particles
 * (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_gradient_extra_kernel(struct part *restrict pi,
                                               struct part *restrict pj,
                                               const float dx[3], const float wi,
                                        const float wj, const float wi_dx,
                                        const float wj_dx) {

#if defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC

  float sph_volume_j = pj->mass / pj->sph_rho;

  //planetary_smoothing_correction_tweak_volume(&volume_j, pj);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      /* Inverse of C matrix (eq 6 in Rosswog 2020) */
      pi->Cinv[i][j] += dx[i] * dx[j] * wi * sph_volume_j;
    }
  }

#endif /* defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC */
    
    
float volume_j = pj->mass / pj->rho_evolved;  
    
      const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;
    
      const float hi = pi->h;
  const float hi_inv = 1.0f / hi;                 /* 1/h */
  const float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
    const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */
    
    
      const float hj = pj->h;
  const float hj_inv = 1.0f / hj;                 /* 1/h */
  const float hj_inv_dim = pow_dimension(hj_inv); /* 1/h^d */
    const float hj_inv_dim_plus_one = hj_inv_dim * hj_inv; /* 1/h^(d+1) */
    
    
    
float wi_term = 0.5f * (wi * hi_inv_dim + wj * hj_inv_dim);    
float wi_dx_term[3];
for (int i = 0; i < 3; i++) {   
    wi_dx_term[i] = dx[i] * r_inv * 0.5f * (wi_dx * hi_inv_dim_plus_one + wj_dx * hj_inv_dim_plus_one);     
}
 

    // h term
float wi_dx_gradhterm[3];
  for (int i = 0; i < 3; i++) {   
    wi_dx_gradhterm[i] = -0.5f *(hydro_dimension * wi + (r / hi) * wi_dx) * hi_inv_dim_plus_one; 
}      


pi->m0 += volume_j * wi_term;
for (int i = 0; i < 3; i++) {
  pi->m1[i] += dx[i] * volume_j * wi_term;

  pi->grad_m0[i] += volume_j * wi_dx_term[i];
    
     pi->grad_m0_gradhterm[i] += volume_j * wi_dx_gradhterm[i];
  for (int j = 0; j < 3; j++) {
    pi->m2[i][j] += dx[i] * dx[j] * volume_j * wi_term;

    pi->grad_m1_term1[i][j] += volume_j * dx[j] * wi_dx_term[i];
      
    pi->grad_m1_term1_gradhterm[i][j] += volume_j * dx[j] * wi_dx_gradhterm[i];  
    if (i == j){
      pi->grad_m1_term2[i][j] += volume_j * wi_term;
    }

    for (int k = 0; k < 3; k++) {
        pi->grad_m2_term1[i][j][k] += volume_j * dx[k] * dx[j] * wi_dx_term[i];
        
        pi->grad_m2_term1_gradhterm[i][j][k] += volume_j * dx[k] * dx[j] * wi_dx_gradhterm[i];

        if (i == j){
          pi->grad_m2_term2[i][j][k] += volume_j * dx[k] * wi_term;
        }
        if (i == k){
          pi->grad_m2_term2[i][j][k] += volume_j * dx[j] * wi_term;
        }
    }
  }
}

}

/**
 * @brief Finishes extra kernel parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_extra_kernel(struct part *restrict p) {

#if defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC
  int i, j;



  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */


  
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      p->Cinv[i][j] *= h_inv_dim;
    }
  }


  /* Find the inverse of the Cinv matrix */

  /* If h=h_max don't do anything fancy. Things like using m/rho to calculate
   * the volume stops working */
if (!p->is_h_max) {
    // Invert C matrix
    float temp_matrix[3][3];

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        temp_matrix[i][j] = p->Cinv[i][j];
      }
    }
    
    invert_dimension_by_dimension_matrix(temp_matrix);
    
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        p->C[i][j] = temp_matrix[i][j];
      }
    }
    }else{
    
    // this shouldn't be used because of later if statements 
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        p->C[i][j] =0.f;
      }
    }
    
    }
#endif


  
  

   
float volume = p->mass / p->rho_evolved;  
    



int k;
p->m0 += volume * kernel_root * h_inv_dim;
p->grad_m1_term2[0][0] += volume * kernel_root * h_inv_dim;
p->grad_m1_term2[1][1] += volume * kernel_root * h_inv_dim;
p->grad_m1_term2[2][2] += volume * kernel_root * h_inv_dim;

const float h_inv_dim_plus_one = h_inv_dim * h_inv; 
p->grad_m0_gradhterm[0] -= 0.5f * volume * hydro_dimension * kernel_root * h_inv_dim_plus_one;
p->grad_m0_gradhterm[1] -= 0.5f * volume * hydro_dimension * kernel_root * h_inv_dim_plus_one;
p->grad_m0_gradhterm[2] -= 0.5f * volume * hydro_dimension * kernel_root * h_inv_dim_plus_one;


for (i = 0; i < 3; i++) {
    p->grad_m0_gradhterm[i] *= p->dh_sphgrad[i];
    p->grad_m0[i] += p->grad_m0_gradhterm[i];
  for (j = 0; j < 3; j++) {
    p->grad_m1_term1_gradhterm[i][j] *= p->dh_sphgrad[i];
    p->grad_m1_term1[i][j] += p->grad_m1_term1_gradhterm[i][j];
    for (k = 0; k < 3; k++) {
        p->grad_m2_term1_gradhterm[i][j][k] *= p->dh_sphgrad[i];
        p->grad_m2_term1[i][j][k] += p->grad_m2_term1_gradhterm[i][j][k];
    }
  }
}


// Combine terms to get final m expressions
float grad_m1[3][3];
float grad_m2[3][3][3];
for (i = 0; i < 3; i++) {
  for (j = 0; j < 3; j++) {
    grad_m1[i][j] = p->grad_m1_term1[i][j] + p->grad_m1_term2[i][j];
    for (k = 0; k < 3; k++) {
        grad_m2[i][j][k] = p->grad_m2_term1[i][j][k] + p->grad_m2_term2[i][j][k];
    }
  }
}


if (!p->is_h_max) {
    // Invert m2 matrix
    float m2_inv[3][3];

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        m2_inv[i][j] = p->m2[i][j];
      }
    }

    invert_dimension_by_dimension_matrix(m2_inv);


    // Calculate A and B
    p->A = p->m0;
    for (i = 0; i < 3; i++) {
      p->B[i] = 0.f;
      for (j = 0; j < 3; j++) {
          p->A -= m2_inv[i][j] * p->m1[i] * p->m1[j];
          p->B[i] -= m2_inv[i][j] * p->m1[j];
        }
      }

    p->A = 1/p->A;


    // Calculate grad_A and grad_B
    int l, m;
    for (i = 0; i < 3; i++) {
        p->grad_A[i] = p->grad_m0[i];
        for (j = 0; j < 3; j++) {
          p->grad_B[i][j] = 0.f;

          for (k = 0; k < 3; k++) {
            p->grad_A[i] += -m2_inv[j][k] * p->m1[k] * grad_m1[i][j] - m2_inv[j][k] * p->m1[j] * grad_m1[i][k];
            p->grad_B[i][j] += -m2_inv[j][k] * grad_m1[i][k];
            for (l = 0; l < 3; l++) {
                for (m = 0; m < 3; m++) {
                  p->grad_A[i] += m2_inv[j][l] * grad_m2[i][l][m] * m2_inv[m][k] * p->m1[k] * p->m1[j];
                  p->grad_B[i][j] += m2_inv[j][l] * grad_m2[i][l][m] * m2_inv[m][k] * p->m1[k];

                }
              }

          }
        }

        p->grad_A[i] *= - p->A * p->A;

      }

         p->vac_term = 1.f;

     float x = p->h * sqrtf(p->B[0] * p->B[0] + p->B[1] * p->B[1] + p->B[2] * p->B[2]); 
    float offset = 0.8f;
    if (x > offset){
        float sigma = 0.2f;
        p->vac_term = expf(-(x - offset) * (x - offset) / (2.f * sigma * sigma));
    }
  
      
  }else{
      // this shouldn't be used because of later if statements 
      p->A = 1.f;
      p->vac_term = 1.f;
      for (i = 0; i < 3; i++) {
          p->B[i] = 0.f;
          p->grad_A[i] = 0.f;
          for (j = 0; j < 3; j++) {
                p->grad_B[i][j] = 0.f;
            }
      }
  
  
  }
   
  
    
  
}





__attribute__((always_inline)) INLINE static void hydro_set_Gi_Gj_forceloop(
    float Gi[3], float Gj[3], const struct part *restrict pi,
    const struct part *restrict pj, const float dx[3], const float wi,
    const float wj, const float wi_dx, const float wj_dx) {

  /* Get r and 1/r. */
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float hi_inv = 1.0f / pi->h;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */

  const float hj_inv = 1.0f / pj->h;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */

  const float wi_dr = hid_inv * wi_dx;
  const float wj_dr = hjd_inv * wj_dx;

#ifdef PLANETARY_MATRIX_INVERSION
  if (!pi->is_h_max && !pj->is_h_max) {
    for (int i = 0; i < 3; ++i) {
      /* eq 4 and 5 in Rosswog 2020. These replace the gradient of the kernel */
      Gi[i] =
          -(pi->C[i][0] * dx[0] + pi->C[i][1] * dx[1] + pi->C[i][2] * dx[2]) *
          wi;
      Gj[i] =
          -(pj->C[i][0] * dx[0] + pj->C[i][1] * dx[1] + pj->C[i][2] * dx[2]) *
          wj;
    }
  } else {
    for (int i = 0; i < 3; ++i) {
      /* If h=h_max use the standard kernel gradients */
      Gi[i] = wi_dr * dx[i] * r_inv;
      Gj[i] = wj_dr * dx[i] * r_inv;
    }
  }

    
//#elif PLANETARY_GDF
  /* Standard GDF kernel gradients, Wadsley+2017 Eqn. 7, in Rosswog2020
   * framework */
  /* Include the dx and r_inv here instead of later */
//  for (int i = 0; i < 3; i++) {
  //  Gi[i] = wi_dr * dx[i] * r_inv * pi->f_gdf;
   // Gj[i] = wj_dr * dx[i] * r_inv * pj->f_gdf;
  //}
#else  /* !PLANETARY_MATRIX_INVERSION, !PLANETARY_GDF */
//elif PLANETARY_CRKSPH
   if (!pi->is_h_max && !pj->is_h_max) {
         int i, j;
      float modified_grad_wi[3];
      float modified_grad_wj[3];

      const float hi_inv_dim = pow_dimension(hi_inv);       /* 1/h^d */
      const float hj_inv_dim = pow_dimension(hj_inv);       /* 1/h^d */
        // kernels with h factors
       float wi_h = wi * hi_inv_dim;
       float wj_h = wj * hj_inv_dim;


        float wi_term = 0.5f * (wi_h + wj_h);    
        float wj_term = 0.5f * (wi_h + wj_h);       
        float wi_dx_term[3], wj_dx_term[3];
        for (i = 0; i < 3; i++) {   
            wi_dx_term[i] = dx[i] * r_inv * 0.5f * (wi_dx * hid_inv + wj_dx * hjd_inv);    
            wj_dx_term[i] = -dx[i] * r_inv * 0.5f * (wi_dx * hid_inv + wj_dx * hjd_inv);           
        }

        // h term

      for (i = 0; i < 3; i++) {   
        wi_dx_term[i] +=-0.5f *(hydro_dimension * wi + (r / pi->h) * wi_dx) * hid_inv * pi->dh_sphgrad[i];
        wj_dx_term[i] += -0.5f *(hydro_dimension * wj + (r / pj->h) * wj_dx) * hjd_inv * pj->dh_sphgrad[i];
    }  
    


      for (i = 0; i < 3; i++) {
        modified_grad_wi[i] = pi->A * wi_dx_term[i] + pi->grad_A[i] * wi_term + pi->A * pi->B[i] * wi_term;
        modified_grad_wj[i] = pj->A * wj_dx_term[i] + pj->grad_A[i] * wj_term + pj->A * pj->B[i] * wj_term;
      }

      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
        modified_grad_wi[i] += pi->A * pi->B[j] * dx[j] * wi_dx_term[i];
        modified_grad_wi[i] += pi->grad_A[i] * pi->B[j] * dx[j] * wi_term;
        modified_grad_wi[i] += pi->A * pi->grad_B[i][j] * dx[j] * wi_term;

        modified_grad_wj[i] += -pj->A * pj->B[j] * dx[j] * wj_dx_term[i];
         modified_grad_wj[i] += -pj->grad_A[i] *  pj->B[j] * dx[j] * wj_term;
          modified_grad_wj[i] += -pj->A * pj->grad_B[i][j] * dx[j]  * wj_term;
      }
      }



          float modified_wi = pi->A * wi_term;
      float modified_wj = pj->A * wj_term;
      modified_wi += pi->A * pi->B[0] * dx[0] * wi_term + pi->A * pi->B[1] * dx[1] * wi_term + pi->A * pi->B[2] * dx[2] * wi_term;
      modified_wj += -(pj->A * pj->B[0] * dx[0] * wj_term + pj->A * pj->B[1] * dx[1] * wj_term + pj->A * pj->B[2] * dx[2] * wj_term);

      float wi_dx_term_vac[3], wj_dx_term_vac[3];
        for (i = 0; i < 3; i++) {   
            wi_dx_term_vac[i] = dx[i] * r_inv * wi_dx * hid_inv;    
            wj_dx_term_vac[i] = -dx[i] * r_inv * wj_dx * hjd_inv;           
        }

        // h term
      for (i = 0; i < 3; i++) {   
        wi_dx_term_vac[i] += -(hydro_dimension * wi + (r / pi->h) * wi_dx) * hid_inv * pi->dh_sphgrad[i];
        wj_dx_term_vac[i] += -(hydro_dimension * wj + (r / pj->h) * wj_dx) * hjd_inv * pj->dh_sphgrad[i];
    }  
       
        for (i = 0; i < 3; i++) {

        modified_grad_wi[i] *= pi->vac_term;  
        modified_grad_wj[i] *= pj->vac_term;

        modified_grad_wi[i] += wi_dx_term_vac[i];
        modified_grad_wj[i] += wj_dx_term_vac[i];

        modified_grad_wi[i] -= pi->vac_term * wi_dx_term_vac[i];
        modified_grad_wj[i] -= pj->vac_term * wj_dx_term_vac[i];     


      } 


      for (i = 0; i < 3; i++) {
        Gi[i] = modified_grad_wi[i];
        Gj[i] = -modified_grad_wj[i];
      }
       

   }else{

       for (int i = 0; i < 3; i++) {
            Gi[i] = wi_dr * dx[i] * r_inv;
           Gj[i] = wj_dr * dx[i] * r_inv;
      }

   }
    
// else    
  /* Variable smoothing length term */
//  const float f_ij = 1.f - pi->force.f / pj->mass;
//  const float f_ji = 1.f - pj->force.f / pi->mass;

//  for (int i = 0; i < 3; i++) {
//    Gi[i] = wi_dr * dx[i] * r_inv * f_ij;
 //   Gj[i] = wj_dr * dx[i] * r_inv * f_ji;
 // }
#endif /* PLANETARY_MATRIX_INVERSION */
}






/**
 * @brief Returns kernel gradient terms used in evolution equations
 */
__attribute__((always_inline)) INLINE static void
hydro_set_kernel_gradient_terms(const float dx[3], float kernel_gradient_i[3],
                                float kernel_gradient_j[3],
                                float Q_kernel_gradient_i[3],
                                float Q_kernel_gradient_j[3], const float Gi[3],
                                const float Gj[3]) {

#ifdef PLANETARY_GDF
  /* In GDF we use average of Gi and Gj. */
  
  kernel_gradient_i[0] = 0.5f * (Gi[0] + Gj[0]);
  kernel_gradient_i[1] = 0.5f * (Gi[1] + Gj[1]);
  kernel_gradient_i[2] = 0.5f * (Gi[2] + Gj[2]);

  kernel_gradient_j[0] = 0.5f * (Gi[0] + Gj[0]);
  kernel_gradient_j[1] = 0.5f * (Gi[1] + Gj[1]);
  kernel_gradient_j[2] = 0.5f * (Gi[2] + Gj[2]);
#else  /* !PLANETARY_GDF */
  kernel_gradient_i[0] = Gi[0];
  kernel_gradient_i[1] = Gi[1];
  kernel_gradient_i[2] = Gi[2];

  kernel_gradient_j[0] = Gj[0];
  kernel_gradient_j[1] = Gj[1];
  kernel_gradient_j[2] = Gj[2];
#endif /* PLANETARY_GDF */

#ifdef PLANETARY_QUAD_VISC
  Q_kernel_gradient_i[0] = kernel_gradient_i[0];
  Q_kernel_gradient_i[1] = kernel_gradient_i[1];
  Q_kernel_gradient_i[2] = kernel_gradient_i[2];

  Q_kernel_gradient_j[0] = kernel_gradient_j[0];
  Q_kernel_gradient_j[1] = kernel_gradient_j[1];
  Q_kernel_gradient_j[2] = kernel_gradient_j[2];
#else  /* !PLANETARY_QUAD_VISC */
  Q_kernel_gradient_i[0] = 0.5f * (kernel_gradient_i[0] + kernel_gradient_j[0]);
  Q_kernel_gradient_i[1] = 0.5f * (kernel_gradient_i[1] + kernel_gradient_j[1]);
  Q_kernel_gradient_i[2] = 0.5f * (kernel_gradient_i[2] + kernel_gradient_j[2]);

  Q_kernel_gradient_j[0] = 0.5f * (kernel_gradient_i[0] + kernel_gradient_j[0]);
  Q_kernel_gradient_j[1] = 0.5f * (kernel_gradient_i[1] + kernel_gradient_j[1]);
  Q_kernel_gradient_j[2] = 0.5f * (kernel_gradient_i[2] + kernel_gradient_j[2]);
#endif /* PLANETARY_QUAD_VISC */
}

#endif /* SWIFT_PLANETARY_HYDRO_KERNELS_ETC_H */