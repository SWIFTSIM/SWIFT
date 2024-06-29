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

    p->m0 = 0.f;
    p->grad_m0_gradhterm = 0.f;
    for (int i = 0; i < 3; i++) {
      p->m1[i] = 0.f;
      p->grad_m0[i] = 0.f;
      p->grad_m1_term1_gradhterm[i] = 0.f;
      for (int j = 0; j < 3; j++) {
        p->m2[i][j] = 0.f;
        p->grad_m1_term1[i][j] = 0.f;
        p->grad_m2_term1_gradhterm[i][j] = 0.f;
        for (int k = 0; k < 3; k++) {
            p->grad_m2_term1[i][j][k] = 0.f;
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
float wi_dx_gradhterm, wj_dx_gradhterm;
wi_dx_gradhterm = -0.5f *(hydro_dimension * wi + (r / hi) * wi_dx) * hi_inv_dim_plus_one;
wj_dx_gradhterm = -0.5f *(hydro_dimension * wj + (r / hj) * wj_dx) * hj_inv_dim_plus_one;



pi->m0 += volume_j * wi_term;
pj->m0 += volume_i * wj_term;

pi->grad_m0_gradhterm += volume_j * wi_dx_gradhterm;
pj->grad_m0_gradhterm += volume_i * wj_dx_gradhterm;
for (int i = 0; i < 3; i++) {
  pi->m1[i] += dx[i] * volume_j * wi_term;
  pj->m1[i] += -dx[i] * volume_i * wj_term;

  pi->grad_m0[i] += volume_j * wi_dx_term[i];
  pj->grad_m0[i] += volume_i * wj_dx_term[i];


  pi->grad_m1_term1_gradhterm[i] += volume_j * dx[i] * wi_dx_gradhterm;
  pj->grad_m1_term1_gradhterm[i] += -volume_i * dx[i] * wj_dx_gradhterm;


  for (int j = 0; j < 3; j++) {
    pi->m2[i][j] += dx[i] * dx[j] * volume_j * wi_term;
    pj->m2[i][j] += dx[i] * dx[j] * volume_i * wj_term;

    pi->grad_m1_term1[i][j] += volume_j * dx[j] * wi_dx_term[i];
    pj->grad_m1_term1[i][j] += -volume_i * dx[j] * wj_dx_term[i];

    pi->grad_m2_term1_gradhterm[i][j] += volume_j * dx[i] * dx[j] * wi_dx_gradhterm;
    pj->grad_m2_term1_gradhterm[i][j] += volume_i * dx[i] * dx[j] * wj_dx_gradhterm;


    for (int k = 0; k < 3; k++) {
        pi->grad_m2_term1[i][j][k] += volume_j * dx[k] * dx[j] * wi_dx_term[i];
        pj->grad_m2_term1[i][j][k] += volume_i * dx[k] * dx[j] * wj_dx_term[i];


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
float wi_dx_gradhterm;
wi_dx_gradhterm = -0.5f *(hydro_dimension * wi + (r / hi) * wi_dx) * hi_inv_dim_plus_one;



pi->m0 += volume_j * wi_term;

pi->grad_m0_gradhterm += volume_j * wi_dx_gradhterm;

for (int i = 0; i < 3; i++) {
  pi->m1[i] += dx[i] * volume_j * wi_term;

  pi->grad_m0[i] += volume_j * wi_dx_term[i];

  pi->grad_m1_term1_gradhterm[i] += volume_j * dx[i] * wi_dx_gradhterm;

  for (int j = 0; j < 3; j++) {
    pi->m2[i][j] += dx[i] * dx[j] * volume_j * wi_term;

    pi->grad_m1_term1[i][j] += volume_j * dx[j] * wi_dx_term[i];

    pi->grad_m2_term1_gradhterm[i][j] += volume_j * dx[i] * dx[j] * wi_dx_gradhterm;

    for (int k = 0; k < 3; k++) {
        pi->grad_m2_term1[i][j][k] += volume_j * dx[k] * dx[j] * wi_dx_term[i];
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


const float h = p->h;
const float h_inv = 1.0f / h;                 /* 1/h */
const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

float volume = p->mass / p->rho_evolved;

int i, j, k;
p->m0 += volume * kernel_root * h_inv_dim;

const float h_inv_dim_plus_one = h_inv_dim * h_inv;
p->grad_m0_gradhterm -= 0.5f * volume * hydro_dimension * kernel_root * h_inv_dim_plus_one;


float grad_m1[3][3];
float grad_m2[3][3][3];
for (i = 0; i < 3; i++) {
    p->grad_m0[i] += p->grad_m0_gradhterm * p->dh_sphgrad[i];
  for (j = 0; j < 3; j++) {
    grad_m1[i][j] = 0.f;
    p->grad_m1_term1[i][j] += p->grad_m1_term1_gradhterm[j] * p->dh_sphgrad[i];
    for (k = 0; k < 3; k++) {
        grad_m2[i][j][k] = 0.f;
        p->grad_m2_term1[i][j][k] += p->grad_m2_term1_gradhterm[j][k] * p->dh_sphgrad[i];
    }
  }
}


// Combine terms to get final m expressions
for (i = 0; i < 3; i++) {
  grad_m1[i][i] +=p->m0;
  for (j = 0; j < 3; j++) {
    grad_m1[i][j] += p->grad_m1_term1[i][j];

    grad_m2[i][i][j] += p->m1[j];
    grad_m2[i][j][i] += p->m1[j];

    for (k = 0; k < 3; k++) {
        grad_m2[i][j][k] += p->grad_m2_term1[i][j][k];
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

  kernel_gradient_i[0] = 0.5f * (Gi[0] + Gj[0]);
  kernel_gradient_i[1] = 0.5f * (Gi[1] + Gj[1]);
  kernel_gradient_i[2] = 0.5f * (Gi[2] + Gj[2]);

  kernel_gradient_j[0] = 0.5f * (Gi[0] + Gj[0]);
  kernel_gradient_j[1] = 0.5f * (Gi[1] + Gj[1]);
  kernel_gradient_j[2] = 0.5f * (Gi[2] + Gj[2]);


  Q_kernel_gradient_i[0] = kernel_gradient_i[0];
  Q_kernel_gradient_i[1] = kernel_gradient_i[1];
  Q_kernel_gradient_i[2] = kernel_gradient_i[2];

  Q_kernel_gradient_j[0] = kernel_gradient_j[0];
  Q_kernel_gradient_j[1] = kernel_gradient_j[1];
  Q_kernel_gradient_j[2] = kernel_gradient_j[2];
}

#endif /* SWIFT_PLANETARY_HYDRO_KERNELS_ETC_H */
