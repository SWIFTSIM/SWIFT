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


    p->m0_no_mean_kernel = 0.f;


    p->grad_m0_no_mean_kernel[0] = 0.f;
  p->grad_m0_no_mean_kernel[1] = 0.f;
  p->grad_m0_no_mean_kernel[2] = 0.f;

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


    pi->m0_no_mean_kernel += pj->mass * wi / pj->rho_evolved;
    pj->m0_no_mean_kernel += pi->mass * wj / pi->rho_evolved;

    for (int i = 0; i < 3; i++) {
      pi->grad_m0_no_mean_kernel[i] += (pj->mass / pj->rho_evolved) * dx[i] * wi_dx * r_inv;
      pj->grad_m0_no_mean_kernel[i] += (pi->mass / pi->rho_evolved) * (-dx[i] * wj_dx * r_inv);
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


    pi->m0_no_mean_kernel += pj->mass * wi / pj->rho_evolved;

        for (int i = 0; i < 3; i++) {
      pi->grad_m0_no_mean_kernel[i] += (pj->mass / pj->rho_evolved) * dx[i] * wi_dx * r_inv;
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


  p->m0_no_mean_kernel += p->mass * kernel_root / p->rho_evolved;
  p->m0_no_mean_kernel *= h_inv_dim;

   for (int i = 0; i < 3; i++) {
      p->grad_m0_no_mean_kernel[i] *= h_inv_dim_plus_one;
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

    zero_sym_matrix(&p->gradient.m2);
    zero_sym_matrix(&p->gradient.grad_m2_term1_x);
    zero_sym_matrix(&p->gradient.grad_m2_term1_y);
    zero_sym_matrix(&p->gradient.grad_m2_term1_z);
    zero_sym_matrix(&p->gradient.grad_m2_term1_gradhterm);


    p->gradient.m0 = 0.f;
    p->gradient.grad_m0_gradhterm = 0.f;
    for (int i = 0; i < 3; i++) {
      p->gradient.m1[i] = 0.f;
      p->gradient.grad_m0[i] = 0.f;
      p->gradient.grad_m1_term1_gradhterm[i] = 0.f;
      for (int j = 0; j < 3; j++) {
        p->gradient.grad_m1_term1[i][j] = 0.f;
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



pi->gradient.m0 += wi_term * volume_j;
pj->gradient.m0 += wj_term * volume_i;

pi->gradient.grad_m0_gradhterm += wi_dx_gradhterm * volume_j;
pj->gradient.grad_m0_gradhterm += wj_dx_gradhterm * volume_i;
for (int i = 0; i < 3; i++) {
  pi->gradient.m1[i] += dx[i] * wi_term * volume_j;
  pj->gradient.m1[i] += -dx[i] * wj_term * volume_i;

  pi->gradient.grad_m0[i] += wi_dx_term[i] * volume_j;
  pj->gradient.grad_m0[i] += wj_dx_term[i] * volume_i;


  pi->gradient.grad_m1_term1_gradhterm[i] += dx[i] * wi_dx_gradhterm * volume_j;
  pj->gradient.grad_m1_term1_gradhterm[i] += -dx[i] * wj_dx_gradhterm * volume_i;

  for (int j = 0; j < 3; j++) {

    pi->gradient.grad_m1_term1[i][j] += dx[j] * wi_dx_term[i] * volume_j;
    pj->gradient.grad_m1_term1[i][j] += -dx[j] * wj_dx_term[i] * volume_i;
  }
}




pi->gradient.m2.xx += dx[0] * dx[0] * wi_term * volume_j;
pi->gradient.m2.yy += dx[1] * dx[1] * wi_term * volume_j;
pi->gradient.m2.zz += dx[2] * dx[2] * wi_term * volume_j;
pi->gradient.m2.xy += dx[0] * dx[1] * wi_term * volume_j;
pi->gradient.m2.xz += dx[0] * dx[2] * wi_term * volume_j;
pi->gradient.m2.yz += dx[1] * dx[2] * wi_term * volume_j;

pj->gradient.m2.xx += dx[0] * dx[0] * wj_term * volume_i;
pj->gradient.m2.yy += dx[1] * dx[1] * wj_term * volume_i;
pj->gradient.m2.zz += dx[2] * dx[2] * wj_term * volume_i;
pj->gradient.m2.xy += dx[0] * dx[1] * wj_term * volume_i;
pj->gradient.m2.xz += dx[0] * dx[2] * wj_term * volume_i;
pj->gradient.m2.yz += dx[1] * dx[2] * wj_term * volume_i;



pi->gradient.grad_m2_term1_gradhterm.xx += dx[0] * dx[0] * wi_dx_gradhterm * volume_j;
pi->gradient.grad_m2_term1_gradhterm.yy += dx[1] * dx[1] * wi_dx_gradhterm * volume_j;
pi->gradient.grad_m2_term1_gradhterm.zz += dx[2] * dx[2] * wi_dx_gradhterm * volume_j;
pi->gradient.grad_m2_term1_gradhterm.xy += dx[0] * dx[1] * wi_dx_gradhterm * volume_j;
pi->gradient.grad_m2_term1_gradhterm.xz += dx[0] * dx[2] * wi_dx_gradhterm * volume_j;
pi->gradient.grad_m2_term1_gradhterm.yz += dx[1] * dx[2] * wi_dx_gradhterm * volume_j;

pj->gradient.grad_m2_term1_gradhterm.xx += dx[0] * dx[0] * wj_dx_gradhterm * volume_i;
pj->gradient.grad_m2_term1_gradhterm.yy += dx[1] * dx[1] * wj_dx_gradhterm * volume_i;
pj->gradient.grad_m2_term1_gradhterm.zz += dx[2] * dx[2] * wj_dx_gradhterm * volume_i;
pj->gradient.grad_m2_term1_gradhterm.xy += dx[0] * dx[1] * wj_dx_gradhterm * volume_i;
pj->gradient.grad_m2_term1_gradhterm.xz += dx[0] * dx[2] * wj_dx_gradhterm * volume_i;
pj->gradient.grad_m2_term1_gradhterm.yz += dx[1] * dx[2] * wj_dx_gradhterm * volume_i;



pi->gradient.grad_m2_term1_x.xx += dx[0] * dx[0] * wi_dx_term[0] * volume_j;
pi->gradient.grad_m2_term1_x.yy += dx[1] * dx[1] * wi_dx_term[0] * volume_j;
pi->gradient.grad_m2_term1_x.zz += dx[2] * dx[2] * wi_dx_term[0] * volume_j;
pi->gradient.grad_m2_term1_x.xy += dx[0] * dx[1] * wi_dx_term[0] * volume_j;
pi->gradient.grad_m2_term1_x.xz += dx[0] * dx[2] * wi_dx_term[0] * volume_j;
pi->gradient.grad_m2_term1_x.yz += dx[1] * dx[2] * wi_dx_term[0] * volume_j;

pj->gradient.grad_m2_term1_x.xx += dx[0] * dx[0] * wj_dx_term[0] * volume_i;
pj->gradient.grad_m2_term1_x.yy += dx[1] * dx[1] * wj_dx_term[0] * volume_i;
pj->gradient.grad_m2_term1_x.zz += dx[2] * dx[2] * wj_dx_term[0] * volume_i;
pj->gradient.grad_m2_term1_x.xy += dx[0] * dx[1] * wj_dx_term[0] * volume_i;
pj->gradient.grad_m2_term1_x.xz += dx[0] * dx[2] * wj_dx_term[0] * volume_i;
pj->gradient.grad_m2_term1_x.yz += dx[1] * dx[2] * wj_dx_term[0] * volume_i;



pi->gradient.grad_m2_term1_y.xx += dx[0] * dx[0] * wi_dx_term[1] * volume_j;
pi->gradient.grad_m2_term1_y.yy += dx[1] * dx[1] * wi_dx_term[1] * volume_j;
pi->gradient.grad_m2_term1_y.zz += dx[2] * dx[2] * wi_dx_term[1] * volume_j;
pi->gradient.grad_m2_term1_y.xy += dx[0] * dx[1] * wi_dx_term[1] * volume_j;
pi->gradient.grad_m2_term1_y.xz += dx[0] * dx[2] * wi_dx_term[1] * volume_j;
pi->gradient.grad_m2_term1_y.yz += dx[1] * dx[2] * wi_dx_term[1] * volume_j;

pj->gradient.grad_m2_term1_y.xx += dx[0] * dx[0] * wj_dx_term[1] * volume_i;
pj->gradient.grad_m2_term1_y.yy += dx[1] * dx[1] * wj_dx_term[1] * volume_i;
pj->gradient.grad_m2_term1_y.zz += dx[2] * dx[2] * wj_dx_term[1] * volume_i;
pj->gradient.grad_m2_term1_y.xy += dx[0] * dx[1] * wj_dx_term[1] * volume_i;
pj->gradient.grad_m2_term1_y.xz += dx[0] * dx[2] * wj_dx_term[1] * volume_i;
pj->gradient.grad_m2_term1_y.yz += dx[1] * dx[2] * wj_dx_term[1] * volume_i;



pi->gradient.grad_m2_term1_z.xx += dx[0] * dx[0] * wi_dx_term[2] * volume_j;
pi->gradient.grad_m2_term1_z.yy += dx[1] * dx[1] * wi_dx_term[2] * volume_j;
pi->gradient.grad_m2_term1_z.zz += dx[2] * dx[2] * wi_dx_term[2] * volume_j;
pi->gradient.grad_m2_term1_z.xy += dx[0] * dx[1] * wi_dx_term[2] * volume_j;
pi->gradient.grad_m2_term1_z.xz += dx[0] * dx[2] * wi_dx_term[2] * volume_j;
pi->gradient.grad_m2_term1_z.yz += dx[1] * dx[2] * wi_dx_term[2] * volume_j;

pj->gradient.grad_m2_term1_z.xx += dx[0] * dx[0] * wj_dx_term[2] * volume_i;
pj->gradient.grad_m2_term1_z.yy += dx[1] * dx[1] * wj_dx_term[2] * volume_i;
pj->gradient.grad_m2_term1_z.zz += dx[2] * dx[2] * wj_dx_term[2] * volume_i;
pj->gradient.grad_m2_term1_z.xy += dx[0] * dx[1] * wj_dx_term[2] * volume_i;
pj->gradient.grad_m2_term1_z.xz += dx[0] * dx[2] * wj_dx_term[2] * volume_i;
pj->gradient.grad_m2_term1_z.yz += dx[1] * dx[2] * wj_dx_term[2] * volume_i;


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



pi->gradient.m0 += wi_term * volume_j;

pi->gradient.grad_m0_gradhterm += wi_dx_gradhterm * volume_j;
for (int i = 0; i < 3; i++) {
  pi->gradient.m1[i] += dx[i] * wi_term * volume_j;

  pi->gradient.grad_m0[i] += wi_dx_term[i] * volume_j;


  pi->gradient.grad_m1_term1_gradhterm[i] += dx[i] * wi_dx_gradhterm * volume_j;

  for (int j = 0; j < 3; j++) {
    pi->gradient.grad_m1_term1[i][j] += dx[j] * wi_dx_term[i] * volume_j;
  }
}




pi->gradient.m2.xx += dx[0] * dx[0] * wi_term * volume_j;
pi->gradient.m2.yy += dx[1] * dx[1] * wi_term * volume_j;
pi->gradient.m2.zz += dx[2] * dx[2] * wi_term * volume_j;
pi->gradient.m2.xy += dx[0] * dx[1] * wi_term * volume_j;
pi->gradient.m2.xz += dx[0] * dx[2] * wi_term * volume_j;
pi->gradient.m2.yz += dx[1] * dx[2] * wi_term * volume_j;



pi->gradient.grad_m2_term1_gradhterm.xx += dx[0] * dx[0] * wi_dx_gradhterm * volume_j;
pi->gradient.grad_m2_term1_gradhterm.yy += dx[1] * dx[1] * wi_dx_gradhterm * volume_j;
pi->gradient.grad_m2_term1_gradhterm.zz += dx[2] * dx[2] * wi_dx_gradhterm * volume_j;
pi->gradient.grad_m2_term1_gradhterm.xy += dx[0] * dx[1] * wi_dx_gradhterm * volume_j;
pi->gradient.grad_m2_term1_gradhterm.xz += dx[0] * dx[2] * wi_dx_gradhterm * volume_j;
pi->gradient.grad_m2_term1_gradhterm.yz += dx[1] * dx[2] * wi_dx_gradhterm * volume_j;


pi->gradient.grad_m2_term1_x.xx += dx[0] * dx[0] * wi_dx_term[0] * volume_j;
pi->gradient.grad_m2_term1_x.yy += dx[1] * dx[1] * wi_dx_term[0] * volume_j;
pi->gradient.grad_m2_term1_x.zz += dx[2] * dx[2] * wi_dx_term[0] * volume_j;
pi->gradient.grad_m2_term1_x.xy += dx[0] * dx[1] * wi_dx_term[0] * volume_j;
pi->gradient.grad_m2_term1_x.xz += dx[0] * dx[2] * wi_dx_term[0] * volume_j;
pi->gradient.grad_m2_term1_x.yz += dx[1] * dx[2] * wi_dx_term[0] * volume_j;



pi->gradient.grad_m2_term1_y.xx += dx[0] * dx[0] * wi_dx_term[1] * volume_j;
pi->gradient.grad_m2_term1_y.yy += dx[1] * dx[1] * wi_dx_term[1] * volume_j;
pi->gradient.grad_m2_term1_y.zz += dx[2] * dx[2] * wi_dx_term[1] * volume_j;
pi->gradient.grad_m2_term1_y.xy += dx[0] * dx[1] * wi_dx_term[1] * volume_j;
pi->gradient.grad_m2_term1_y.xz += dx[0] * dx[2] * wi_dx_term[1] * volume_j;
pi->gradient.grad_m2_term1_y.yz += dx[1] * dx[2] * wi_dx_term[1] * volume_j;



pi->gradient.grad_m2_term1_z.xx += dx[0] * dx[0] * wi_dx_term[2] * volume_j;
pi->gradient.grad_m2_term1_z.yy += dx[1] * dx[1] * wi_dx_term[2] * volume_j;
pi->gradient.grad_m2_term1_z.zz += dx[2] * dx[2] * wi_dx_term[2] * volume_j;
pi->gradient.grad_m2_term1_z.xy += dx[0] * dx[1] * wi_dx_term[2] * volume_j;
pi->gradient.grad_m2_term1_z.xz += dx[0] * dx[2] * wi_dx_term[2] * volume_j;
pi->gradient.grad_m2_term1_z.yz += dx[1] * dx[2] * wi_dx_term[2] * volume_j;

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

int i, j;//, k;
p->gradient.m0 += volume * kernel_root * h_inv_dim;

const float h_inv_dim_plus_one = h_inv_dim * h_inv;
p->gradient.grad_m0_gradhterm -= 0.5f * volume * hydro_dimension * kernel_root * h_inv_dim_plus_one;


for (i = 0; i < 3; i++) {
    p->gradient.grad_m0[i] += p->gradient.grad_m0_gradhterm * p->dh_norm_kernel[i];
  for (j = 0; j < 3; j++) {
    p->gradient.grad_m1_term1[i][j] += p->gradient.grad_m1_term1_gradhterm[j] * p->dh_norm_kernel[i];
  }
}


p->gradient.grad_m2_term1_x.xx += p->gradient.grad_m2_term1_gradhterm.xx * p->dh_norm_kernel[0];
p->gradient.grad_m2_term1_x.yy += p->gradient.grad_m2_term1_gradhterm.yy * p->dh_norm_kernel[0];
p->gradient.grad_m2_term1_x.zz += p->gradient.grad_m2_term1_gradhterm.zz * p->dh_norm_kernel[0];
p->gradient.grad_m2_term1_x.xy += p->gradient.grad_m2_term1_gradhterm.xy * p->dh_norm_kernel[0];
p->gradient.grad_m2_term1_x.xz += p->gradient.grad_m2_term1_gradhterm.xz * p->dh_norm_kernel[0];
p->gradient.grad_m2_term1_x.yz += p->gradient.grad_m2_term1_gradhterm.yz * p->dh_norm_kernel[0];

p->gradient.grad_m2_term1_y.xx += p->gradient.grad_m2_term1_gradhterm.xx * p->dh_norm_kernel[1];
p->gradient.grad_m2_term1_y.yy += p->gradient.grad_m2_term1_gradhterm.yy * p->dh_norm_kernel[1];
p->gradient.grad_m2_term1_y.zz += p->gradient.grad_m2_term1_gradhterm.zz * p->dh_norm_kernel[1];
p->gradient.grad_m2_term1_y.xy += p->gradient.grad_m2_term1_gradhterm.xy * p->dh_norm_kernel[1];
p->gradient.grad_m2_term1_y.xz += p->gradient.grad_m2_term1_gradhterm.xz * p->dh_norm_kernel[1];
p->gradient.grad_m2_term1_y.yz += p->gradient.grad_m2_term1_gradhterm.yz * p->dh_norm_kernel[1];

p->gradient.grad_m2_term1_z.xx += p->gradient.grad_m2_term1_gradhterm.xx * p->dh_norm_kernel[2];
p->gradient.grad_m2_term1_z.yy += p->gradient.grad_m2_term1_gradhterm.yy * p->dh_norm_kernel[2];
p->gradient.grad_m2_term1_z.zz += p->gradient.grad_m2_term1_gradhterm.zz * p->dh_norm_kernel[2];
p->gradient.grad_m2_term1_z.xy += p->gradient.grad_m2_term1_gradhterm.xy * p->dh_norm_kernel[2];
p->gradient.grad_m2_term1_z.xz += p->gradient.grad_m2_term1_gradhterm.xz * p->dh_norm_kernel[2];
p->gradient.grad_m2_term1_z.yz += p->gradient.grad_m2_term1_gradhterm.yz * p->dh_norm_kernel[2];

}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_force_extra_kernel(struct part *restrict p) {

  if (!p->is_h_max) {

    float grad_m1_x[3], grad_m1_y[3], grad_m1_z[3];

    grad_m1_x[0] = p->gradient.grad_m1_term1[0][0] + p->gradient.m0;
    grad_m1_x[1] = p->gradient.grad_m1_term1[0][1];
    grad_m1_x[2] = p->gradient.grad_m1_term1[0][2];

    grad_m1_y[0] = p->gradient.grad_m1_term1[1][0];
    grad_m1_y[1] = p->gradient.grad_m1_term1[1][1] + p->gradient.m0;
    grad_m1_y[2] = p->gradient.grad_m1_term1[1][2];

    grad_m1_z[0] = p->gradient.grad_m1_term1[2][0];
    grad_m1_z[1] = p->gradient.grad_m1_term1[2][1];
    grad_m1_z[2] = p->gradient.grad_m1_term1[2][2] + p->gradient.m0;



    struct sym_matrix grad_m2_x;
    struct sym_matrix grad_m2_y;
    struct sym_matrix grad_m2_z;

    grad_m2_x.xx = p->gradient.grad_m2_term1_x.xx+ 2.f * p->gradient.m1[0];
    grad_m2_x.yy = p->gradient.grad_m2_term1_x.yy;
    grad_m2_x.zz = p->gradient.grad_m2_term1_x.zz;
    grad_m2_x.xy = p->gradient.grad_m2_term1_x.xy+ p->gradient.m1[1];
    grad_m2_x.xz = p->gradient.grad_m2_term1_x.xz+ p->gradient.m1[2];
    grad_m2_x.yz = p->gradient.grad_m2_term1_x.yz;

    grad_m2_y.xx = p->gradient.grad_m2_term1_y.xx;
    grad_m2_y.yy = p->gradient.grad_m2_term1_y.yy + 2.f * p->gradient.m1[1];
    grad_m2_y.zz = p->gradient.grad_m2_term1_y.zz;
    grad_m2_y.xy = p->gradient.grad_m2_term1_y.xy + p->gradient.m1[0];
    grad_m2_y.xz = p->gradient.grad_m2_term1_y.xz;
    grad_m2_y.yz = p->gradient.grad_m2_term1_y.yz + p->gradient.m1[2];

    grad_m2_z.xx = p->gradient.grad_m2_term1_z.xx;
    grad_m2_z.yy = p->gradient.grad_m2_term1_z.yy;
    grad_m2_z.zz = p->gradient.grad_m2_term1_z.zz + 2.f * p->gradient.m1[2];
    grad_m2_z.xy = p->gradient.grad_m2_term1_z.xy;
    grad_m2_z.xz = p->gradient.grad_m2_term1_z.xz + p->gradient.m1[0];
    grad_m2_z.yz = p->gradient.grad_m2_term1_z.yz + p->gradient.m1[1];


    struct sym_matrix m2_inv;
    sym_matrix_invert(&m2_inv, &p->gradient.m2);



    // Vector and symetric matrices needed for muliplication
    float m2_inv_mult_m1[3];
    sym_matrix_multiply_by_vector(m2_inv_mult_m1, &m2_inv, p->gradient.m1);

    float m2_inv_mult_grad_m1_x[3], m2_inv_mult_grad_m1_y[3], m2_inv_mult_grad_m1_z[3];
    sym_matrix_multiply_by_vector(m2_inv_mult_grad_m1_x, &m2_inv, grad_m1_x);
    sym_matrix_multiply_by_vector(m2_inv_mult_grad_m1_y, &m2_inv, grad_m1_y);
    sym_matrix_multiply_by_vector(m2_inv_mult_grad_m1_z, &m2_inv, grad_m1_z);

    struct sym_matrix m2_inv_mult_grad_m2_mult_m2_inv_x;
    sym_matrix_multiplication_ABA(&m2_inv_mult_grad_m2_mult_m2_inv_x, &m2_inv, &grad_m2_x);
    struct sym_matrix m2_inv_mult_grad_m2_mult_m2_inv_y;
    sym_matrix_multiplication_ABA(&m2_inv_mult_grad_m2_mult_m2_inv_y, &m2_inv, &grad_m2_y);
    struct sym_matrix m2_inv_mult_grad_m2_mult_m2_inv_z;
    sym_matrix_multiplication_ABA(&m2_inv_mult_grad_m2_mult_m2_inv_z, &m2_inv, &grad_m2_z);

    float ABA_mult_m1_x[3], ABA_mult_m1_y[3], ABA_mult_m1_z[3];
    sym_matrix_multiply_by_vector(ABA_mult_m1_x, &m2_inv_mult_grad_m2_mult_m2_inv_x, p->gradient.m1);
    sym_matrix_multiply_by_vector(ABA_mult_m1_y, &m2_inv_mult_grad_m2_mult_m2_inv_y, p->gradient.m1);
    sym_matrix_multiply_by_vector(ABA_mult_m1_z, &m2_inv_mult_grad_m2_mult_m2_inv_z, p->gradient.m1);

    float A, B[3], grad_A[3], grad_B[3][3];

    // Calculata A
    A = p->gradient.m0;
    A -= m2_inv_mult_m1[0] * p->gradient.m1[0] + m2_inv_mult_m1[1] * p->gradient.m1[1] + m2_inv_mult_m1[2] * p->gradient.m1[2];
    A = 1 / A;

    // Calculate B
    B[0] = -m2_inv_mult_m1[0];
    B[1] = -m2_inv_mult_m1[1];
    B[2] = -m2_inv_mult_m1[2];


    // Calculate grad A
    grad_A[0] = p->gradient.grad_m0[0];
    grad_A[1] = p->gradient.grad_m0[1];
    grad_A[2] = p->gradient.grad_m0[2];

    grad_A[0] -= 2 * (m2_inv_mult_m1[0] * grad_m1_x[0] +
                         m2_inv_mult_m1[1] * grad_m1_x[1] +
                         m2_inv_mult_m1[2] * grad_m1_x[2]);
    grad_A[1] -= 2 * (m2_inv_mult_m1[0] * grad_m1_y[0] +
                         m2_inv_mult_m1[1] * grad_m1_y[1] +
                         m2_inv_mult_m1[2] * grad_m1_y[2]);
    grad_A[2] -= 2 * (m2_inv_mult_m1[0] * grad_m1_z[0] +
                         m2_inv_mult_m1[1] * grad_m1_z[1] +
                         m2_inv_mult_m1[2] * grad_m1_z[2]);

    grad_A[0] += ABA_mult_m1_x[0] * p->gradient.m1[0] + ABA_mult_m1_x[1] * p->gradient.m1[1] + ABA_mult_m1_x[2] * p->gradient.m1[2];
    grad_A[1] += ABA_mult_m1_y[0] * p->gradient.m1[0] + ABA_mult_m1_y[1] * p->gradient.m1[1] + ABA_mult_m1_y[2] * p->gradient.m1[2];
    grad_A[2] += ABA_mult_m1_z[0] * p->gradient.m1[0] + ABA_mult_m1_z[1] * p->gradient.m1[1] + ABA_mult_m1_z[2] * p->gradient.m1[2];

    grad_A[0] *= -A * A;
    grad_A[1] *= -A * A;
    grad_A[2] *= -A * A;


    // Calculate grad B
    for (int i = 0; i < 3; i++) {
      grad_B[0][i] = -m2_inv_mult_grad_m1_x[i];
      grad_B[1][i] = -m2_inv_mult_grad_m1_y[i];
      grad_B[2][i] = -m2_inv_mult_grad_m1_z[i];

      grad_B[0][i] += ABA_mult_m1_x[i];
      grad_B[1][i] += ABA_mult_m1_y[i];
      grad_B[2][i] += ABA_mult_m1_z[i];
    }


     p->vac_switch = 1.f;

     float x = p->h * sqrtf(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
    float offset = 0.8f;
    if (x > offset){
        float sigma = 0.2f;
        p->vac_switch = expf(-(x - offset) * (x - offset) / (2.f * sigma * sigma));
    }


    p->force.A = A;
    for (int i = 0; i < 3; i++) {
        p->force.B[i] = B[i];
        p->force.grad_A[i] = grad_A[i];
        for (int j = 0; j < 3; j++) {
              p->force.grad_B[i][j] = grad_B[i][j];
          }
    }



    }else{
        // this shouldn't be used because of later if statements
        p->force.A = 1.f;
        p->vac_switch = 1.f;
        for (int i = 0; i < 3; i++) {
            p->force.B[i] = 0.f;
            p->force.grad_A[i] = 0.f;
            for (int j = 0; j < 3; j++) {
                  p->force.grad_B[i][j] = 0.f;
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
        wi_dx_term[i] +=-0.5f *(hydro_dimension * wi + (r / pi->h) * wi_dx) * hid_inv * pi->dh_norm_kernel[i];
        wj_dx_term[i] += -0.5f *(hydro_dimension * wj + (r / pj->h) * wj_dx) * hjd_inv * pj->dh_norm_kernel[i];
    }



      for (i = 0; i < 3; i++) {
        modified_grad_wi[i] = pi->force.A * wi_dx_term[i] + pi->force.grad_A[i] * wi_term + pi->force.A * pi->force.B[i] * wi_term;
        modified_grad_wj[i] = pj->force.A * wj_dx_term[i] + pj->force.grad_A[i] * wj_term + pj->force.A * pj->force.B[i] * wj_term;
      }

      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
        modified_grad_wi[i] += pi->force.A * pi->force.B[j] * dx[j] * wi_dx_term[i];
        modified_grad_wi[i] += pi->force.grad_A[i] * pi->force.B[j] * dx[j] * wi_term;
        modified_grad_wi[i] += pi->force.A * pi->force.grad_B[i][j] * dx[j] * wi_term;

        modified_grad_wj[i] += -pj->force.A * pj->force.B[j] * dx[j] * wj_dx_term[i];
         modified_grad_wj[i] += -pj->force.grad_A[i] *  pj->force.B[j] * dx[j] * wj_term;
          modified_grad_wj[i] += -pj->force.A * pj->force.grad_B[i][j] * dx[j]  * wj_term;
      }
      }



          float modified_wi = pi->force.A * wi_term;
      float modified_wj = pj->force.A * wj_term;
      modified_wi += pi->force.A * pi->force.B[0] * dx[0] * wi_term + pi->force.A * pi->force.B[1] * dx[1] * wi_term + pi->force.A * pi->force.B[2] * dx[2] * wi_term;
      modified_wj += -(pj->force.A * pj->force.B[0] * dx[0] * wj_term + pj->force.A * pj->force.B[1] * dx[1] * wj_term + pj->force.A * pj->force.B[2] * dx[2] * wj_term);

      float wi_dx_term_vac[3], wj_dx_term_vac[3];
        for (i = 0; i < 3; i++) {
            wi_dx_term_vac[i] = dx[i] * r_inv * wi_dx * hid_inv;
            wj_dx_term_vac[i] = -dx[i] * r_inv * wj_dx * hjd_inv;
        }

        // h term
      for (i = 0; i < 3; i++) {
        wi_dx_term_vac[i] += -(hydro_dimension * wi + (r / pi->h) * wi_dx) * hid_inv * pi->dh_norm_kernel[i];
        wj_dx_term_vac[i] += -(hydro_dimension * wj + (r / pj->h) * wj_dx) * hjd_inv * pj->dh_norm_kernel[i];
    }

        for (i = 0; i < 3; i++) {

        modified_grad_wi[i] *= pi->vac_switch;
        modified_grad_wj[i] *= pj->vac_switch;

        modified_grad_wi[i] += wi_dx_term_vac[i];
        modified_grad_wj[i] += wj_dx_term_vac[i];

        modified_grad_wi[i] -= pi->vac_switch * wi_dx_term_vac[i];
        modified_grad_wj[i] -= pj->vac_switch * wj_dx_term_vac[i];


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
