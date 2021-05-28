/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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

#ifndef SWIFT_KERNEL_DARK_MATTER_H
#define SWIFT_KERNEL_DARK_MATTER_H

/**
 * @file kernel_dark_matter.h
 * @brief Kernel functions from SPH (scalar and vector version).
 *
 * All constants and kernel coefficients are taken from table 1 of
 * Dehnen & Aly, MNRAS, 425, pp. 1062-1082 (2012).
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "error.h"
#include "inline.h"
#include "minmax.h"

/* ------------------------------------------------------------------------- */
/* Coefficients for the kernel. */
#define dm_kernel_name "Cubic spline"
#define dm_kernel_degree 3 /*!< Degree of the polynomial */
#define dm_kernel_ivals 2  /*!< Number of branches */

/* For 3D dimension */
#define dm_kernel_gamma ((float)(1.825742))
#define dm_kernel_constant ((float)(16. * M_1_PI))

static const float dm_kernel_coeffs[(dm_kernel_degree + 1) * (dm_kernel_ivals + 1)]
    __attribute__((aligned(16))) = {3.f,  -3.f, 0.f,  0.5f, /* 0 < u < 0.5 */
                                    -1.f, 3.f,  -3.f, 1.f,  /* 0.5 < u < 1 */
                                    0.f,  0.f,  0.f,  0.f}; /* 1 < u */

/* First some powers of gamma = H/h */
#define dm_kernel_gamma_inv ((float)(1. / dm_kernel_gamma))
#define dm_kernel_gamma2 ((float)(dm_kernel_gamma * dm_kernel_gamma))
#define dm_kernel_gamma3 ((float)(dm_kernel_gamma * dm_kernel_gamma * dm_kernel_gamma))

/* define gamma^d, gamma^(d+1), 1/gamma^d and 1/gamma^(d+1) */
#define dm_kernel_gamma_dim ((float)(dm_kernel_gamma * dm_kernel_gamma * dm_kernel_gamma))
#define dm_kernel_gamma_dim_plus_one \
((float)(dm_kernel_gamma * dm_kernel_gamma * dm_kernel_gamma * dm_kernel_gamma))
#define dm_kernel_gamma_inv_dim \
((float)(1. / (dm_kernel_gamma * dm_kernel_gamma * dm_kernel_gamma)))
#define dm_kernel_gamma_inv_dim_plus_one \
((float)(1. / (dm_kernel_gamma * dm_kernel_gamma * dm_kernel_gamma * dm_kernel_gamma)))

/* The number of branches (floating point conversion) */
#define dm_kernel_ivals_f ((float)(dm_kernel_ivals))

/* Kernel self contribution (i.e. W(0,h)) */
#define dm_kernel_root                                          \
((float)(dm_kernel_coeffs[dm_kernel_degree]) * dm_kernel_constant * \
dm_kernel_gamma_inv_dim)


/* Kernel normalisation constant (volume term) */
#define dimension_unit_sphere ((float)(4. * M_PI / 3.))
#define dm_kernel_norm ((float)(dm_dimension_unit_sphere * dm_kernel_gamma_dim))

/* ------------------------------------------------------------------------- */

/**
 * @brief Computes the kernel function and its derivative.
 *
 * The kernel function needs to be mutliplied by \f$h^{-d}\f$ and the gradient
 * by \f$h^{-(d+1)}\f$, where \f$d\f$ is the dimensionality of the problem.
 *
 * Returns 0 if \f$u > \gamma = H/h\f$.
 *
 * @param u The ratio of the distance to the smoothing length \f$u = x/h\f$.
 * @param W (return) The value of the kernel function \f$W(x,h)\f$.
 * @param dW_dx (return) The norm of the gradient of \f$|\nabla W(x,h)|\f$.
 */
__attribute__((always_inline)) INLINE static void dm_kernel_deval(
                                                               float u, float *restrict W, float *restrict dW_dx) {
    
    /* Go to the range [0,1[ from [0,H[ */
    const float x = u * dm_kernel_gamma_inv;
    
    /* Pick the correct branch of the kernel */
    const int temp = (int)(x * dm_kernel_ivals_f);
    const int ind = temp > dm_kernel_ivals ? dm_kernel_ivals : temp;
    const float *const coeffs = &dm_kernel_coeffs[ind * (dm_kernel_degree + 1)];
    
    /* First two terms of the polynomial ... */
    float w = coeffs[0] * x + coeffs[1];
    float dw_dx = coeffs[0];
    
    /* ... and the rest of them */
    for (int k = 2; k <= dm_kernel_degree; k++) {
        dw_dx = dw_dx * x + w;
        w = x * w + coeffs[k];
    }
    
    w = max(w, 0.f);
    dw_dx = min(dw_dx, 0.f);
    
    /* Return everything */
    *W = w * dm_kernel_constant * dm_kernel_gamma_inv_dim;
    *dW_dx = dw_dx * dm_kernel_constant * dm_kernel_gamma_inv_dim_plus_one;
}

/**
 * @brief Computes the kernel function.
 *
 * The kernel function needs to be mutliplied by \f$h^{-d}\f$,
 * where \f$d\f$ is the dimensionality of the problem.
 *
 * Returns 0 if \f$u > \gamma = H/h\f$
 *
 * @param u The ratio of the distance to the smoothing length \f$u = x/h\f$.
 * @param W (return) The value of the kernel function \f$W(x,h)\f$.
 */
__attribute__((always_inline)) INLINE static void dm_kernel_eval(
                                                              float u, float *restrict W) {
    
    /* Go to the range [0,1[ from [0,H[ */
    const float x = u * dm_kernel_gamma_inv;
    
    /* Pick the correct branch of the kernel */
    const int temp = (int)(x * dm_kernel_ivals_f);
    const int ind = temp > dm_kernel_ivals ? dm_kernel_ivals : temp;
    const float *const coeffs = &dm_kernel_coeffs[ind * (dm_kernel_degree + 1)];
    
    /* First two terms of the polynomial ... */
    float w = coeffs[0] * x + coeffs[1];
    
    /* ... and the rest of them */
    for (int k = 2; k <= dm_kernel_degree; k++) w = x * w + coeffs[k];
    
    w = max(w, 0.f);
    
    /* Return everything */
    *W = w * dm_kernel_constant * dm_kernel_gamma_inv_dim;
}



#endif  // SWIFT_KERNEL_DARK_MATTER_H
