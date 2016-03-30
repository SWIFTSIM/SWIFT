/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_KERNEL_HYDRO_H
#define SWIFT_KERNEL_HYDRO_H

/* Includes. */
#include "const.h"
#include "error.h"
#include "inline.h"
#include "vector.h"

/* ------------------------------------------------------------------------- */
#if defined(CUBIC_SPLINE_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Cubic spline (M4)"
#define kernel_degree 3
#define kernel_ivals 2
#define kernel_gamma 1.825742
#define kernel_constant 16. * M_1_PI
static float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {3.f,  -3.f, 0.f,  0.5f, /* 0 < u < 0.5 */
                                    -1.f, 3.f,  -3.f, 1.f,  /* 0.5 < u < 1 */
                                    0.f,  0.f,  0.f,  0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(QUARTIC_SPLINE_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Quartic spline (M5)"
#define kernel_degree 4
#define kernel_ivals 5
#define kernel_gamma 2.018932
#define kernel_constant 15625. * M_1_PI / 512.
static float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        6.f,  0.f,  -2.4f, 0.f,   0.368f, /* 0 < u < 0.2 */
        -4.f, 8.f,  -4.8f, 0.32f, 0.352f, /* 0.2 < u < 0.4 */
        -4.f, 8.f,  -4.8f, 0.32f, 0.352f, /* 0.4 < u < 0.6 */
        1.f,  -4.f, 6.f,   -4.f,  1.f,    /* 0.6 < u < 0.8 */
        1.f,  -4.f, 6.f,   -4.f,  1.f,    /* 0.8 < u < 1 */
        0.f,  0.f,  0.f,   0.f,   0.f};   /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(QUINTIC_SPLINE_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Quintic spline (M6)"
#define kernel_degree 5
#define kernel_ivals 3
#define kernel_gamma 2.195775
#define kernel_constant 2187. * M_1_PI / 40.
static float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        -10.f,        10.f,      0.f,
        -2.2222222f,  0.f,       0.271604938f, /* 0 < u < 1/3 */
        5.f,          -15.f,     16.666667f,
        -7.77777777f, 0.925925f, 0.209876543f, /* 1/3 < u < 2/3 */
        -1.f,         5.f,       -10.f,
        10.f,         -5.f,      1.f, /* 2/3 < u < 1. */
        0.f,          0.f,       0.f,
        0.f,          0.f,       0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(WENDLAND_C2_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Wendland C2"
#define kernel_degree 5
#define kernel_ivals 1
#define kernel_gamma 1.936492
#define kernel_constant 21. * M_1_PI / 2.
static float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        4.f, -15.f, 20.f, -10.f, 1.f,  /* 0 < u < 1 */
        0.f, 0.f,   0.f,  0.f,   0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(WENDLAND_C4_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Wendland C4"
#define kernel_degree 8
#define kernel_ivals 1
#define kernel_gamma 2.207940
#define kernel_constant 495. * M_1_PI / 32.
static float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        11.666667f, -64.f,       140.f, -149.333333f, 70.f,
        0.f,        -9.3333333f, 0.f,   1.f, /* 0 < u < 1 */
        0.f,        0.f,         0.f,   0.f,          0.f,
        0.f,        0.f,         0.f,   0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(WENDLAND_C6_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Wendland C6"
#define kernel_degree 11
#define kernel_ivals 1
#define kernel_gamma 2.449490
#define kernel_constant 1365. * M_1_PI / 64.
static float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        32.f, -231.f, 704.f, -1155.f, 1056.f, -462.f,
        0.f,  66.f,   0.f,   -11.f,   0.f,    1.f, /* 0 < u < 1 */
        0.f,  0.f,    0.f,   0.f,     0.f,    0.f,
        0.f,  0.f,    0.f,   0.f,     0.f,    0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#else

#error "A kernel function must be chosen in const.h !!"

/* ------------------------------------------------------------------------- */
#endif

/* Ok, now comes the real deal. */

/* First some powers of gamma = H/h */
#define kernel_gamma2 kernel_gamma *kernel_gamma
#define kernel_gamma3 kernel_gamma2 *kernel_gamma
#define kernel_gamma4 kernel_gamma3 *kernel_gamma
#define kernel_igamma 1. / kernel_gamma
#define kernel_igamma2 kernel_igamma *kernel_igamma
#define kernel_igamma3 kernel_igamma2 *kernel_igamma
#define kernel_igamma4 kernel_igamma3 *kernel_igamma

/* Some powers of eta */
#define kernel_eta3 const_eta_kernel *const_eta_kernel *const_eta_kernel

/* The number of neighbours (i.e. N_ngb) */
#define kernel_nwneigh 4.0 * M_PI *kernel_gamma3 *kernel_eta3 / 3.0

/* Kernel self contribution (i.e. W(0,h)) */
#define kernel_root \
  (kernel_coeffs[kernel_degree]) * kernel_constant *kernel_igamma3

/**
 * @brief Computes the kernel function and its derivative.
 *
 * Return 0 if $u > \\gamma = H/h$
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param W (return) The value of the kernel function $W(x,h)$.
 * @param dW_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 */
__attribute__((always_inline)) INLINE static void kernel_deval(
    float u, float *const W, float *const dW_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  const float x = u * (float)kernel_igamma;

  /* Pick the correct branch of the kernel */
  const int ind = (int)fminf(x * (float)kernel_ivals, kernel_ivals);
  const float *const coeffs = &kernel_coeffs[ind * (kernel_degree + 1)];

  /* First two terms of the polynomial ... */
  float w = coeffs[0] * x + coeffs[1];
  float dw_dx = coeffs[0];

  /* ... and the rest of them */
  for (int k = 2; k <= kernel_degree; k++) {
    dw_dx = dw_dx * x + w;
    w = x * w + coeffs[k];
  }

  /* Return everything */
  *W = w * (float)kernel_constant * (float)kernel_igamma3;
  *dW_dx = dw_dx * (float)kernel_constant * (float)kernel_igamma4;
}

/**
 * @brief Computes the kernel function.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param W (return) The value of the kernel function $W(x,h)$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval(float x,
                                                              float *const W) {
  const int ind = fminf(x * 0.5f, kernel_ivals);
  const float *const coeffs = &kernel_coeffs[ind * (kernel_degree + 1)];
  float w = coeffs[0] * x + coeffs[1];
  for (int k = 2; k <= kernel_degree; k++) w = x * w + coeffs[k];
  *W = w;
}

/* Some cross-check functions */
void hydro_kernel_dump(int N);

#endif  // SWIFT_KERNEL_HYDRO_H
