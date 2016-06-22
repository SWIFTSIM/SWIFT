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

#include <math.h>

/* Includes. */
#include "const.h"
#include "error.h"
#include "inline.h"
#include "vector.h"

/* ------------------------------------------------------------------------- */
#if defined(CUBIC_SPLINE_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Cubic spline (M4)"
#define kernel_degree 3 /* Degree of the polynomial */
#define kernel_ivals 2  /* Number of branches */
#define kernel_gamma ((float)(1.825742))
#define kernel_constant ((float)(16. * M_1_PI))
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {3.f,  -3.f, 0.f,  0.5f, /* 0 < u < 0.5 */
                                    -1.f, 3.f,  -3.f, 1.f,  /* 0.5 < u < 1 */
                                    0.f,  0.f,  0.f,  0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(QUARTIC_SPLINE_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Quartic spline (M5)"
#define kernel_degree 4
#define kernel_ivals 5
#define kernel_gamma ((float)(2.018932))
#define kernel_constant ((float)(15625. * M_1_PI / 512.))
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
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
#define kernel_gamma ((float)(2.195775))
#define kernel_constant ((float)(2187. * M_1_PI / 40.))
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
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
#define kernel_gamma ((float)(1.936492))
#define kernel_constant ((float)(21. * M_1_PI / 2.))
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        4.f, -15.f, 20.f, -10.f, 0.f, 1.f,  /* 0 < u < 1 */
        0.f, 0.f,   0.f,  0.f,   0.f, 0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(WENDLAND_C4_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Wendland C4"
#define kernel_degree 8
#define kernel_ivals 1
#define kernel_gamma ((float)(2.207940))
#define kernel_constant ((float)(495. * M_1_PI / 32.))
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
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
#define kernel_gamma ((float)(2.449490))
#define kernel_constant ((float)(1365. * M_1_PI / 64.))
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
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
#define kernel_gamma2 ((float)(kernel_gamma * kernel_gamma))
#define kernel_gamma3 ((float)(kernel_gamma * kernel_gamma * kernel_gamma))
#define kernel_gamma4 \
  ((float)(kernel_gamma * kernel_gamma * kernel_gamma * kernel_gamma))
#define kernel_igamma ((float)(1. / kernel_gamma))
#define kernel_igamma2 ((float)(kernel_igamma * kernel_igamma))
#define kernel_igamma3 ((float)(kernel_igamma * kernel_igamma * kernel_igamma))
#define kernel_igamma4 \
  ((float)(kernel_igamma * kernel_igamma * kernel_igamma * kernel_igamma))

/* The number of branches */
#define kernel_ivals_f ((float)(kernel_ivals))

/* Kernel self contribution (i.e. W(0,h)) */
#define kernel_root \
  ((float)(kernel_coeffs[kernel_degree]) * kernel_constant * kernel_igamma3)

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
    float u, float *restrict W, float *restrict dW_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  const float x = u * kernel_igamma;

#if kernel_ivals == 1
  /* Only one branch in this case */
  const float *const coeffs = &kernel_coeffs[0];
#else
  /* Pick the correct branch of the kernel */
  const int temp = (int)(x * kernel_ivals_f);
  const int ind = temp > kernel_ivals ? kernel_ivals : temp;
  const float *const coeffs = &kernel_coeffs[ind * (kernel_degree + 1)];
#endif

  /* First two terms of the polynomial ... */
  float w = coeffs[0] * x + coeffs[1];
  float dw_dx = coeffs[0];

  /* ... and the rest of them */
  for (int k = 2; k <= kernel_degree; k++) {
    dw_dx = dw_dx * x + w;
    w = x * w + coeffs[k];
  }

  /* Return everything */
  *W = w * kernel_constant * kernel_igamma3;
  *dW_dx = dw_dx * kernel_constant * kernel_igamma4;
}

/**
 * @brief Computes the kernel function.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param W (return) The value of the kernel function $W(x,h)$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval(
    float u, float *restrict W) {
  /* Go to the range [0,1[ from [0,H[ */
  const float x = u * kernel_igamma;

#if kernel_ivals == 1
  /* Only one branch in this case */
  const float *const coeffs = &kernel_coeffs[0];
#else
  /* Pick the correct branch of the kernel */
  const int temp = (int)(x * kernel_ivals_f);
  const int ind = temp > kernel_ivals ? kernel_ivals : temp;
  const float *const coeffs = &kernel_coeffs[ind * (kernel_degree + 1)];
#endif

  /* First two terms of the polynomial ... */
  float w = coeffs[0] * x + coeffs[1];

  /* ... and the rest of them */
  for (int k = 2; k <= kernel_degree; k++) w = x * w + coeffs[k];

  /* Return everything */
  *W = w * kernel_constant * kernel_igamma3;
}

#ifdef VECTORIZE

static const vector kernel_igamma_vec = FILL_VEC((float)kernel_igamma);

static const vector kernel_ivals_vec = FILL_VEC((float)kernel_ivals);

static const vector kernel_constant_vec = FILL_VEC((float)kernel_constant);

static const vector kernel_igamma3_vec = FILL_VEC((float)kernel_igamma3);

static const vector kernel_igamma4_vec = FILL_VEC((float)kernel_igamma4);

/**
 * @brief Computes the kernel function and its derivative (Vectorised version).
 *
 * Return 0 if $u > \\gamma = H/h$
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param w (return) The value of the kernel function $W(x,h)$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 */
__attribute__((always_inline)) INLINE static void kernel_deval_vec(
    vector *u, vector *w, vector *dw_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x;
  x.v = u->v * kernel_igamma_vec.v;

  /* Load x and get the interval id. */
  vector ind;
  ind.m = vec_ftoi(vec_fmin(x.v * kernel_ivals_vec.v, kernel_ivals_vec.v));

  /* load the coefficients. */
  vector c[kernel_degree + 1];
  for (int k = 0; k < VEC_SIZE; k++)
    for (int j = 0; j < kernel_degree + 1; j++)
      c[j].f[k] = kernel_coeffs[ind.i[k] * (kernel_degree + 1) + j];

  /* Init the iteration for Horner's scheme. */
  w->v = (c[0].v * x.v) + c[1].v;
  dw_dx->v = c[0].v;

  /* And we're off! */
  for (int k = 2; k <= kernel_degree; k++) {
    dw_dx->v = (dw_dx->v * x.v) + w->v;
    w->v = (x.v * w->v) + c[k].v;
  }

  /* Return everything */
  w->v = w->v * kernel_constant_vec.v * kernel_igamma3_vec.v;
  dw_dx->v = dw_dx->v * kernel_constant_vec.v * kernel_igamma4_vec.v;
}

#endif

/* Some cross-check functions */
void hydro_kernel_dump(int N);

#endif  // SWIFT_KERNEL_HYDRO_H
