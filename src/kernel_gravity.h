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
#ifndef SWIFT_KERNEL_GRAVITY_H
#define SWIFT_KERNEL_GRAVITY_H

/* Includes. */
#include "const.h"
#include "inline.h"
#include "vector.h"

/* #define const_epsilon 1. */
/* #define const_iepsilon (1. / const_epsilon) */
/* #define const_iepsilon2 (const_iepsilon *const_iepsilon) */
/* #define const_iepsilon3 (const_iepsilon2 *const_iepsilon) */
/* #define const_iepsilon4 (const_iepsilon2 *const_iepsilon2) */
/* #define const_iepsilon5 (const_iepsilon3 *const_iepsilon2) */
/* #define const_iepsilon6 (const_iepsilon3 *const_iepsilon3) */

/* The gravity kernel is defined as a degree 6 polynomial in the distance
   r. The resulting value should be post-multiplied with r^-3, resulting
   in a polynomial with terms ranging from r^-3 to r^3, which are
   sufficient to model both the direct potential as well as the splines
   near the origin. */


/* Coefficients for the kernel. */
#define kernel_grav_name "Gadget-2 softening kernel"
#define kernel_grav_degree 6 /* Degree of the polynomial */
#define kernel_grav_ivals 2  /* Number of branches */
static const float kernel_grav_coeffs[(kernel_grav_degree + 1) * (kernel_grav_ivals + 1)]
    __attribute__((aligned(16))) = {
  32.f,  -38.4f, 0.f,  10.66666667f, 0.f, 0.f, 0.f,/* 0 < u < 0.5 */
  -10.66666667f, 38.4f,  -48.f, 21.3333333, 0.f, 0.f, 0.66666667f,  /* 0.5 < u < 1 */
  0.f,  0.f,  0.f,  0.f, 0.f, 0.f, 0.f}; /* 1 < u */


/**
 * @brief Computes the kernel function.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param W (return) The value of the kernel function $W(x,h)$.
 */
__attribute__((always_inline)) INLINE static void kernel_grav_eval(float u,
                                                              float *const W) {
  /* Go to the range [0,1[ from [0,H[ */
  const float x = u * (float)kernel_grav_igamma;

  /* Pick the correct branch of the kernel */
  const int ind = (int)fminf(x * (float)kernel_grav_ivals, kernel_grav_ivals);
  const float *const coeffs = &kernel_grav_coeffs[ind * (kernel_grav_degree + 1)];

  /* First two terms of the polynomial ... */
  float w = coeffs[0] * x + coeffs[1];

  /* ... and the rest of them */
  for (int k = 2; k <= kernel_grav_degree; k++) w = x * w + coeffs[k];

  /* Return everything */
  *W = w * (float)kernel_grav_igamma3;
}





#if 0


/* Coefficients for the gravity kernel. */
#define kernel_grav_degree 6
#define kernel_grav_ivals 2
#define kernel_grav_scale (2 * const_iepsilon)
static float kernel_grav_coeffs
    [(kernel_grav_degree + 1) * (kernel_grav_ivals + 1)] = {
        32.0f * const_iepsilon6,         -192.0f / 5.0f * const_iepsilon5,
        0.0f,                            32.0f / 3.0f * const_iepsilon3,
        0.0f,                            0.0f,
        0.0f,                            -32.0f / 3.0f * const_iepsilon6,
        192.0f / 5.0f * const_iepsilon5, -48.0f * const_iepsilon4,
        64.0f / 3.0f * const_iepsilon3,  0.0f,
        0.0f,                            -1.0f / 15.0f,
        0.0f,                            0.0f,
        0.0f,                            0.0f,
        0.0f,                            0.0f,
        1.0f};

/**
 * @brief Computes the gravity cubic spline for a given distance x.
 */

__attribute__((always_inline)) INLINE static void kernel_grav_eval(float x,
                                                                   float *W) {
  int ind = fmin(x * kernel_grav_scale, kernel_grav_ivals);
  float *coeffs = &kernel_grav_coeffs[ind * (kernel_grav_degree + 1)];
  float w = coeffs[0] * x + coeffs[1];
  for (int k = 2; k <= kernel_grav_degree; k++) w = x * w + coeffs[k];
  *W = w;
}

#ifdef VECTORIZE

/**
 * @brief Computes the gravity cubic spline for a given distance x (Vectorized
 * version).
 */

__attribute__((always_inline))
    INLINE static void kernel_grav_eval_vec(vector *x, vector *w) {

  vector ind, c[kernel_grav_degree + 1];
  int j, k;

  /* Load x and get the interval id. */
  ind.m = vec_ftoi(vec_fmin(x->v * vec_set1(kernel_grav_scale),
                            vec_set1((float)kernel_grav_ivals)));

  /* load the coefficients. */
  for (k = 0; k < VEC_SIZE; k++)
    for (j = 0; j < kernel_grav_degree + 1; j++)
      c[j].f[k] = kernel_grav_coeffs[ind.i[k] * (kernel_grav_degree + 1) + j];

  /* Init the iteration for Horner's scheme. */
  w->v = (c[0].v * x->v) + c[1].v;

  /* And we're off! */
  for (int k = 2; k <= kernel_grav_degree; k++) w->v = (x->v * w->v) + c[k].v;
}

#endif

/* Blending function stuff
 * --------------------------------------------------------------------------------------------
 */

/* Coefficients for the blending function. */
#define blender_degree 3
#define blender_ivals 3
#define blender_scale 4.0f
static float blender_coeffs[(blender_degree + 1) * (blender_ivals + 1)] = {
    0.0f,   0.0f,  0.0f,   1.0f,  -32.0f, 24.0f, -6.0f, 1.5f,
    -32.0f, 72.0f, -54.0f, 13.5f, 0.0f,   0.0f,  0.0f,  0.0f};

/**
 * @brief Computes the cubic spline blender for a given distance x.
 */

__attribute__((always_inline)) INLINE static void blender_eval(float x,
                                                               float *W) {
  int ind = fmin(x * blender_scale, blender_ivals);
  float *coeffs = &blender_coeffs[ind * (blender_degree + 1)];
  float w = coeffs[0] * x + coeffs[1];
  for (int k = 2; k <= blender_degree; k++) w = x * w + coeffs[k];
  *W = w;
}

/**
 * @brief Computes the cubic spline blender and its derivative for a given
 * distance x.
 */

__attribute__((always_inline)) INLINE static void blender_deval(float x,
                                                                float *W,
                                                                float *dW_dx) {
  int ind = fminf(x * blender_scale, blender_ivals);
  float *coeffs = &blender_coeffs[ind * (blender_degree + 1)];
  float w = coeffs[0] * x + coeffs[1];
  float dw_dx = coeffs[0];
  for (int k = 2; k <= blender_degree; k++) {
    dw_dx = dw_dx * x + w;
    w = x * w + coeffs[k];
  }
  *W = w;
  *dW_dx = dw_dx;
}

#ifdef VECTORIZE

/**
 * @brief Computes the cubic spline blender and its derivative for a given
 * distance x (Vectorized version). Gives a sensible answer only if x<2.
 */

__attribute__((always_inline)) INLINE static void blender_eval_vec(vector *x,
                                                                   vector *w) {

  vector ind, c[blender_degree + 1];
  int j, k;

  /* Load x and get the interval id. */
  ind.m = vec_ftoi(
      vec_fmin(x->v * vec_set1(blender_scale), vec_set1((float)blender_ivals)));

  /* load the coefficients. */
  for (k = 0; k < VEC_SIZE; k++)
    for (j = 0; j < blender_degree + 1; j++)
      c[j].f[k] = blender_coeffs[ind.i[k] * (blender_degree + 1) + j];

  /* Init the iteration for Horner's scheme. */
  w->v = (c[0].v * x->v) + c[1].v;

  /* And we're off! */
  for (int k = 2; k <= blender_degree; k++) w->v = (x->v * w->v) + c[k].v;
}

/**
 * @brief Computes the cubic spline blender and its derivative for a given
 * distance x (Vectorized version). Gives a sensible answer only if x<2.
 */

__attribute__((always_inline))
    INLINE static void blender_deval_vec(vector *x, vector *w, vector *dw_dx) {

  vector ind, c[blender_degree + 1];
  int j, k;

  /* Load x and get the interval id. */
  ind.m = vec_ftoi(
      vec_fmin(x->v * vec_set1(blender_scale), vec_set1((float)blender_ivals)));

  /* load the coefficients. */
  for (k = 0; k < VEC_SIZE; k++)
    for (j = 0; j < blender_degree + 1; j++)
      c[j].f[k] = blender_coeffs[ind.i[k] * (blender_degree + 1) + j];

  /* Init the iteration for Horner's scheme. */
  w->v = (c[0].v * x->v) + c[1].v;
  dw_dx->v = c[0].v;

  /* And we're off! */
  for (int k = 2; k <= blender_degree; k++) {
    dw_dx->v = (dw_dx->v * x->v) + w->v;
    w->v = (x->v * w->v) + c[k].v;
  }
}

#endif
#endif

void gravity_kernel_dump(float r_max, int N);

#endif  // SWIFT_KERNEL_GRAVITY_H
