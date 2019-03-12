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

/**
 * @file kernel_hydro.h
 * @brief Kernel functions for SPH (scalar and vector version).
 *
 * All constants and kernel coefficients are taken from table 1 of
 * Dehnen & Aly, MNRAS, 425, pp. 1062-1082 (2012).
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "dimension.h"
#include "error.h"
#include "inline.h"
#include "minmax.h"
#include "vector.h"

/* ------------------------------------------------------------------------- */
#if defined(CUBIC_SPLINE_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Cubic spline (M4)"
#define kernel_degree 3 /*!< Degree of the polynomial */
#define kernel_ivals 2  /*!< Number of branches */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(1.825742))
#define kernel_constant ((float)(16. * M_1_PI))
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(1.778002))
#define kernel_constant ((float)(80. * M_1_PI / 7.))
#elif defined(HYDRO_DIMENSION_1D)
#define kernel_gamma ((float)(1.732051))
#define kernel_constant ((float)(8. / 3.))
#endif
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {3.f,  -3.f, 0.f,  0.5f, /* 0 < u < 0.5 */
                                    -1.f, 3.f,  -3.f, 1.f,  /* 0.5 < u < 1 */
                                    0.f,  0.f,  0.f,  0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(QUARTIC_SPLINE_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Quartic spline (M5)"
#define kernel_degree 4 /* Degree of the polynomial */
#define kernel_ivals 5  /* Number of branches */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(2.018932))
#define kernel_constant ((float)(15625. * M_1_PI / 512.))
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(1.977173))
#define kernel_constant ((float)(46875. * M_1_PI / 2398.))
#elif defined(HYDRO_DIMENSION_1D)
#define kernel_gamma ((float)(1.936492))
#define kernel_constant ((float)(3125. / 768.))
#endif
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
#define kernel_degree 5 /* Degree of the polynomial */
#define kernel_ivals 3  /* Number of branches */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(2.195775))
#define kernel_constant ((float)(2187. * M_1_PI / 40.))
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(2.158131))
#define kernel_constant ((float)(15309. * M_1_PI / 478.))
#elif defined(HYDRO_DIMENSION_1D)
#define kernel_gamma ((float)(2.121321))
#define kernel_constant ((float)(243. / 40.))
#endif
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
#define kernel_degree 5 /* Degree of the polynomial */
#define kernel_ivals 1  /* Number of branches */
#if defined(HYDRO_DIMENSION_1D)
/* Wendland C* have different form in 1D than 2D/3D */
#define kernel_gamma ((float)(1.620185))
#define kernel_constant ((float)(5. / 4.))
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        0.f, -3.f, 8.f, -6.f, 0.f, 1.f, /* 0 < u < 1 */
        0.f, 0.f,  0.f, 0.f,  0.f, 0.f};
#else
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(1.936492))
#define kernel_constant ((float)(21. * M_1_PI / 2.))
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(1.897367))
#define kernel_constant ((float)(7. * M_1_PI))
#endif
static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        4.f, -15.f, 20.f, -10.f, 0.f, 1.f,  /* 0 < u < 1 */
        0.f, 0.f,   0.f,  0.f,   0.f, 0.f}; /* 1 < u */
#endif
/* ------------------------------------------------------------------------- */
#elif defined(WENDLAND_C4_KERNEL)

/* Coefficients for the kernel. */
#define kernel_name "Wendland C4"
#define kernel_degree 8 /* Degree of the polynomial */
#define kernel_ivals 1  /* Number of branches */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(2.207940))
#define kernel_constant ((float)(495. * M_1_PI / 32.))
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(2.171239))
#define kernel_constant ((float)(9. * M_1_PI))
#elif defined(HYDRO_DIMENSION_1D)
#error "Wendland C4 kernel not defined in 1D."
#endif
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
#define kernel_degree 11 /* Degree of the polynomial */
#define kernel_ivals 1   /* Number of branches */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma ((float)(2.449490))
#define kernel_constant ((float)(1365. * M_1_PI / 64.))
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma ((float)(2.415230))
#define kernel_constant ((float)(78. * M_1_PI / 7.))
#elif defined(HYDRO_DIMENSION_1D)
#error "Wendland C6 kernel not defined in 1D."
#endif
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
#define kernel_gamma_inv ((float)(1. / kernel_gamma))
#define kernel_gamma2 ((float)(kernel_gamma * kernel_gamma))

/* define gamma^d, gamma^(d+1), 1/gamma^d and 1/gamma^(d+1) */
#if defined(HYDRO_DIMENSION_3D)
#define kernel_gamma_dim ((float)(kernel_gamma * kernel_gamma * kernel_gamma))
#define kernel_gamma_dim_plus_one \
  ((float)(kernel_gamma * kernel_gamma * kernel_gamma * kernel_gamma))
#define kernel_gamma_inv_dim \
  ((float)(1. / (kernel_gamma * kernel_gamma * kernel_gamma)))
#define kernel_gamma_inv_dim_plus_one \
  ((float)(1. / (kernel_gamma * kernel_gamma * kernel_gamma * kernel_gamma)))
#elif defined(HYDRO_DIMENSION_2D)
#define kernel_gamma_dim ((float)(kernel_gamma * kernel_gamma))
#define kernel_gamma_dim_plus_one \
  ((float)(kernel_gamma * kernel_gamma * kernel_gamma))
#define kernel_gamma_inv_dim ((float)(1. / (kernel_gamma * kernel_gamma)))
#define kernel_gamma_inv_dim_plus_one \
  ((float)(1. / (kernel_gamma * kernel_gamma * kernel_gamma)))
#elif defined(HYDRO_DIMENSION_1D)
#define kernel_gamma_dim ((float)(kernel_gamma))
#define kernel_gamma_dim_plus_one ((float)(kernel_gamma * kernel_gamma))
#define kernel_gamma_inv_dim ((float)(1. / (kernel_gamma)))
#define kernel_gamma_inv_dim_plus_one \
  ((float)(1. / (kernel_gamma * kernel_gamma)))
#endif

/* The number of branches (floating point conversion) */
#define kernel_ivals_f ((float)(kernel_ivals))

/* Kernel self contribution (i.e. W(0,h)) */
#define kernel_root                                          \
  ((float)(kernel_coeffs[kernel_degree]) * kernel_constant * \
   kernel_gamma_inv_dim)

/* Kernel normalisation constant (volume term) */
#define kernel_norm ((float)(hydro_dimension_unit_sphere * kernel_gamma_dim))

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
__attribute__((always_inline)) INLINE static void kernel_deval(
    float u, float *restrict W, float *restrict dW_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  const float x = u * kernel_gamma_inv;

  /* Pick the correct branch of the kernel */
  const int temp = (int)(x * kernel_ivals_f);
  const int ind = temp > kernel_ivals ? kernel_ivals : temp;
  const float *const coeffs = &kernel_coeffs[ind * (kernel_degree + 1)];

  /* First two terms of the polynomial ... */
  float w = coeffs[0] * x + coeffs[1];
  float dw_dx = coeffs[0];

  /* ... and the rest of them */
  for (int k = 2; k <= kernel_degree; k++) {
    dw_dx = dw_dx * x + w;
    w = x * w + coeffs[k];
  }

  w = max(w, 0.f);
  dw_dx = min(dw_dx, 0.f);

  /* Return everything */
  *W = w * kernel_constant * kernel_gamma_inv_dim;
  *dW_dx = dw_dx * kernel_constant * kernel_gamma_inv_dim_plus_one;
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
__attribute__((always_inline)) INLINE static void kernel_eval(
    float u, float *restrict W) {

  /* Go to the range [0,1[ from [0,H[ */
  const float x = u * kernel_gamma_inv;

  /* Pick the correct branch of the kernel */
  const int temp = (int)(x * kernel_ivals_f);
  const int ind = temp > kernel_ivals ? kernel_ivals : temp;
  const float *const coeffs = &kernel_coeffs[ind * (kernel_degree + 1)];

  /* First two terms of the polynomial ... */
  float w = coeffs[0] * x + coeffs[1];

  /* ... and the rest of them */
  for (int k = 2; k <= kernel_degree; k++) w = x * w + coeffs[k];

  w = max(w, 0.f);

  /* Return everything */
  *W = w * kernel_constant * kernel_gamma_inv_dim;
}

/**
 * @brief Computes the kernel function derivative.
 *
 * The kernel function needs to be mutliplied by \f$h^{-d}\f$ and the gradient
 * by \f$h^{-(d+1)}\f$, where \f$d\f$ is the dimensionality of the problem.
 *
 * Returns 0 if \f$u > \gamma = H/h\f$.
 *
 * @param u The ratio of the distance to the smoothing length \f$u = x/h\f$.
 * @param dW_dx (return) The norm of the gradient of \f$|\nabla W(x,h)|\f$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_dWdx(
    float u, float *restrict dW_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  const float x = u * kernel_gamma_inv;

  /* Pick the correct branch of the kernel */
  const int temp = (int)(x * kernel_ivals_f);
  const int ind = temp > kernel_ivals ? kernel_ivals : temp;
  const float *const coeffs = &kernel_coeffs[ind * (kernel_degree + 1)];

  /* First two terms of the polynomial ... */
  float dw_dx = ((float)kernel_degree * coeffs[0] * x) +
                (float)(kernel_degree - 1) * coeffs[1];

  /* ... and the rest of them */
  for (int k = 2; k < kernel_degree; k++)
    dw_dx = dw_dx * x + (float)(kernel_degree - k) * coeffs[k];

  dw_dx = min(dw_dx, 0.f);

  /* Return everything */
  *dW_dx = dw_dx * kernel_constant * kernel_gamma_inv_dim_plus_one;
}

  /* -------------------------------------------------------------------------
   */

#ifdef WITH_OLD_VECTORIZATION
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
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);

  /* Load x and get the interval id. */
  vector ind;
  ind.m =
      vec_ftoi(vec_fmin(vec_mul(x.v, kernel_ivals_vec.v), kernel_ivals_vec.v));

  /* load the coefficients. */
  vector c[kernel_degree + 1];
  for (int k = 0; k < VEC_SIZE; k++)
    for (int j = 0; j < kernel_degree + 1; j++)
      c[j].f[k] = kernel_coeffs[ind.i[k] * (kernel_degree + 1) + j];

  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(c[0].v, x.v, c[1].v);
  dw_dx->v = c[0].v;

  /* And we're off! */
  for (int k = 2; k <= kernel_degree; k++) {
    dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
    w->v = vec_fma(x.v, w->v, c[k].v);
  }

  /* Return everything */
  w->v =
      vec_mul(w->v, vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
}
#endif

#ifdef WITH_VECTORIZATION

static const vector kernel_gamma_inv_vec = FILL_VEC((float)kernel_gamma_inv);

static const vector kernel_ivals_vec = FILL_VEC((float)kernel_ivals);

static const vector kernel_constant_vec = FILL_VEC((float)kernel_constant);

static const vector kernel_gamma_inv_dim_vec =
    FILL_VEC((float)kernel_gamma_inv_dim);

static const vector kernel_gamma_inv_dim_plus_one_vec =
    FILL_VEC((float)kernel_gamma_inv_dim_plus_one);

/* Define constant vectors for the Wendland C2 and Cubic Spline kernel
 * coefficients. */
#ifdef WENDLAND_C2_KERNEL
static const vector wendland_const_c0 = FILL_VEC(4.f);
static const vector wendland_const_c1 = FILL_VEC(-15.f);
static const vector wendland_const_c2 = FILL_VEC(20.f);
static const vector wendland_const_c3 = FILL_VEC(-10.f);
static const vector wendland_const_c4 = FILL_VEC(0.f);
static const vector wendland_const_c5 = FILL_VEC(1.f);

static const vector wendland_dwdx_const_c0 = FILL_VEC(20.f);
static const vector wendland_dwdx_const_c1 = FILL_VEC(-60.f);
static const vector wendland_dwdx_const_c2 = FILL_VEC(60.f);
static const vector wendland_dwdx_const_c3 = FILL_VEC(-20.f);
#elif defined(CUBIC_SPLINE_KERNEL)
/* First region 0 < u < 0.5 */
static const vector cubic_1_const_c0 = FILL_VEC(3.f);
static const vector cubic_1_const_c1 = FILL_VEC(-3.f);
static const vector cubic_1_const_c2 = FILL_VEC(0.f);
static const vector cubic_1_const_c3 = FILL_VEC(0.5f);
static const vector cubic_1_dwdx_const_c0 = FILL_VEC(9.f);
static const vector cubic_1_dwdx_const_c1 = FILL_VEC(-6.f);
static const vector cubic_1_dwdx_const_c2 = FILL_VEC(0.f);

/* Second region 0.5 <= u < 1 */
static const vector cubic_2_const_c0 = FILL_VEC(-1.f);
static const vector cubic_2_const_c1 = FILL_VEC(3.f);
static const vector cubic_2_const_c2 = FILL_VEC(-3.f);
static const vector cubic_2_const_c3 = FILL_VEC(1.f);
static const vector cubic_2_dwdx_const_c0 = FILL_VEC(-3.f);
static const vector cubic_2_dwdx_const_c1 = FILL_VEC(6.f);
static const vector cubic_2_dwdx_const_c2 = FILL_VEC(-3.f);
static const vector cond = FILL_VEC(0.5f);
#endif

/**
 * @brief Computes the kernel function and its derivative for two particles
 * using vectors. The return value is undefined if $u > \\gamma = H/h$.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param w (return) The value of the kernel function $W(x,h)$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 */
__attribute__((always_inline)) INLINE static void kernel_deval_1_vec(
    vector *u, vector *w, vector *dw_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(wendland_const_c0.v, x.v, wendland_const_c1.v);
  dw_dx->v = wendland_const_c0.v;

  /* Calculate the polynomial interleaving vector operations */
  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c3.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  w->v = vec_mul(x.v, w->v); /* wendland_const_c4 is zero. */

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c5.v);
#elif defined(CUBIC_SPLINE_KERNEL)
  vector w2, dw_dx2;
  mask_t mask_reg;

  /* Form a mask for one part of the kernel. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  vec_create_mask(mask_reg, vec_cmp_gte(x.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using a mask. */

  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(cubic_1_const_c0.v, x.v, cubic_1_const_c1.v);
  w2.v = vec_fma(cubic_2_const_c0.v, x.v, cubic_2_const_c1.v);
  dw_dx->v = cubic_1_const_c0.v;
  dw_dx2.v = cubic_2_const_c0.v;

  /* Calculate the polynomial interleaving vector operations. */
  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2.v = vec_fma(dw_dx2.v, x.v, w2.v);
  w->v = vec_mul(x.v, w->v); /* cubic_1_const_c2 is zero. */
  w2.v = vec_fma(x.v, w2.v, cubic_2_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2.v = vec_fma(dw_dx2.v, x.v, w2.v);
  w->v = vec_fma(x.v, w->v, cubic_1_const_c3.v);
  w2.v = vec_fma(x.v, w2.v, cubic_2_const_c3.v);

  /* Blend both kernel regions into one vector (mask out unneeded values). */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  w->v = vec_blend(mask_reg, w->v, w2.v);
  dw_dx->v = vec_blend(mask_reg, dw_dx->v, dw_dx2.v);

#else
#error "Vectorisation not supported for this kernel!!!"
#endif

  /* Return everyting */
  w->v =
      vec_mul(w->v, vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
}

/**
 * @brief Computes the kernel function and its derivative for two particles
 * using interleaved vectors. The return value is undefined if $u > \\gamma =
 * H/h$.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param w (return) The value of the kernel function $W(x,h)$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 * @param u2 The ratio of the distance to the smoothing length $u = x/h$ for
 * second particle.
 * @param w2 (return) The value of the kernel function $W(x,h)$ for second
 * particle.
 * @param dw_dx2 (return) The norm of the gradient of $|\\nabla W(x,h)|$ for
 * second particle.
 */
__attribute__((always_inline)) INLINE static void kernel_deval_2_vec(
    vector *u, vector *w, vector *dw_dx, vector *u2, vector *w2,
    vector *dw_dx2) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x, x2;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);
  x2.v = vec_mul(u2->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(wendland_const_c0.v, x.v, wendland_const_c1.v);
  w2->v = vec_fma(wendland_const_c0.v, x2.v, wendland_const_c1.v);
  dw_dx->v = wendland_const_c0.v;
  dw_dx2->v = wendland_const_c0.v;

  /* Calculate the polynomial interleaving vector operations */
  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c2.v);
  w2->v = vec_fma(x2.v, w2->v, wendland_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c3.v);
  w2->v = vec_fma(x2.v, w2->v, wendland_const_c3.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  w->v = vec_mul(x.v, w->v);    /* wendland_const_c4 is zero. */
  w2->v = vec_mul(x2.v, w2->v); /* wendland_const_c4 is zero. */

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  w->v = vec_fma(x.v, w->v, wendland_const_c5.v);
  w2->v = vec_fma(x2.v, w2->v, wendland_const_c5.v);

  /* Return everything */
  w->v =
      vec_mul(w->v, vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  w2->v = vec_mul(w2->v,
                  vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
  dw_dx2->v = vec_mul(dw_dx2->v, vec_mul(kernel_constant_vec.v,
                                         kernel_gamma_inv_dim_plus_one_vec.v));
#elif defined(CUBIC_SPLINE_KERNEL)
  vector w_2, dw_dx_2;
  vector w2_2, dw_dx2_2;
  mask_t mask_reg, mask_reg_v2;

  /* Form a mask for one part of the kernel for each vector. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  vec_create_mask(mask_reg, vec_cmp_gte(x.v, cond.v));     /* 0.5 < x < 1 */
  vec_create_mask(mask_reg_v2, vec_cmp_gte(x2.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using masks. */

  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(cubic_1_const_c0.v, x.v, cubic_1_const_c1.v);
  w2->v = vec_fma(cubic_1_const_c0.v, x2.v, cubic_1_const_c1.v);
  w_2.v = vec_fma(cubic_2_const_c0.v, x.v, cubic_2_const_c1.v);
  w2_2.v = vec_fma(cubic_2_const_c0.v, x2.v, cubic_2_const_c1.v);
  dw_dx->v = cubic_1_const_c0.v;
  dw_dx2->v = cubic_1_const_c0.v;
  dw_dx_2.v = cubic_2_const_c0.v;
  dw_dx2_2.v = cubic_2_const_c0.v;

  /* Calculate the polynomial interleaving vector operations. */
  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  dw_dx_2.v = vec_fma(dw_dx_2.v, x.v, w_2.v);
  dw_dx2_2.v = vec_fma(dw_dx2_2.v, x2.v, w2_2.v);
  w->v = vec_mul(x.v, w->v);    /* cubic_1_const_c2 is zero. */
  w2->v = vec_mul(x2.v, w2->v); /* cubic_1_const_c2 is zero. */
  w_2.v = vec_fma(x.v, w_2.v, cubic_2_const_c2.v);
  w2_2.v = vec_fma(x2.v, w2_2.v, cubic_2_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, w->v);
  dw_dx2->v = vec_fma(dw_dx2->v, x2.v, w2->v);
  dw_dx_2.v = vec_fma(dw_dx_2.v, x.v, w_2.v);
  dw_dx2_2.v = vec_fma(dw_dx2_2.v, x2.v, w2_2.v);
  w->v = vec_fma(x.v, w->v, cubic_1_const_c3.v);
  w2->v = vec_fma(x2.v, w2->v, cubic_1_const_c3.v);
  w_2.v = vec_fma(x.v, w_2.v, cubic_2_const_c3.v);
  w2_2.v = vec_fma(x2.v, w2_2.v, cubic_2_const_c3.v);

  /* Blend both kernel regions into one vector (mask out unneeded values). */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  w->v = vec_blend(mask_reg, w->v, w_2.v);
  w2->v = vec_blend(mask_reg_v2, w2->v, w2_2.v);
  dw_dx->v = vec_blend(mask_reg, dw_dx->v, dw_dx_2.v);
  dw_dx2->v = vec_blend(mask_reg_v2, dw_dx2->v, dw_dx2_2.v);

  /* Return everything */
  w->v =
      vec_mul(w->v, vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  w2->v = vec_mul(w2->v,
                  vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
  dw_dx2->v = vec_mul(dw_dx2->v, vec_mul(kernel_constant_vec.v,
                                         kernel_gamma_inv_dim_plus_one_vec.v));

#endif
}

/**
 * @brief Computes the kernel function for two particles
 * using vectors. The return value is undefined if $u > \\gamma = H/h$.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param w (return) The value of the kernel function $W(x,h)$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_W_vec(vector *u,
                                                                    vector *w) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(wendland_const_c0.v, x.v, wendland_const_c1.v);

  /* Calculate the polynomial interleaving vector operations */
  w->v = vec_fma(x.v, w->v, wendland_const_c2.v);
  w->v = vec_fma(x.v, w->v, wendland_const_c3.v);
  w->v = vec_mul(x.v, w->v); /* wendland_const_c4 is zero.*/
  w->v = vec_fma(x.v, w->v, wendland_const_c5.v);
#elif defined(CUBIC_SPLINE_KERNEL)
  vector w2;
  mask_t mask_reg;

  /* Form a mask for each part of the kernel. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  vec_create_mask(mask_reg, vec_cmp_gte(x.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using masks. */

  /* Init the iteration for Horner's scheme. */
  w->v = vec_fma(cubic_1_const_c0.v, x.v, cubic_1_const_c1.v);
  w2.v = vec_fma(cubic_2_const_c0.v, x.v, cubic_2_const_c1.v);

  /* Calculate the polynomial interleaving vector operations. */
  w->v = vec_mul(x.v, w->v); /* cubic_1_const_c2 is zero */
  w2.v = vec_fma(x.v, w2.v, cubic_2_const_c2.v);

  w->v = vec_fma(x.v, w->v, cubic_1_const_c3.v);
  w2.v = vec_fma(x.v, w2.v, cubic_2_const_c3.v);

  /* Mask out unneeded values. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  w->v = vec_blend(mask_reg, w->v, w2.v);

#else
#error "Vectorisation not supported for this kernel!!!"
#endif

  /* Return everything */
  w->v =
      vec_mul(w->v, vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_vec.v));
}

/**
 * @brief Computes the kernel function derivative for two particles
 * using vectors. The return value is undefined if $u > \\gamma = H/h$.
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_dWdx_vec(
    vector *u, vector *dw_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(wendland_dwdx_const_c0.v, x.v, wendland_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations */
  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c3.v);

  dw_dx->v = vec_mul(dw_dx->v, x.v);

#elif defined(CUBIC_SPLINE_KERNEL)
  vector dw_dx2;
  mask_t mask_reg1, mask_reg2;

  /* Form a mask for each part of the kernel. */
  vec_create_mask(mask_reg1, vec_cmp_lt(x.v, cond.v));  /* 0 < x < 0.5 */
  vec_create_mask(mask_reg2, vec_cmp_gte(x.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using masks. */

  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(cubic_1_dwdx_const_c0.v, x.v, cubic_1_dwdx_const_c1.v);
  dw_dx2.v = vec_fma(cubic_2_dwdx_const_c0.v, x.v, cubic_2_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations. */
  dw_dx->v = vec_mul(dw_dx->v, x.v); /* cubic_1_dwdx_const_c2 is zero. */
  dw_dx2.v = vec_fma(dw_dx2.v, x.v, cubic_2_dwdx_const_c2.v);

  /* Mask out unneeded values. */
  dw_dx->v = vec_and_mask(dw_dx->v, mask_reg1);
  dw_dx2.v = vec_and_mask(dw_dx2.v, mask_reg2);

  /* Added both dwdx and dwdx2 together to form complete result. */
  dw_dx->v = vec_add(dw_dx->v, dw_dx2.v);
#else
#error "Vectorisation not supported for this kernel!!!"
#endif

  /* Return everything */
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
}

/**
 * @brief Computes the kernel function derivative for two particles
 * using vectors.
 *
 * Return 0 if $u > \\gamma = H/h$
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_dWdx_force_vec(
    vector *u, vector *dw_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(wendland_dwdx_const_c0.v, x.v, wendland_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations */
  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c3.v);

  dw_dx->v = vec_mul(dw_dx->v, x.v);

#elif defined(CUBIC_SPLINE_KERNEL)
  vector dw_dx2;
  mask_t mask_reg;

  /* Form a mask for each part of the kernel. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  vec_create_mask(mask_reg, vec_cmp_gte(x.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using masks. */

  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(cubic_1_dwdx_const_c0.v, x.v, cubic_1_dwdx_const_c1.v);
  dw_dx2.v = vec_fma(cubic_2_dwdx_const_c0.v, x.v, cubic_2_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations. */
  dw_dx->v = vec_mul(dw_dx->v, x.v); /* cubic_1_dwdx_const_c2 is zero. */
  dw_dx2.v = vec_fma(dw_dx2.v, x.v, cubic_2_dwdx_const_c2.v);

  /* Mask out unneeded values. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  dw_dx->v = vec_blend(mask_reg, dw_dx->v, dw_dx2.v);

#else
#error "Vectorisation not supported for this kernel!!!"
#endif

  /* Mask out result for particles that lie outside of the kernel function. */
  mask_t mask;
  vec_create_mask(mask, vec_cmp_lt(x.v, vec_set1(1.f))); /* x < 1 */

  dw_dx->v = vec_and_mask(dw_dx->v, mask);

  /* Return everything */
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
}

/**
 * @brief Computes the kernel function derivative for two particles
 * using interleaved vectors.
 *
 * Return 0 if $u > \\gamma = H/h$
 *
 * @param u The ratio of the distance to the smoothing length $u = x/h$.
 * @param dw_dx (return) The norm of the gradient of $|\\nabla W(x,h)|$.
 * @param u_2 The ratio of the distance to the smoothing length $u = x/h$ for
 * second particle.
 * @param dw_dx_2 (return) The norm of the gradient of $|\\nabla W(x,h)|$ for
 * second particle.
 */
__attribute__((always_inline)) INLINE static void kernel_eval_dWdx_force_2_vec(
    vector *u, vector *dw_dx, vector *u_2, vector *dw_dx_2) {

  /* Go to the range [0,1[ from [0,H[ */
  vector x, x_2;
  x.v = vec_mul(u->v, kernel_gamma_inv_vec.v);
  x_2.v = vec_mul(u_2->v, kernel_gamma_inv_vec.v);

#ifdef WENDLAND_C2_KERNEL
  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(wendland_dwdx_const_c0.v, x.v, wendland_dwdx_const_c1.v);
  dw_dx_2->v =
      vec_fma(wendland_dwdx_const_c0.v, x_2.v, wendland_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations */
  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c2.v);
  dw_dx_2->v = vec_fma(dw_dx_2->v, x_2.v, wendland_dwdx_const_c2.v);

  dw_dx->v = vec_fma(dw_dx->v, x.v, wendland_dwdx_const_c3.v);
  dw_dx_2->v = vec_fma(dw_dx_2->v, x_2.v, wendland_dwdx_const_c3.v);

  dw_dx->v = vec_mul(dw_dx->v, x.v);
  dw_dx_2->v = vec_mul(dw_dx_2->v, x_2.v);

#elif defined(CUBIC_SPLINE_KERNEL)
  vector dw_dx2, dw_dx2_2;
  mask_t mask_reg;
  mask_t mask_reg_v2;

  /* Form a mask for one part of the kernel. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  vec_create_mask(mask_reg, vec_cmp_gte(x.v, cond.v));      /* 0.5 < x < 1 */
  vec_create_mask(mask_reg_v2, vec_cmp_gte(x_2.v, cond.v)); /* 0.5 < x < 1 */

  /* Work out w for both regions of the kernel and combine the results together
   * using masks. */

  /* Init the iteration for Horner's scheme. */
  dw_dx->v = vec_fma(cubic_1_dwdx_const_c0.v, x.v, cubic_1_dwdx_const_c1.v);
  dw_dx_2->v = vec_fma(cubic_1_dwdx_const_c0.v, x_2.v, cubic_1_dwdx_const_c1.v);
  dw_dx2.v = vec_fma(cubic_2_dwdx_const_c0.v, x.v, cubic_2_dwdx_const_c1.v);
  dw_dx2_2.v = vec_fma(cubic_2_dwdx_const_c0.v, x_2.v, cubic_2_dwdx_const_c1.v);

  /* Calculate the polynomial interleaving vector operations. */
  dw_dx->v = vec_mul(dw_dx->v, x.v);       /* cubic_1_dwdx_const_c2 is zero. */
  dw_dx_2->v = vec_mul(dw_dx_2->v, x_2.v); /* cubic_1_dwdx_const_c2 is zero. */
  dw_dx2.v = vec_fma(dw_dx2.v, x.v, cubic_2_dwdx_const_c2.v);
  dw_dx2_2.v = vec_fma(dw_dx2_2.v, x_2.v, cubic_2_dwdx_const_c2.v);

  /* Mask out unneeded values. */
  /* Only need the mask for one region as the vec_blend defaults to the vector
   * when the mask is 0.*/
  dw_dx->v = vec_blend(mask_reg, dw_dx->v, dw_dx2.v);
  dw_dx_2->v = vec_blend(mask_reg_v2, dw_dx_2->v, dw_dx2_2.v);

#else
#error "Vectorisation not supported for this kernel!!!"
#endif

  /* Mask out result for particles that lie outside of the kernel function. */
  mask_t mask, mask_2;
  vec_create_mask(mask, vec_cmp_lt(x.v, vec_set1(1.f)));     /* x < 1 */
  vec_create_mask(mask_2, vec_cmp_lt(x_2.v, vec_set1(1.f))); /* x < 1 */

  dw_dx->v = vec_and_mask(dw_dx->v, mask);
  dw_dx_2->v = vec_and_mask(dw_dx_2->v, mask_2);

  /* Return everything */
  dw_dx->v = vec_mul(dw_dx->v, vec_mul(kernel_constant_vec.v,
                                       kernel_gamma_inv_dim_plus_one_vec.v));
  dw_dx_2->v = vec_mul(
      dw_dx_2->v,
      vec_mul(kernel_constant_vec.v, kernel_gamma_inv_dim_plus_one_vec.v));
}

#endif /* WITH_VECTORIZATION */

/* Some cross-check functions */
void hydro_kernel_dump(int N);

#endif  // SWIFT_KERNEL_HYDRO_H
