/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Katy Proctor (katy.proctor@fysik.su.se)
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
#ifndef SWIFT_KERNEL_SIDM_H
#define SWIFT_KERNEL_SIDM_H

/**
 * @file kernel_sidm.h
 * @brief Kernel functions for SIDM.
 *
 * Two kernels are defined here, both chosen at configure time independently
 * of the hydro kernel:
 *  - The density kernel (--with-sidm-kernel), used for the SIDM density loop
 *    and hence the smoothing-length iteration and the neighbour search radius
 *    H = sidm_kernel_gamma * h.
 *  - The overlap kernel (--with-sidm-overlap-kernel), used to compute the
 *    pair overlap integral in the SIDM force (scattering) loop.
 *
 * All names carry a sidm_ prefix so that they do not clash with the hydro
 * kernel defined in kernel_hydro.h.
 *
 * All density-kernel constants and coefficients are taken from table 1 of
 * Dehnen & Aly, MNRAS, 425, pp. 1068-1082 (2012).
 *
 * The overlap kernel gives the pair scattering rate per unit (sigma/m) m v,ß
 *
 *   Lambda_ij = int W(|x - x_i|, h_i) W(|x - x_j|, h_j) d^3x,
 *
 * Because each W integrates to 1, so does
 * Lambda_ij over all separations r = |x_i - x_j|. The force loop only sees
 * pairs with r < sidm_kernel_gamma * h_max, so both overlap kernels below are
 * truncated there and divided by
 *
 *   N(s) = int_0^{sidm_kernel_gamma} 4 pi r^2 Lambda(r, s) dr   (h_max = 1),
 *
 * with s = h_min / h_max, which makes the truncated Lambda integrate to 1.
 */

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "dimension.h"
#include "error.h"
#include "inline.h"
#include "minmax.h"

/* ========================================================================= */
/* Density / neighbour-search kernel                                         */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
#if defined(SIDM_CUBIC_SPLINE_KERNEL)

/* Coefficients for the kernel. */
#define sidm_kernel_name "Cubic spline (M4)"
#define sidm_kernel_degree 3 /*!< Degree of the polynomial */
#define sidm_kernel_ivals 2  /*!< Number of branches */
#if defined(HYDRO_DIMENSION_3D)
#define sidm_kernel_gamma ((float)(1.825742))
#define sidm_kernel_constant ((float)(16. * M_1_PI))
#elif defined(HYDRO_DIMENSION_2D)
#define sidm_kernel_gamma ((float)(1.778002))
#define sidm_kernel_constant ((float)(80. * M_1_PI / 7.))
#elif defined(HYDRO_DIMENSION_1D)
#define sidm_kernel_gamma ((float)(1.732051))
#define sidm_kernel_constant ((float)(8. / 3.))
#endif
static const float
    sidm_kernel_coeffs[(sidm_kernel_degree + 1) * (sidm_kernel_ivals + 1)]
    __attribute__((aligned(16))) = {3.f,  -3.f, 0.f,  0.5f, /* 0 < u < 0.5 */
                                    -1.f, 3.f,  -3.f, 1.f,  /* 0.5 < u < 1 */
                                    0.f,  0.f,  0.f,  0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(SIDM_QUARTIC_SPLINE_KERNEL)

/* Coefficients for the kernel. */
#define sidm_kernel_name "Quartic spline (M5)"
#define sidm_kernel_degree 4 /* Degree of the polynomial */
#define sidm_kernel_ivals 5  /* Number of branches */
#if defined(HYDRO_DIMENSION_3D)
#define sidm_kernel_gamma ((float)(2.018932))
#define sidm_kernel_constant ((float)(15625. * M_1_PI / 512.))
#elif defined(HYDRO_DIMENSION_2D)
#define sidm_kernel_gamma ((float)(1.977173))
#define sidm_kernel_constant ((float)(46875. * M_1_PI / 2398.))
#elif defined(HYDRO_DIMENSION_1D)
#define sidm_kernel_gamma ((float)(1.936492))
#define sidm_kernel_constant ((float)(3125. / 768.))
#endif
static const float
    sidm_kernel_coeffs[(sidm_kernel_degree + 1) * (sidm_kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        6.f,  0.f,  -2.4f, 0.f,   0.368f, /* 0 < u < 0.2 */
        -4.f, 8.f,  -4.8f, 0.32f, 0.352f, /* 0.2 < u < 0.4 */
        -4.f, 8.f,  -4.8f, 0.32f, 0.352f, /* 0.4 < u < 0.6 */
        1.f,  -4.f, 6.f,   -4.f,  1.f,    /* 0.6 < u < 0.8 */
        1.f,  -4.f, 6.f,   -4.f,  1.f,    /* 0.8 < u < 1 */
        0.f,  0.f,  0.f,   0.f,   0.f};   /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(SIDM_QUINTIC_SPLINE_KERNEL)

/* Coefficients for the kernel. */
#define sidm_kernel_name "Quintic spline (M6)"
#define sidm_kernel_degree 5 /* Degree of the polynomial */
#define sidm_kernel_ivals 3  /* Number of branches */
#if defined(HYDRO_DIMENSION_3D)
#define sidm_kernel_gamma ((float)(2.195775))
#define sidm_kernel_constant ((float)(2187. * M_1_PI / 40.))
#elif defined(HYDRO_DIMENSION_2D)
#define sidm_kernel_gamma ((float)(2.158131))
#define sidm_kernel_constant ((float)(15309. * M_1_PI / 478.))
#elif defined(HYDRO_DIMENSION_1D)
#define sidm_kernel_gamma ((float)(2.121321))
#define sidm_kernel_constant ((float)(243. / 40.))
#endif
static const float
    sidm_kernel_coeffs[(sidm_kernel_degree + 1) * (sidm_kernel_ivals + 1)]
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
#elif defined(SIDM_WENDLAND_C2_KERNEL)

/* Coefficients for the kernel. */
#define sidm_kernel_name "Wendland C2"
#define sidm_kernel_degree 5 /* Degree of the polynomial */
#define sidm_kernel_ivals 1  /* Number of branches */
#if defined(HYDRO_DIMENSION_1D)
/* Wendland C* have different form in 1D than 2D/3D */
#define sidm_kernel_gamma ((float)(1.620185))
#define sidm_kernel_constant ((float)(5. / 4.))
static const float
    sidm_kernel_coeffs[(sidm_kernel_degree + 1) * (sidm_kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        0.f, -3.f, 8.f, -6.f, 0.f, 1.f, /* 0 < u < 1 */
        0.f, 0.f,  0.f, 0.f,  0.f, 0.f};
#else
#if defined(HYDRO_DIMENSION_3D)
#define sidm_kernel_gamma ((float)(1.936492))
#define sidm_kernel_constant ((float)(21. * M_1_PI / 2.))
#elif defined(HYDRO_DIMENSION_2D)
#define sidm_kernel_gamma ((float)(1.897367))
#define sidm_kernel_constant ((float)(7. * M_1_PI))
#endif
static const float
    sidm_kernel_coeffs[(sidm_kernel_degree + 1) * (sidm_kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        4.f, -15.f, 20.f, -10.f, 0.f, 1.f,  /* 0 < u < 1 */
        0.f, 0.f,   0.f,  0.f,   0.f, 0.f}; /* 1 < u */
#endif
/* ------------------------------------------------------------------------- */
#elif defined(SIDM_WENDLAND_C4_KERNEL)

/* Coefficients for the kernel. */
#define sidm_kernel_name "Wendland C4"
#define sidm_kernel_degree 8 /* Degree of the polynomial */
#define sidm_kernel_ivals 1  /* Number of branches */
#if defined(HYDRO_DIMENSION_3D)
#define sidm_kernel_gamma ((float)(2.207940))
#define sidm_kernel_constant ((float)(495. * M_1_PI / 32.))
#elif defined(HYDRO_DIMENSION_2D)
#define sidm_kernel_gamma ((float)(2.171239))
#define sidm_kernel_constant ((float)(9. * M_1_PI))
#elif defined(HYDRO_DIMENSION_1D)
#error "Wendland C4 SIDM kernel not defined in 1D."
#endif
static const float
    sidm_kernel_coeffs[(sidm_kernel_degree + 1) * (sidm_kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        11.666667f, -64.f,       140.f, -149.333333f, 70.f,
        0.f,        -9.3333333f, 0.f,   1.f, /* 0 < u < 1 */
        0.f,        0.f,         0.f,   0.f,          0.f,
        0.f,        0.f,         0.f,   0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#elif defined(SIDM_WENDLAND_C6_KERNEL)

/* Coefficients for the kernel. */
#define sidm_kernel_name "Wendland C6"
#define sidm_kernel_degree 11 /* Degree of the polynomial */
#define sidm_kernel_ivals 1   /* Number of branches */
#if defined(HYDRO_DIMENSION_3D)
#define sidm_kernel_gamma ((float)(2.449490))
#define sidm_kernel_constant ((float)(1365. * M_1_PI / 64.))
#elif defined(HYDRO_DIMENSION_2D)
#define sidm_kernel_gamma ((float)(2.415230))
#define sidm_kernel_constant ((float)(78. * M_1_PI / 7.))
#elif defined(HYDRO_DIMENSION_1D)
#error "Wendland C6 SIDM kernel not defined in 1D."
#endif
static const float
    sidm_kernel_coeffs[(sidm_kernel_degree + 1) * (sidm_kernel_ivals + 1)]
    __attribute__((aligned(16))) = {
        32.f, -231.f, 704.f, -1155.f, 1056.f, -462.f,
        0.f,  66.f,   0.f,   -11.f,   0.f,    1.f, /* 0 < u < 1 */
        0.f,  0.f,    0.f,   0.f,     0.f,    0.f,
        0.f,  0.f,    0.f,   0.f,     0.f,    0.f}; /* 1 < u */

/* ------------------------------------------------------------------------- */
#else

#error "An SIDM kernel function must be chosen at configure time !!"

/* ------------------------------------------------------------------------- */
#endif

/* Ok, now comes the real deal. */

/* First some powers of gamma = H/h */
#define sidm_kernel_gamma_inv ((float)(1. / sidm_kernel_gamma))
#define sidm_kernel_gamma2 ((float)(sidm_kernel_gamma * sidm_kernel_gamma))

/* define gamma^d, gamma^(d+1), 1/gamma^d and 1/gamma^(d+1) */
#if defined(HYDRO_DIMENSION_3D)
#define sidm_kernel_gamma_dim \
  ((float)(sidm_kernel_gamma * sidm_kernel_gamma * sidm_kernel_gamma))
#define sidm_kernel_gamma_dim_plus_one                                 \
  ((float)(sidm_kernel_gamma * sidm_kernel_gamma * sidm_kernel_gamma * \
           sidm_kernel_gamma))
#define sidm_kernel_gamma_inv_dim \
  ((float)(1. / (sidm_kernel_gamma * sidm_kernel_gamma * sidm_kernel_gamma)))
#define sidm_kernel_gamma_inv_dim_plus_one                                   \
  ((float)(1. / (sidm_kernel_gamma * sidm_kernel_gamma * sidm_kernel_gamma * \
                 sidm_kernel_gamma)))
#elif defined(HYDRO_DIMENSION_2D)
#define sidm_kernel_gamma_dim ((float)(sidm_kernel_gamma * sidm_kernel_gamma))
#define sidm_kernel_gamma_dim_plus_one \
  ((float)(sidm_kernel_gamma * sidm_kernel_gamma * sidm_kernel_gamma))
#define sidm_kernel_gamma_inv_dim \
  ((float)(1. / (sidm_kernel_gamma * sidm_kernel_gamma)))
#define sidm_kernel_gamma_inv_dim_plus_one \
  ((float)(1. / (sidm_kernel_gamma * sidm_kernel_gamma * sidm_kernel_gamma)))
#elif defined(HYDRO_DIMENSION_1D)
#define sidm_kernel_gamma_dim ((float)(sidm_kernel_gamma))
#define sidm_kernel_gamma_dim_plus_one \
  ((float)(sidm_kernel_gamma * sidm_kernel_gamma))
#define sidm_kernel_gamma_inv_dim ((float)(1. / (sidm_kernel_gamma)))
#define sidm_kernel_gamma_inv_dim_plus_one \
  ((float)(1. / (sidm_kernel_gamma * sidm_kernel_gamma)))
#endif

/* The number of branches (floating point conversion) */
#define sidm_kernel_ivals_f ((float)(sidm_kernel_ivals))

/* Kernel self contribution (i.e. W(0,h)) */
#define sidm_kernel_root                                                    \
  ((float)(sidm_kernel_coeffs[sidm_kernel_degree]) * sidm_kernel_constant * \
   sidm_kernel_gamma_inv_dim)

/* Kernel normalisation constant (volume term) */
#define sidm_kernel_norm \
  ((float)(hydro_dimension_unit_sphere * sidm_kernel_gamma_dim))

/* ------------------------------------------------------------------------- */

/**
 * @brief Computes the SIDM kernel function and its derivative.
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
__attribute__((always_inline)) INLINE static void sidm_kernel_deval(
    float u, float *restrict W, float *restrict dW_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  const float x = u * sidm_kernel_gamma_inv;

  /* Pick the correct branch of the kernel */
  const int temp = (int)(x * sidm_kernel_ivals_f);
  const int ind = temp > sidm_kernel_ivals ? sidm_kernel_ivals : temp;
  const float *const coeffs =
      &sidm_kernel_coeffs[ind * (sidm_kernel_degree + 1)];

  /* First two terms of the polynomial ... */
  float w = coeffs[0] * x + coeffs[1];
  float dw_dx = coeffs[0];

  /* ... and the rest of them */
  for (int k = 2; k <= sidm_kernel_degree; k++) {
    dw_dx = dw_dx * x + w;
    w = x * w + coeffs[k];
  }

  w = max(w, 0.f);
  dw_dx = min(dw_dx, 0.f);

  /* Return everything */
  *W = w * sidm_kernel_constant * sidm_kernel_gamma_inv_dim;
  *dW_dx = dw_dx * sidm_kernel_constant * sidm_kernel_gamma_inv_dim_plus_one;
}

/**
 * @brief Computes the SIDM kernel function.
 *
 * The kernel function needs to be mutliplied by \f$h^{-d}\f$,
 * where \f$d\f$ is the dimensionality of the problem.
 *
 * Returns 0 if \f$u > \gamma = H/h\f$
 *
 * @param u The ratio of the distance to the smoothing length \f$u = x/h\f$.
 * @param W (return) The value of the kernel function \f$W(x,h)\f$.
 */
__attribute__((always_inline)) INLINE static void sidm_kernel_eval(
    float u, float *restrict W) {

  /* Go to the range [0,1[ from [0,H[ */
  const float x = u * sidm_kernel_gamma_inv;

  /* Pick the correct branch of the kernel */
  const int temp = (int)(x * sidm_kernel_ivals_f);
  const int ind = temp > sidm_kernel_ivals ? sidm_kernel_ivals : temp;
  const float *const coeffs =
      &sidm_kernel_coeffs[ind * (sidm_kernel_degree + 1)];

  /* First two terms of the polynomial ... */
  float w = coeffs[0] * x + coeffs[1];

  /* ... and the rest of them */
  for (int k = 2; k <= sidm_kernel_degree; k++) w = x * w + coeffs[k];

  w = max(w, 0.f);

  /* Return everything */
  *W = w * sidm_kernel_constant * sidm_kernel_gamma_inv_dim;
}

/* ========================================================================= */
/* Overlap (interaction) kernel                                              */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
#if defined(SIDM_OVERLAP_TOPHAT)

#if !defined(HYDRO_DIMENSION_3D)
#error "The top-hat SIDM overlap kernel is only defined in 3D."
#endif

#define sidm_overlap_kernel_name "Top-hat"

/**
 * @brief Normalisation N(s) of the truncated top-hat overlap (h_max = 1).
 *
 * The untruncated overlap integrates to 1, so N = 1 - M, where M is the part
 * at r > sidm_kernel_gamma. M is non-zero only when the lens reaches past the
 * search radius, i.e. 1 + s > sidm_kernel_gamma.
 *
 * With Lambda = V_lens / (V_i V_j), V_i = 4 pi / 3, V_j = 4 pi s^3 / 3 and
 * the substitution t = 1 + s - r, the integrand 4 pi r^2 Lambda becomes
 * 3 / (16 s^3) t^2 (1 + s - t) (t^2 - 4 (1 + s) t + 12 s), so that
 *
 *   M = 3 / (16 s^3) int_0^T t^2 (1 + s - t) (t^2 - 4 (1 + s) t + 12 s) dt
 *     = 3 / (16 s^3) T^3 [4 a s - (a^2 + 3 s) T + a T^2 - T^3 / 6],
 *
 * with a = 1 + s and T = a - sidm_kernel_gamma.
 *
 * @param s The ratio h_min / h_max (in [0, 1]).
 */
INLINE static float sidm_kernel_overlap_tophat_norm(const float s) {

  const float a = 1.f + s;
  const float T = a - sidm_kernel_gamma;
  if (T <= 0.f) return 1.f;

  const float T3 = T * T * T;
  const float integral =
      T3 * (4.f * a * s - (a * a + 3.f * s) * T + a * T * T - T3 / 6.f);

  return 1.f - 3.f * integral / (16.f * s * s * s);
}

/**
 * @brief Kernel overlap of two top-hat kernels of radius h_i and h_j.
 *
 * For W = 1 / V inside a sphere of radius h (and 0 outside), the overlap
 * integral is Lambda_ij = V_lens / (V_i V_j), with V_lens the volume of the
 * intersection of the two spheres:
 *
 *   V_lens = pi (h_i + h_j - r)^2 [r^2 + 2 r (h_i + h_j) - 3 (h_i - h_j)^2]
 *            / (12 r)                                 for |h_i - h_j| < r,
 *   V_lens = 4 pi h_min^3 / 3                         for r <= |h_i - h_j|,
 *
 * and 0 for r >= h_i + h_j. The result is truncated at the search radius and divided by
 * N(s) (see sidm_kernel_overlap_tophat_norm()).
 *
 *
 * @param r Comoving distance between the two particles.ß
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 */
INLINE static float sidm_kernel_overlap(const float r, const float hi,
                                        const float hj) {

  const float h_max = max(hi, hj);
  const float h_min = min(hi, hj);

  /* Outside the search radius or the two spheres do not intersect */
  if (r >= sidm_kernel_gamma * h_max) return 0.f;
  if (r >= hi + hj) return 0.f;

  const float sphere_volume = (4.f / 3.f) * M_PI;
  const float V_i = sphere_volume * hi * hi * hi;
  const float V_j = sphere_volume * hj * hj * hj;

  float V_lens;
  if (r <= h_max - h_min) {
    /* The smaller sphere lies entirely inside the larger one */
    V_lens = sphere_volume * h_min * h_min * h_min;
  } else {
    /* Lens-shaped intersection */
    const float depth = hi + hj - r;
    const float dh = hi - hj;
    V_lens = M_PI * depth * depth * (r * r + 2.f * r * (hi + hj) - 3.f * dh * dh) /
             (12.f * r);
  }

  const float N = sidm_kernel_overlap_tophat_norm(h_min / h_max);

  return V_lens / (V_i * V_j * N);
}

/* ------------------------------------------------------------------------- */
#elif defined(SIDM_OVERLAP_CUBIC_SPLINE)

#if !defined(HYDRO_DIMENSION_3D)
#error "The cubic spline SIDM overlap kernel is only defined in 3D."
#endif

#define sidm_overlap_kernel_name "Cubic spline (M4)"

/*! Support of the cubic spline in units of h (H = gamma h), from table 1 of
 * Dehnen & Aly (2012). */
#define sidm_overlap_gamma_double 1.825742
#define sidm_overlap_gamma ((float)(sidm_overlap_gamma_double))

/*! Normalisation of the 3D cubic spline, W = C / gamma^3 P(u / gamma) */
#define sidm_overlap_constant_double (16. * M_1_PI)

/*! Number of table entries in s = h_min / h_max (over [0, 1]) */
#define sidm_overlap_table_n_s 128

/*! Number of table entries in x = r / (gamma (hi + hj)) (over [0, 1]) */
#define sidm_overlap_table_n_x 256

/*! Tabulated overlap F(x, s) for h_max = 1. Filled by
 * sidm_kernel_overlap_init(). Row-major in s. */
extern float
    sidm_overlap_table[sidm_overlap_table_n_s * sidm_overlap_table_n_x];

/*! Inverse of the normalisation 1 / N(s) of the overlap truncated at the
 * neighbour search radius. Filled by sidm_kernel_overlap_init(). */
extern float sidm_overlap_norm_inv[sidm_overlap_table_n_s];

/**
 * @brief Kernel overlap assuming cubic spline kernels.
 *
 * Returns \f$\Lambda_{ij} = \int W(|x - x_i|, h_i) W(|x - x_j|, h_j) d^3x\f$
 * with W the (normalised) cubic spline. The full overlap is non-zero for
 * r < gamma (hi + hj), beyond the search radius sidm_kernel_gamma h_max of
 * the force loop. It is therefore truncated at that radius and divided by
 * N(s), its integral over the search sphere, so that it integrates to 1.
 *
 * Lambda only depends on r / h_max and s = h_min / h_max, up to a factor
 * h_max^{-3}: \f$\Lambda_{ij} = F(x, s) / (N(s) h_{max}^3)\f$ with
 * x = r / (gamma (hi + hj)) in [0, 1]. F and 1 / N are tabulated by
 * sidm_kernel_overlap_init() (see kernel_sidm.c for the integrals) and read
 * here by bilinear interpolation.
 *
 * @param r Comoving distance between the two particles.
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 */
INLINE static float sidm_kernel_overlap(const float r, const float hi,
                                        const float hj) {

  const float h_max = max(hi, hj);

  /* Truncate at the neighbour search radius */
  if (r >= sidm_kernel_gamma * h_max) return 0.f;

  const float x = r / (sidm_overlap_gamma * (hi + hj));
  if (x >= 1.f) return 0.f;

  const float h_min = min(hi, hj);
  const float s = h_min / h_max;

  /* Locate the cell of the table */
  const float fs = s * (float)(sidm_overlap_table_n_s - 1);
  const float fx = x * (float)(sidm_overlap_table_n_x - 1);
  int is = (int)fs;
  int ix = (int)fx;
  if (is > sidm_overlap_table_n_s - 2) is = sidm_overlap_table_n_s - 2;
  if (ix > sidm_overlap_table_n_x - 2) ix = sidm_overlap_table_n_x - 2;
  const float ds = fs - (float)is;
  const float dx = fx - (float)ix;

  /* Bilinear interpolation */
  const float *row0 = &sidm_overlap_table[is * sidm_overlap_table_n_x];
  const float *row1 = row0 + sidm_overlap_table_n_x;
  const float F0 = (1.f - dx) * row0[ix] + dx * row0[ix + 1];
  const float F1 = (1.f - dx) * row1[ix] + dx * row1[ix + 1];
  const float F = (1.f - ds) * F0 + ds * F1;

  /* Normalisation of the truncated overlap */
  const float N_inv = (1.f - ds) * sidm_overlap_norm_inv[is] +
                      ds * sidm_overlap_norm_inv[is + 1];

  const float h_max_inv = 1.f / h_max;
  return F * N_inv * h_max_inv * h_max_inv * h_max_inv;
}

double sidm_kernel_overlap_exact(const double r, const double hi,
                                 const double hj);

/* ------------------------------------------------------------------------- */
#else

#error "An SIDM overlap kernel must be chosen at configure time !!"

/* ------------------------------------------------------------------------- */
#endif

void sidm_kernel_overlap_init(void);

#endif /* SWIFT_KERNEL_SIDM_H */
