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

#include <config.h>

#include "kernel_sidm.h"
#include <math.h>

#if defined(SIDM_OVERLAP_CUBIC_SPLINE)

/* The overlap table (see kernel_sidm.h) */
float sidm_overlap_table[sidm_overlap_table_n_s * sidm_overlap_table_n_x];

/* The normalisation table of the truncated overlap (see kernel_sidm.h) */
float sidm_overlap_norm_inv[sidm_overlap_table_n_s];

/*! Number of Simpson intervals used to compute each table entry (even).
 * The integrands have kinks (the spline branches at u = gamma / 2 and
 * |r - s t| = 0), which limits Simpson's rule to second order there; with
 * 512 intervals Lambda agrees with direct 3D quadrature to ~1e-10. */
#define sidm_overlap_n_integration 512

#if (sidm_overlap_n_integration % 2) != 0
#error "Simpson's rule needs an even number of intervals."
#endif

/**
 * @brief Weight of point k of an N-interval composite Simpson rule (x 3/dx).
 */
static double simpson_weight(const int k, const int N) {
  if (k == 0 || k == N) return 1.;
  return (k % 2) ? 4. : 2.;
}

/**
 * @brief Cubic spline kernel W(u) for h = 1 (3D).
 *
 * W(u) = C / gamma^3 P(u / gamma), with P(x) = 3 x^3 - 3 x^2 + 1/2 for
 * x < 1/2, (1 - x)^3 for 1/2 <= x < 1 and 0 otherwise (Dehnen & Aly 2012).
 *
 * @param u The distance in units of h.
 */
static double sidm_overlap_w(const double u) {

  const double gamma = sidm_overlap_gamma_double;
  const double x = u / gamma;

  double P;
  if (x < 0.5)
    P = 3. * x * x * x - 3. * x * x + 0.5;
  else if (x < 1.)
    P = (1. - x) * (1. - x) * (1. - x);
  else
    P = 0.;

  return sidm_overlap_constant_double * P / (gamma * gamma * gamma);
}

/**
 * @brief G(b) = int_0^b u W(u) du for the cubic spline with h = 1 (3D).
 *
 * Substituting u = gamma y gives G(b) = C / gamma Q(b / gamma), with
 *
 *   Q(y) = int_0^y x P(x) dx
 *        = 3 y^5 / 5 - 3 y^4 / 4 + y^2 / 4                for y < 1/2,
 *        = Q(1/2) + R(y) - R(1/2)                        for 1/2 <= y < 1,
 *
 * where R(x) = x^2 / 2 - x^3 + 3 x^4 / 4 - x^5 / 5 is a primitive of
 * x (1 - x)^3 and Q(1/2) = 11 / 320. G is constant for b >= gamma.
 *
 * @param b The upper limit in units of h.
 */
static double sidm_overlap_G(const double b) {

  const double gamma = sidm_overlap_gamma_double;
  const double y = fmin(b / gamma, 1.);

  double Q;
  if (y < 0.5) {
    const double y2 = y * y;
    Q = 0.6 * y2 * y2 * y - 0.75 * y2 * y2 + 0.25 * y2;
  } else {
#define R(x) ((x) * (x) * (0.5 + (x) * (-1. + (x) * (0.75 - 0.2 * (x)))))
    Q = 11. / 320. + R(y) - R(0.5);
#undef R
  }

  return sidm_overlap_constant_double * Q / gamma;
}

/**
 * @brief Overlap integral Lambda(r, s) of two cubic splines with h_max = 1.
 *
 * Particle j has the smaller kernel, h_j = s. Integrate over spherical shells
 * of radius rho = s t around it (t = distance in units of h_j): the shell
 * holds the weight W_j = W(t) / s^3 over a volume 4 pi s^3 t^2 dt, and the
 * other kernel, W(|x - x_i|) with h_i = 1, averages over that shell to
 *
 *   <W>_shell = [G(r + rho) - G(|r - rho|)] / (2 r rho)
 *
 * (spherical mean of a radial function, from the law of cosines). So
 *
 *   Lambda(r, s) = 2 pi / (r s) int_0^gamma t W(t) [G(r + s t) - G(|r - s t|)] dt.
 *
 * Limits: for r -> 0 the shell mean is W(s t), giving
 * Lambda = 4 pi int_0^gamma t^2 W(t) W(s t) dt; for s -> 0 kernel j is a
 * delta function and Lambda = W(r).
 *
 * @param r The distance in units of h_max.
 * @param s The ratio h_min / h_max (in [0, 1]).
 */
static double sidm_overlap_integral(const double r, const double s) {

  if (s == 0.) return sidm_overlap_w(r);

  const int N = sidm_overlap_n_integration;
  const double dt = sidm_overlap_gamma_double / N;
  double integral = 0.;

  for (int k = 0; k <= N; k++) {
    const double t = k * dt;

    double integrand;
    if (r == 0.)
      integrand = t * t * sidm_overlap_w(t) * sidm_overlap_w(s * t);
    else
      integrand = t * sidm_overlap_w(t) *
                  (sidm_overlap_G(r + s * t) - sidm_overlap_G(fabs(r - s * t)));

    integral += simpson_weight(k, N) * integrand;
  }
  integral *= dt / 3.;

  if (r == 0.) return 4. * M_PI * integral;
  return 2. * M_PI * integral / (r * s);
}

/**
 * @brief Normalisation N(s) of the truncated overlap (h_max = 1).
 *
 * N(s) = int_0^r_max 4 pi r^2 Lambda(r, s) dr, with r_max the smaller of the
 * search radius sidm_kernel_gamma and the support gamma (1 + s) of Lambda.
 * N = 1 when the search sphere contains the whole support.
 *
 * @param s The ratio h_min / h_max (in [0, 1]).
 */
static double sidm_overlap_norm(const double s) {

  const double r_max =
      fmin((double)sidm_kernel_gamma, sidm_overlap_gamma_double * (1. + s));
  const int N = sidm_overlap_n_integration;
  const double dr = r_max / N;
  double integral = 0.;

  for (int k = 0; k <= N; k++) {
    const double r = k * dr;
    integral += simpson_weight(k, N) * r * r * sidm_overlap_integral(r, s);
  }

  return 4. * M_PI * integral * dr / 3.;
}

/**
 * @brief Overlap of two cubic spline kernels, computed directly (no table).
 *
 * Truncated at the search radius and normalised as sidm_kernel_overlap().
 * Slow; used to build and test the table.
 *
 * @param r Distance between the two particles.
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 */
double sidm_kernel_overlap_exact(const double r, const double hi,
                                 const double hj) {

  const double h_max = fmax(hi, hj);
  if (r >= sidm_kernel_gamma * h_max) return 0.;

  const double s = fmin(hi, hj) / h_max;
  return sidm_overlap_integral(r / h_max, s) /
         (sidm_overlap_norm(s) * h_max * h_max * h_max);
}

/**
 * @brief Fill the cubic spline overlap table.
 *
 * Entry (is, ix) holds F(x, s) with s = is / (n_s - 1) and
 * x = ix / (n_x - 1) = r / (gamma (1 + s)) for h_max = 1.
 */
void sidm_kernel_overlap_init(void) {

  static int initialised = 0;
  if (initialised) return;

  for (int is = 0; is < sidm_overlap_table_n_s; is++) {
    const double s = (double)is / (sidm_overlap_table_n_s - 1);
    sidm_overlap_norm_inv[is] = 1. / sidm_overlap_norm(s);
    for (int ix = 0; ix < sidm_overlap_table_n_x; ix++) {
      const double x = (double)ix / (sidm_overlap_table_n_x - 1);
      const double r = x * sidm_overlap_gamma_double * (1. + s);
      sidm_overlap_table[is * sidm_overlap_table_n_x + ix] =
          (ix == sidm_overlap_table_n_x - 1) ? 0.f
                                             : sidm_overlap_integral(r, s);
    }
  }

  initialised = 1;
}

#else

/**
 * @brief Nothing to initialise for this overlap kernel.
 */
void sidm_kernel_overlap_init(void) {}

#endif
