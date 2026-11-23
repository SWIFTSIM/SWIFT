/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/*
 * Vector weights of the GEAR mechanical feedback (Hopkins et al. 2018, Appendix
 * E). For a star with gas on both sides of every axis, the normalized weights
 * w_j_bar of its neighbours obey sum_j |w_j_bar| = 1 (they share the star's
 * budget) and sum_j w_j_bar = 0 (the momenta they carry cancel). Both hold up
 * to rounding. The normalization sum_j |w_j| is stored in single precision
 * (#feedback_spart_data.enrichment_weight), so sum_j |w_j_bar| = 1 holds to
 * the float epsilon times the number of terms. The vector sum divides every
 * term by the same normalization, so it holds to the double epsilon.
 */

#ifdef FEEDBACK_GEAR_MECHANICAL

#define N_GAS 200

/* Rounding bar of sum_j |w_j_bar| - 1: N_GAS float roundings. */
static const double norm_sum_bar = N_GAS * FLT_EPSILON;

/* Rounding bar of |sum_j w_j_bar|: N_GAS double roundings of terms up to 1. */
static const double vector_sum_bar = 64. * N_GAS * DBL_EPSILON;

/**
 * @brief Check the two sums of a set of normalized weights.
 *
 * @param w The normalized weights, n x 3.
 * @param n The number of weights.
 * @param norm_err (return) |sum_j |w_j| - 1|.
 * @param vec_err (return) |sum_j w_j|.
 * @return 1 if a sum is not finite or above the bar, 0 otherwise.
 */
static int check_weight_sums(const double (*w)[3], const int n,
                             double *norm_err, double *vec_err) {
  double norm_sum = 0.;
  double vec_sum[3] = {0., 0., 0.};
  for (int j = 0; j < n; j++) {
    norm_sum += sqrt(w[j][0] * w[j][0] + w[j][1] * w[j][1] + w[j][2] * w[j][2]);
    for (int k = 0; k < 3; k++) vec_sum[k] += w[j][k];
  }
  *norm_err = fabs(norm_sum - 1.);
  *vec_err = sqrt(vec_sum[0] * vec_sum[0] + vec_sum[1] * vec_sum[1] +
                  vec_sum[2] * vec_sum[2]);
  /* A NaN fails every comparison, so test finiteness first. */
  if (!(*norm_err <= DBL_MAX) || !(*vec_err <= DBL_MAX)) return 1;
  return *norm_err > norm_sum_bar || *vec_err > vector_sum_bar;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  unsigned int seed = 12345;
  const float h_star = 1.f;

  static struct part parts[N_GAS];
  struct spart star;
  memset(parts, 0, sizeof(parts));
  memset(&star, 0, sizeof(star));
  star.h = h_star;
  star.density.wcount = 30.f;

  /* Gas inside the kernel of the star, on both sides of every axis. */
  float dx[N_GAS][3];
  for (int j = 0; j < N_GAS; j++) {
    float r2;
    do {
      for (int k = 0; k < 3; k++)
        dx[j][k] = (2.f * rand_r(&seed) / (float)RAND_MAX - 1.f) * h_star *
                   kernel_gamma;
      r2 = dx[j][0] * dx[j][0] + dx[j][1] * dx[j][1] + dx[j][2] * dx[j][2];
    } while (r2 >= h_star * h_star * kernel_gamma2 || r2 < 1e-4f);
    parts[j].h = 0.5f + rand_r(&seed) / (float)RAND_MAX;
    parts[j].feedback_data.density.wcount =
        20.f + 20.f * rand_r(&seed) / (float)RAND_MAX;
  }

  /* Loop 2 of the module: the sums of the f_plus and f_minus factors. */
  for (int j = 0; j < N_GAS; j++) {
    const float r2 =
        dx[j][0] * dx[j][0] + dx[j][1] * dx[j][1] + dx[j][2] * dx[j][2];
    double plus[3], minus[3], omega;
    feedback_compute_scalar_weight(r2, dx[j], h_star, parts[j].h, &star,
                                   &parts[j], plus, minus, &omega);
    for (int k = 0; k < 3; k++) {
      star.feedback_data.f_sum_plus_term[k] += omega * fabs(plus[k]);
      star.feedback_data.f_sum_minus_term[k] += omega * fabs(minus[k]);
    }
  }

  /* Loop 3: the normalization, sum_j |w_j|. */
  for (int j = 0; j < N_GAS; j++) {
    const float r2 =
        dx[j][0] * dx[j][0] + dx[j][1] * dx[j][1] + dx[j][2] * dx[j][2];
    double f_plus[3], f_minus[3], w[3];
    feedback_compute_vector_weight_non_normalized(
        r2, dx[j], h_star, parts[j].h, &star, &parts[j], f_plus, f_minus, w);
    star.feedback_data.enrichment_weight +=
        sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
  }

  /* Feedback loop: the normalized weights. */
  static double w_bar[N_GAS][3];
  for (int j = 0; j < N_GAS; j++) {
    const float r2 =
        dx[j][0] * dx[j][0] + dx[j][1] * dx[j][1] + dx[j][2] * dx[j][2];
    feedback_compute_vector_weight_normalized(r2, dx[j], h_star, parts[j].h,
                                              &star, &parts[j], w_bar[j]);
  }

  int fail = 0;
  double norm_err, vec_err;
  const int bad =
      check_weight_sums((const double (*)[3])w_bar, N_GAS, &norm_err, &vec_err);
  printf(
      "weights: |sum |w_bar| - 1| = %.3e (bar %.3e), |sum w_bar| = %.3e "
      "(bar %.3e)\n",
      norm_err, norm_sum_bar, vec_err, vector_sum_bar);
  if (bad) {
    printf("FAILED: the weight sums exceed the rounding bar\n");
    fail = 1;
  }

  /* Control: one weight changed by 1e-6 of its size must fail the check. */
  static double w_bad[N_GAS][3];
  memcpy(w_bad, w_bar, sizeof(w_bar));
  for (int k = 0; k < 3; k++) w_bad[7][k] *= 1. + 1e-6;
  if (!check_weight_sums((const double (*)[3])w_bad, N_GAS, &norm_err,
                         &vec_err)) {
    printf("FAILED: the check passes a weight changed by 1e-6\n");
    fail = 1;
  }

  /* Control: all weights scaled by 1 + 1e-3 must fail the norm check. */
  for (int j = 0; j < N_GAS; j++)
    for (int k = 0; k < 3; k++) w_bad[j][k] = w_bar[j][k] * (1. + 1e-3);
  if (!check_weight_sums((const double (*)[3])w_bad, N_GAS, &norm_err,
                         &vec_err)) {
    printf("FAILED: the check passes weights scaled by 1 + 1e-3\n");
    fail = 1;
  }

  if (!fail) printf("PASSED\n");
  return fail;
}

#else

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  printf("skip: GEAR mechanical feedback is not configured\n");
  return 77;
}

#endif /* FEEDBACK_GEAR_MECHANICAL */
