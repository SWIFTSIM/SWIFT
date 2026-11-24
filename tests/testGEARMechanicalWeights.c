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
#include "feedback_iact.h"
#include "swift.h"

/*
 * Vector weights of the GEAR mechanical feedback (Hopkins et al. 2018, Appendix
 * E), computed by the module's own prep loops. Gas on the six half axes at one
 * distance and smoothing length must get w_j_bar = +/- e_k / 6 by symmetry, and
 * gas on three positive half axes e_k / 3. Gas on the star, beyond the maximal
 * radius or outside both kernels gets nothing and leaves the other weights
 * unchanged. For random gas, sum_j |w_j_bar| = 1 and sum_j w_j_bar = 0.
 */

#ifdef FEEDBACK_GEAR_MECHANICAL

#define N_MAX 200

/**
 * @brief Gas neighbours of one star.
 */
struct neighbours {
  int n;
  float dx[N_MAX][3];
  struct part parts[N_MAX];
};

static struct feedback_props props;
static struct cosmology cosmo;

/**
 * @brief Add a gas neighbour at the given position relative to the star.
 */
static void add_gas(struct neighbours *ng, const float x, const float y,
                    const float z, const float h) {
  const int j = ng->n++;
  memset(&ng->parts[j], 0, sizeof(struct part));
  ng->dx[j][0] = -x;
  ng->dx[j][1] = -y;
  ng->dx[j][2] = -z;
  ng->parts[j].h = h;
  ng->parts[j].feedback_data.density.wcount = 30.f;
}

/**
 * @brief Run the prep loops 2 and 3 of the module over the neighbours, then
 * the normalized weights as the feedback loop computes them (zero for gas the
 * feedback loop skips).
 */
static void compute_weights(struct spart *si, struct neighbours *ng,
                            double (*w_bar)[3]) {
  memset(si->feedback_data.f_sum_plus_term, 0, 3 * sizeof(double));
  memset(si->feedback_data.f_sum_minus_term, 0, 3 * sizeof(double));
  si->feedback_data.enrichment_weight = 0.f;

  for (int j = 0; j < ng->n; j++) {
    const float *dx = ng->dx[j];
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    runner_iact_nonsym_feedback_prep2(r2, dx, si->h, ng->parts[j].h, si,
                                      &ng->parts[j], NULL, &cosmo, &props, 0);
  }
  for (int j = 0; j < ng->n; j++) {
    const float *dx = ng->dx[j];
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    runner_iact_nonsym_feedback_prep3(r2, dx, si->h, ng->parts[j].h, si,
                                      &ng->parts[j], NULL, &cosmo, &props, 0);
  }
  for (int j = 0; j < ng->n; j++) {
    const float *dx = ng->dx[j];
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    w_bar[j][0] = w_bar[j][1] = w_bar[j][2] = 0.;
    if (r2 > props.r_max * props.r_max) continue;
    if (si->feedback_data.enrichment_weight <= 0.f) continue;
    feedback_compute_vector_weight_normalized(r2, dx, si->h, ng->parts[j].h, si,
                                              &ng->parts[j], w_bar[j]);
  }
}

/**
 * @brief Compare the weights with a reference, component by component.
 *
 * @return The number of components further than @p bar from the reference,
 * or not finite.
 */
static int check_reference(const char *name, const double (*w)[3],
                           const double (*ref)[3], const int n,
                           const double bar) {
  int bad = 0;
  double worst = 0.;
  for (int j = 0; j < n; j++)
    for (int k = 0; k < 3; k++) {
      const double d = fabs(w[j][k] - ref[j][k]);
      /* A NaN fails every comparison, so test finiteness first. */
      if (!(d <= DBL_MAX) || d > bar) bad++;
      if (d > worst) worst = d;
    }
  printf("%s: largest difference to the reference %.3e (bar %.3e)%s\n", name,
         worst, bar, bad ? " FAILED" : "");
  return bad;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  memset(&props, 0, sizeof(props));
  memset(&cosmo, 0, sizeof(cosmo));
  props.r_max = 10.f;

  struct spart star;
  memset(&star, 0, sizeof(star));
  star.h = 1.f;
  star.density.wcount = 30.f;
  star.feedback_data.supernovae.energy_ejected = 1.f;
  star.feedback_data.supernovae.mass_ejected = 1.f;

  static struct neighbours ng;
  static double w[N_MAX][3], ref[N_MAX][3];
  int fail = 0;

  /* The normalization sum_j |w_j| is a float: up to N roundings of
     FLT_EPSILON relative, twice for margin, times the weight 1/N. */
  const double six_bar = 2. * 6. * FLT_EPSILON / 6.;
  const double three_bar = 2. * 3. * FLT_EPSILON / 3.;

  /* Six half axes. */
  ng.n = 0;
  const float d = 0.5f;
  for (int k = 0; k < 3; k++) {
    float x[3] = {0.f, 0.f, 0.f};
    x[k] = d;
    add_gas(&ng, x[0], x[1], x[2], 1.f);
    x[k] = -d;
    add_gas(&ng, x[0], x[1], x[2], 1.f);
  }
  memset(ref, 0, sizeof(ref));
  for (int k = 0; k < 3; k++) {
    ref[2 * k][k] = 1. / 6.;
    ref[2 * k + 1][k] = -1. / 6.;
  }
  compute_weights(&star, &ng, w);
  fail += check_reference("six half axes", (const double (*)[3])w,
                          (const double (*)[3])ref, ng.n, six_bar);

  /* The same, plus gas on the star, beyond the maximal radius and outside
     both kernels: these get nothing and change nothing. */
  add_gas(&ng, 0.f, 0.f, 0.f, 1.f);
  add_gas(&ng, 11.f, 0.f, 0.f, 1.f);
  add_gas(&ng, 0.f, 3.f, 0.f, 1.f);
  compute_weights(&star, &ng, w);
  fail += check_reference("six half axes and three empty neighbours",
                          (const double (*)[3])w, (const double (*)[3])ref,
                          ng.n, six_bar);

  /* One-sided gas: three positive half axes. */
  ng.n = 0;
  memset(ref, 0, sizeof(ref));
  for (int k = 0; k < 3; k++) {
    float x[3] = {0.f, 0.f, 0.f};
    x[k] = d;
    add_gas(&ng, x[0], x[1], x[2], 1.f);
    ref[k][k] = 1. / 3.;
  }
  compute_weights(&star, &ng, w);
  fail += check_reference("three positive half axes", (const double (*)[3])w,
                          (const double (*)[3])ref, ng.n, three_bar);

  /* Random gas inside the kernel of the star. */
  ng.n = 0;
  unsigned int seed = 12345;
  while (ng.n < N_MAX) {
    float x[3];
    for (int k = 0; k < 3; k++)
      x[k] = (2.f * rand_r(&seed) / (float)RAND_MAX - 1.f) * kernel_gamma;
    const float r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    if (r2 >= kernel_gamma2 || r2 < 1e-4f) continue;
    add_gas(&ng, x[0], x[1], x[2], 0.5f + rand_r(&seed) / (float)RAND_MAX);
    ng.parts[ng.n - 1].feedback_data.density.wcount =
        20.f + 20.f * rand_r(&seed) / (float)RAND_MAX;
  }
  compute_weights(&star, &ng, w);
  double norm_sum = 0., vec[3] = {0., 0., 0.};
  for (int j = 0; j < ng.n; j++) {
    norm_sum += sqrt(w[j][0] * w[j][0] + w[j][1] * w[j][1] + w[j][2] * w[j][2]);
    for (int k = 0; k < 3; k++) vec[k] += w[j][k];
  }
  const double norm_err = fabs(norm_sum - 1.);
  const double vec_err =
      sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  /* The normalization is a float: N roundings of FLT_EPSILON. */
  const double norm_bar = N_MAX * FLT_EPSILON;
  /* Per axis, f+ S+ and f- S- agree to about 8 double roundings (products and
     square roots of the f factors); the sum of N terms of total size at most
     sum_j |w_j_bar| = 1 adds N roundings; sqrt(3) for the norm. */
  const double vec_bar = sqrt(3.) * (N_MAX + 8) * DBL_EPSILON;
  printf(
      "random gas: |sum |w_bar| - 1| = %.3e (bar %.3e), |sum w_bar| = %.3e "
      "(bar %.3e)\n",
      norm_err, norm_bar, vec_err, vec_bar);
  if (!(norm_err <= norm_bar) || !(vec_err <= vec_bar)) {
    printf("FAILED: random gas\n");
    fail++;
  }

#if FEEDBACK_GEAR_MECHANICAL_MODE == 2
  /* A star with a supernova whose normalization is zero: loop 4 must leave its
     accumulators exactly zero, without a 1/0 before its zero tests. */
  ng.n = 0;
  add_gas(&ng, d, 0.f, 0.f, 1.f);
  hydro_set_mass(&ng.parts[0], 1.f);
  ng.parts[0].v[0] = 1.f;
  compute_weights(&star, &ng, w);
  star.feedback_data.enrichment_weight = 0.f;
  memset(&star.feedback_data.accumulator_sn, 0,
         sizeof(star.feedback_data.accumulator_sn));
  memset(&star.feedback_data.accumulator_winds, 0,
         sizeof(star.feedback_data.accumulator_winds));
  runner_iact_nonsym_feedback_prep4(d * d, ng.dx[0], star.h, ng.parts[0].h,
                                    &star, &ng.parts[0], NULL, &cosmo, &props,
                                    0);
  static const char zeros[sizeof(star.feedback_data.accumulator_sn)];
  if (memcmp(&star.feedback_data.accumulator_sn, zeros, sizeof(zeros)) ||
      memcmp(&star.feedback_data.accumulator_winds, zeros, sizeof(zeros))) {
    printf(
        "FAILED: loop 4 changed the accumulators of a star with no weight\n");
    fail++;
  } else {
    printf("loop 4 with no weight: accumulators unchanged\n");
  }
#endif

  if (!fail) printf("PASSED\n");
  return fail ? 1 : 0;
}

#else

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  printf("skip: GEAR mechanical feedback is not configured\n");
  return 77;
}

#endif /* FEEDBACK_GEAR_MECHANICAL */
