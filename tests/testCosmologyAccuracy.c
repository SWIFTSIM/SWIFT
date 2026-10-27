/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/**
 * @file testCosmologyAccuracy.c
 * @brief Sign, monotonicity and accuracy checks of the cosmology module.
 *
 * The time-integration operators (drift, kicks, delta-time) and the
 * time-since-Big-Bang function are differences of linearly interpolated
 * tables. This test verifies that, for every ordered pair of inputs, these
 * differences are never negative because of rounding, that the quantities
 * updated by cosmology_update() are positive, that the inverse tables stay
 * within their domain, and it bounds the accuracy of the table-based
 * integrals against a direct quadrature of the same integrands.
 *
 * Several cosmologies are exercised, including the production-like
 * a_begin = 1/128, a run ending after a = 1, massive neutrinos and a
 * dynamical dark-energy equation of state.
 */

/* Some standard headers. */
#include <config.h>
#include <float.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <stdlib.h>
#include <time.h>

/* Includes. */
#include "swift.h"

/* Integrands defined in cosmology.c (not exposed in the header). */
double drift_integrand(double a, void *param);
double time_integrand(double a, void *param);

/* Number of entries in the interpolation tables (defined in cosmology.c) */
extern const int cosmology_table_length;

/*! Number of random samples per check */
#define N_RANDOM 200000

/*! Number of samples along the time-line in cosmology_update() sweeps */
#define N_SWEEP 20000

/*! Smallest time-bin for which the operators must be strictly positive.
 * Below this, integer times near the end of the run are not resolved by a
 * double (max_nr_timesteps = 2^57 > 2^53) and the operators can be 0. */
#define MIN_RESOLVED_BIN 8

/*! Tolerance on a(r(a)) round trip (measured 3e-8). */
#define TOL_DISTANCE_ROUND_TRIP 1e-7

/*! Relative tolerance on |t(a=1) - age today| and on |lookback(a=1)| */
#define TOL_AGE_AT_PRESENT 1e-9

/*! Tolerance on |r(a=1)| relative to r(a_begin) */
#define TOL_DISTANCE_AT_PRESENT 1e-9

/*! Relative accuracy of the table-based operators versus direct quadrature.
 * For a step shorter than one table bin the operator is the bin-averaged
 * integrand (measured 1.2e-4); it improves as the step grows. */
#define TOL_OPERATOR_SUB_BIN 3e-4
#define TOL_OPERATOR_BIN_44 2e-5
#define TOL_OPERATOR_BIN_52 2e-7

/*! Residual allowed for equal inputs, relative to the full-run value. The two
 * table lookups may be rounded differently when inlined (fast-math + LTO). */
#define TOL_EQUAL_INPUTS 1e-13

/*! Relative decrease tolerated between adjacent doubles (same reason) */
#define TOL_MONOTONIC 1e-14

/*! Tolerance on H(a=1)/H0 - 1 (closure of the density parameters) */
#define TOL_CLOSURE 1e-12

/* ------------------------------------------------------------------------ */

/* Uniform in [0, 1) */
static double random_unit(void) { return rand() / ((double)RAND_MAX + 1.); }

/* Uniform integer in [0, range) for ranges up to 2^62 */
static unsigned long long random_integer(const unsigned long long range) {
  const unsigned long long r =
      ((unsigned long long)rand() << 31) | (unsigned long long)rand();
  return r % range;
}

static gsl_integration_workspace *workspace;

/* Direct quadrature of a cosmology.c integrand with the same settings as
 * cosmology_init_tables() */
static double direct_integral(double (*f)(double, void *), struct cosmology *c,
                              const double a1, const double a2) {
  gsl_function F = {f, c};
  double result, abserr;
  const int status =
      gsl_integration_qag(&F, a1, a2, 0, 1.0e-10, 100000, GSL_INTEG_GAUSS61,
                          workspace, &result, &abserr);
  if (status) error("GSL quadrature failed (status=%d)", status);
  return result;
}

typedef double (*operator_fn)(const struct cosmology *, integertime_t,
                              integertime_t);

/* ------------------------------------------------------------------------ */

struct test_cosmology {
  const char *name;
  const char *params[16];

  /*! Tolerance on the a(t(a)) round trip. Limited by the uniform-in-t inverse
   * table at early times (measured 4e-5 at a_begin = 1/128, 2e-8 at 0.1). */
  double tol_time_round_trip;
};

static const struct test_cosmology cosmologies[] = {
    {"Planck-like, a = 1/128 -> 1",
     {"Cosmology:Omega_cdm:0.2589", "Cosmology:Omega_lambda:0.6910",
      "Cosmology:Omega_b:0.0486", "Cosmology:h:0.6774",
      "Cosmology:a_begin:0.0078125", "Cosmology:a_end:1.0", NULL},
     1e-4},
    {"Planck-like, a = 0.1 -> 1",
     {"Cosmology:Omega_cdm:0.2589", "Cosmology:Omega_lambda:0.6910",
      "Cosmology:Omega_b:0.0486", "Cosmology:h:0.6774", "Cosmology:a_begin:0.1",
      "Cosmology:a_end:1.0", NULL},
     1e-7},
    {"Planck-like, a = 1/128 -> 1.1",
     {"Cosmology:Omega_cdm:0.2589", "Cosmology:Omega_lambda:0.6910",
      "Cosmology:Omega_b:0.0486", "Cosmology:h:0.6774",
      "Cosmology:a_begin:0.0078125", "Cosmology:a_end:1.1", NULL},
     1e-4},
    {"Three massive neutrinos, a = 0.01 -> 1",
     {"Cosmology:Omega_cdm:0.241989232", "Cosmology:Omega_lambda:0.685365",
      "Cosmology:Omega_b:0.049199", "Cosmology:h:0.6737",
      "Cosmology:a_begin:1e-2", "Cosmology:a_end:1.0", "Cosmology:N_nu:3",
      "Cosmology:M_nu_eV:0.2,0.3,0.4", NULL},
     1e-4},
    {"w0-wa dark energy with radiation, a = 1/128 -> 1",
     {"Cosmology:Omega_cdm:0.2589", "Cosmology:Omega_lambda:0.6910",
      "Cosmology:Omega_b:0.0486", "Cosmology:h:0.6774",
      "Cosmology:Omega_r:8.5e-5", "Cosmology:w_0:-0.9", "Cosmology:w_a:0.1",
      "Cosmology:a_begin:0.0078125", "Cosmology:a_end:1.0", NULL},
     1e-4},
};

/* ------------------------------------------------------------------------ */

/**
 * @brief cosmology_get_delta_time_from_scale_factors() must be >= 0 for
 * a_start <= a_end and exactly 0 for a_start == a_end.
 */
static void test_delta_time_from_scale_factors(struct cosmology *c) {

  const double log_a_begin = c->log_a_begin;
  const double log_a_range = c->log_a_end - c->log_a_begin;
  const double delta_log_a = log_a_range / cosmology_table_length;

  /* Random ordered pairs */
  for (int k = 0; k < N_RANDOM; k++) {
    const double x1 = log_a_begin + random_unit() * log_a_range;
    const double x2 = log_a_begin + random_unit() * log_a_range;
    const double a1 = exp(min(x1, x2));
    const double a2 = exp(max(x1, x2));
    const double dt = cosmology_get_delta_time_from_scale_factors(c, a1, a2);
    if (dt < 0.)
      error("Negative delta time %e for a_start=%.17g a_end=%.17g", dt, a1, a2);
  }

  /* Adjacent doubles */
  for (int k = 0; k < N_RANDOM; k++) {
    const double a1 = exp(log_a_begin + random_unit() * log_a_range);
    const double a2 = nextafter(a1, 2.0);
    const double dt = cosmology_get_delta_time_from_scale_factors(c, a1, a2);
    if (dt < 0.)
      error("Negative delta time %e for adjacent a_start=%.17g", dt, a1);
  }

  /* Pairs straddling every table-bin edge, where the interpolation switches
   * between two different linear pieces */
  for (int i = 0; i <= cosmology_table_length; i++) {
    const double a_edge = exp(log_a_begin + i * delta_log_a);
    double lo = a_edge, hi = a_edge;
    for (int k = 0; k < 4; k++) {
      lo = nextafter(lo, 0.0);
      hi = nextafter(hi, 2.0);
      const double dt1 = cosmology_get_delta_time_from_scale_factors(c, lo, hi);
      const double dt2 =
          cosmology_get_delta_time_from_scale_factors(c, lo, a_edge);
      const double dt3 =
          cosmology_get_delta_time_from_scale_factors(c, a_edge, hi);
      if (dt1 < 0. || dt2 < 0. || dt3 < 0.)
        error("Negative delta time across table edge %d: %e %e %e", i, dt1, dt2,
              dt3);
    }
  }

  /* Equal inputs: zero up to a rounding residual of the table lookups */
  const double dt_full = c->time_end - c->time_begin;
  for (int k = 0; k < N_RANDOM / 10; k++) {
    const double a = exp(log_a_begin + random_unit() * log_a_range);
    const double dt = cosmology_get_delta_time_from_scale_factors(c, a, a);
    if (fabs(dt) > TOL_EQUAL_INPUTS * dt_full)
      error("Non-zero delta time %e for a_start == a_end", dt);
  }

  /* Monotonicity of the underlying time function */
  for (int k = 0; k < N_RANDOM; k++) {
    const double a1 = exp(log_a_begin + random_unit() * log_a_range);
    const double a2 = nextafter(a1, 2.0);
    if (cosmology_get_time_since_big_bang(c, a2) <
        cosmology_get_time_since_big_bang(c, a1) * (1. - TOL_MONOTONIC))
      error(
          "Time since Big Bang decreasing between a=%.17g and the next "
          "double",
          a1);
  }

  /* Stars born before a_begin: result must be finite and non-negative
   * (the accuracy of this extrapolation is not guaranteed, see below) */
  const double dt_early =
      cosmology_get_delta_time_from_scale_factors(c, 0.5 * c->a_begin, c->a);
  if (!(dt_early >= 0.) || !isfinite(dt_early))
    error("Invalid delta time %e for a_start < a_begin", dt_early);

  message("  delta_time_from_scale_factors: no negative values");
}

/**
 * @brief The integer-time operators must be >= 0, exactly 0 for equal times
 * and strictly positive for resolved time-bins.
 */
static void test_integer_time_operators(struct cosmology *c) {

  const operator_fn fns[5] = {
      cosmology_get_drift_factor, cosmology_get_grav_kick_factor,
      cosmology_get_hydro_kick_factor, cosmology_get_corr_kick_factor,
      cosmology_get_delta_time};
  const char *names[5] = {"drift", "grav_kick", "hydro_kick", "corr_kick",
                          "delta_time"};

  for (int bin = 0; bin <= num_time_bins; bin++) {
    const integertime_t step = 1LL << bin;
    const unsigned long long range =
        (unsigned long long)(max_nr_timesteps - step + 1);

    for (int k = 0; k < N_RANDOM / 20; k++) {
      const integertime_t ti_start = (integertime_t)random_integer(range);
      const integertime_t ti_end = ti_start + step;

      for (int f = 0; f < 5; f++) {
        const double v = fns[f](c, ti_start, ti_end);
        if (v < 0.)
          error("Negative %s factor %e for bin %d ti_start=%lld", names[f], v,
                bin, ti_start);
        if (bin >= MIN_RESOLVED_BIN && v == 0.)
          error("Zero %s factor for resolved bin %d ti_start=%lld", names[f],
                bin, ti_start);
      }
    }
  }

  /* Integer times straddling the table-bin edges */
  for (int i = 1; i < cosmology_table_length; i++) {
    const integertime_t ti_edge =
        (integertime_t)((double)i / cosmology_table_length *
                        (double)max_nr_timesteps);
    for (integertime_t d = -3; d <= 3; d++) {
      for (integertime_t e = d + 1; e <= 4; e++) {
        for (int f = 0; f < 5; f++) {
          const double v = fns[f](c, ti_edge + d, ti_edge + e);
          if (v < 0.)
            error("Negative %s factor %e across table edge %d", names[f], v, i);
        }
      }
    }
  }

  /* Equal times: zero up to a rounding residual of the table lookups */
  for (int k = 0; k < N_RANDOM / 20; k++) {
    const integertime_t ti =
        (integertime_t)random_integer((unsigned long long)max_nr_timesteps);
    for (int f = 0; f < 5; f++) {
      const double full = fns[f](c, 0, max_nr_timesteps);
      if (fabs(fns[f](c, ti, ti)) > TOL_EQUAL_INPUTS * full)
        error("Non-zero %s factor for ti_start == ti_end", names[f]);
    }
  }

  message(
      "  integer-time operators: no negative values, positive for bins "
      ">= %d",
      MIN_RESOLVED_BIN);
}

/**
 * @brief Accuracy of the table-based operators against direct quadrature.
 */
static void test_operator_accuracy(struct cosmology *c) {

  const int bins[3] = {30, 44, 52};
  const double tols[3] = {TOL_OPERATOR_SUB_BIN, TOL_OPERATOR_BIN_44,
                          TOL_OPERATOR_BIN_52};

  for (int b = 0; b < 3; b++) {
    const integertime_t step = 1LL << bins[b];
    const unsigned long long range =
        (unsigned long long)(max_nr_timesteps - step + 1);
    double max_err_drift = 0., max_err_time = 0.;

    for (int k = 0; k < 100; k++) {
      const integertime_t ti_start = (integertime_t)random_integer(range);
      const integertime_t ti_end = ti_start + step;
      const double a1 = c->a_begin * exp(ti_start * c->time_base);
      const double a2 = c->a_begin * exp(ti_end * c->time_base);

      const double drift_ref = direct_integral(drift_integrand, c, a1, a2);
      const double time_ref = direct_integral(time_integrand, c, a1, a2);
      const double err_drift = fabs(
          cosmology_get_drift_factor(c, ti_start, ti_end) / drift_ref - 1.);
      const double err_time =
          fabs(cosmology_get_delta_time(c, ti_start, ti_end) / time_ref - 1.);
      max_err_drift = max(max_err_drift, err_drift);
      max_err_time = max(max_err_time, err_time);
    }

    message("  bin %2d: max rel. error drift=%.2e delta_time=%.2e (tol %.0e)",
            bins[b], max_err_drift, max_err_time, tols[b]);
    if (max_err_drift > tols[b] || max_err_time > tols[b])
      error("Operator accuracy worse than %e for bin %d", tols[b], bins[b]);
  }
}

/**
 * @brief Quantities set by cosmology_update() along the whole time-line.
 */
static void test_update_sweep(struct cosmology *c,
                              const struct phys_const *phys_const) {

  double time_prev = -DBL_MAX;

  for (int k = 0; k <= N_SWEEP; k++) {
    const integertime_t ti =
        (integertime_t)((double)max_nr_timesteps * k / N_SWEEP);
    cosmology_update(c, phys_const, ti);

    if (!(c->a > 0.) || !(c->H > 0.) || !(c->critical_density > 0.) ||
        !(c->a_dot > 0.) || !(c->mean_density > 0.))
      error("Non-positive a/H/rho_crit/a_dot/rho_mean at ti=%lld", ti);
    if (!(c->time >= 0.)) error("Negative time %e at ti=%lld", c->time, ti);
    if (!(c->Omega_nu >= 0.)) error("Negative Omega_nu at ti=%lld", ti);
    if (!(c->overdensity_BN98 > 0.))
      error("Non-positive BN98 overdensity at ti=%lld", ti);
    if (c->time < time_prev)
      error("Time decreasing along the time-line at ti=%lld", ti);
    time_prev = c->time;
  }

  /* At ti = 0 we must recover the start of the run exactly */
  cosmology_update(c, phys_const, 0);
  if (c->a != c->a_begin || c->time != c->time_begin)
    error("cosmology_update(ti=0) does not recover a_begin / time_begin");

  /* Checks at a = 1 if the run reaches it */
  if (c->a_end >= 1.) {
    const integertime_t ti_1 =
        (integertime_t)llround(-c->log_a_begin / c->time_base);
    cosmology_update(c, phys_const, ti_1);
    if (fabs(c->a - 1.) > 1e-12)
      error("Could not land on a=1 (got %.17g)", c->a);

    const double closure = fabs(c->H / c->H0 - 1.);
    const double lookback =
        fabs(c->lookback_time) / c->universe_age_at_present_day;
    const double age_err = fabs(cosmology_get_time_since_big_bang(c, 1.) -
                                c->universe_age_at_present_day) /
                           c->universe_age_at_present_day;
    const double dist_err = fabs(cosmology_get_comoving_distance(c, 1.)) /
                            cosmology_get_comoving_distance(c, c->a_begin);

    message(
        "  at a=1: |H/H0-1|=%.1e |lookback|/age=%.1e |t(1)-age|/age=%.1e "
        "|r(1)|/r(a_begin)=%.1e",
        closure, lookback, age_err, dist_err);
    if (closure > TOL_CLOSURE) error("Density parameters do not close");
    if (lookback > TOL_AGE_AT_PRESENT || age_err > TOL_AGE_AT_PRESENT)
      error("Age at a=1 inconsistent with universe_age_at_present_day");
    if (dist_err > TOL_DISTANCE_AT_PRESENT)
      error("Comoving distance at a=1 is not zero");
  }

  message("  cosmology_update sweep: all quantities positive, time monotonic");
}

/**
 * @brief Inverse tables: domain and round-trip accuracy.
 */
static void test_inverse_tables(struct cosmology *c,
                                const double tol_time_round_trip) {

  double max_err_time = 0., max_err_dist = 0.;
  double r_prev = DBL_MAX;

  for (int k = 0; k <= N_SWEEP; k++) {
    double a = exp(c->log_a_begin +
                   (c->log_a_end - c->log_a_begin) * k / (double)N_SWEEP);
    a = fmin(fmax(a, c->a_begin), c->a_end);

    /* a(t(a)) */
    const double t = cosmology_get_time_since_big_bang(c, a);
    const double a_t = cosmology_get_scale_factor(c, t);
    if (a_t < c->a_begin) error("a(t) below a_begin at a=%.17g", a);
    max_err_time = max(max_err_time, fabs(a_t / a - 1.));

    /* a(r(a)) */
    const double r = cosmology_get_comoving_distance(c, a);
    if (r > r_prev) error("Comoving distance increasing with a at a=%.17g", a);
    r_prev = r;
    if (a < 1. && r < 0.)
      error("Negative comoving distance %e at a=%.17g < 1", r, a);
    const double a_r = cosmology_scale_factor_at_comoving_distance(c, r);
    if (a_r < c->a_begin) error("a(r) below a_begin at a=%.17g", a);
    max_err_dist = max(max_err_dist, fabs(a_r / a - 1.));
  }

  message("  inverse tables: max round-trip error a(t(a))=%.2e a(r(a))=%.2e",
          max_err_time, max_err_dist);
  if (max_err_time > tol_time_round_trip)
    error("a(t(a)) round trip worse than %e", tol_time_round_trip);
  if (max_err_dist > TOL_DISTANCE_ROUND_TRIP)
    error("a(r(a)) round trip worse than %e", TOL_DISTANCE_ROUND_TRIP);
}

/**
 * @brief Neutrino density parameter: sign and closure.
 */
static void test_neutrinos(struct cosmology *c) {

  if (c->N_nu == 0) return;

  for (int k = 0; k <= N_SWEEP; k++) {
    const double a =
        exp(log(1e-9) + (c->log_a_end - log(1e-9)) * k / (double)N_SWEEP);
    const double Omega_nu = cosmology_get_neutrino_density(c, a);
    if (!(Omega_nu >= 0.)) error("Negative Omega_nu(a=%e)", a);
  }
  if (!(c->Omega_nu_0 > 0.)) error("Non-positive Omega_nu_0");

  const double sum = c->Omega_cdm + c->Omega_b + c->Omega_r + c->Omega_lambda +
                     c->Omega_nu_0 + c->Omega_k;
  if (fabs(sum - 1.) > 1e-14) error("Density parameters do not sum to 1");

  message("  neutrinos: Omega_nu(a) >= 0 and parameters close");
}

/* ------------------------------------------------------------------------ */

int main(int argc, char *argv[]) {

  /* Seed the random number generator with the time, or with the seed given
   * on the command line to reproduce an earlier run */
  const unsigned int seed = (argc > 1)
                                ? (unsigned int)strtoul(argv[1], NULL, 10)
                                : (unsigned int)time(NULL);
  srand(seed);
  message("Random seed: %u (pass it as argument to reproduce this run)", seed);

  workspace = gsl_integration_workspace_alloc(100000);

  /* initialization of unit system and constants */
  struct unit_system us;
  units_init_cgs(&us);

  const int n_cosmo = sizeof(cosmologies) / sizeof(cosmologies[0]);
  for (int i = 0; i < n_cosmo; i++) {

    message("Cosmology %d/%d: %s", i + 1, n_cosmo, cosmologies[i].name);

    struct swift_params params;
    parser_init("", &params);
    for (int j = 0; cosmologies[i].params[j] != NULL; j++)
      parser_set_param(&params, cosmologies[i].params[j]);

    struct phys_const phys_const;
    phys_const_init(&us, &params, &phys_const);

    struct cosmology cosmo;
    cosmology_init(&params, &us, &phys_const, &cosmo);

    test_delta_time_from_scale_factors(&cosmo);
    test_integer_time_operators(&cosmo);
    test_operator_accuracy(&cosmo);
    test_update_sweep(&cosmo, &phys_const);
    test_inverse_tables(&cosmo, cosmologies[i].tol_time_round_trip);
    test_neutrinos(&cosmo);

    cosmology_clean(&cosmo);
  }

  gsl_integration_workspace_free(workspace);

  message(
      "Everything seems fine with the accuracy of the cosmology functions.");
  return 0;
}
