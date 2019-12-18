/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

/* Some standard headers */
#include <fenv.h>
#include <math.h>

/* Local headers. */
#include "swift.h"

#ifdef HAVE_LIBGSL
#include <gsl/gsl_integration.h>
#endif

static INLINE double w_tilde(double a, double w0, double wa) {
  return (a - 1.) * wa - (1. + w0 + wa) * log(a);
}

static INLINE double E(double Or, double Om, double Ok, double Ol, double w0,
                       double wa, double a) {
  const double a_inv = 1. / a;
  return sqrt(Or * a_inv * a_inv * a_inv * a_inv + Om * a_inv * a_inv * a_inv +
              Ok * a_inv * a_inv + Ol * exp(3. * w_tilde(a, w0, wa)));
}

static INLINE double drift_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double a_inv = 1. / a;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;
  return (1. / H) * a_inv * a_inv * a_inv;
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Fake parameter file with Planck+13 cosmology
   * and the usual cosmological system of units */
  struct swift_params *params =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  parser_set_param(params, "InternalUnitSystem:UnitMass_in_cgs:1.98848e43");
  parser_set_param(params,
                   "InternalUnitSystem:UnitLength_in_cgs:3.08567758e24");
  parser_set_param(params, "InternalUnitSystem:UnitVelocity_in_cgs:1e5");
  parser_set_param(params, "InternalUnitSystem:UnitCurrent_in_cgs:1.");
  parser_set_param(params, "InternalUnitSystem:UnitTemp_in_cgs:1.");
  parser_set_param(params, "Cosmology:Omega_m:0.307");
  parser_set_param(params, "Cosmology:Omega_lambda:0.693");
  parser_set_param(params, "Cosmology:Omega_b:0.0482519");
  parser_set_param(params, "Cosmology:h:0.6777");
  parser_set_param(params, "Cosmology:a_begin:0.0078125");
  parser_set_param(params, "Cosmology:a_end:1.");

  /* Initialise everything */
  struct unit_system units;
  units_init_from_params(&units, params, "InternalUnitSystem");

  struct phys_const phys_const;
  phys_const_init(&units, params, &phys_const);

  struct cosmology cosmo;
  cosmology_init(params, &units, &phys_const, &cosmo);

  /* Initalise the GSL workspace */
  size_t workspace_size = 100000;
  gsl_integration_workspace *workspace =
      gsl_integration_workspace_alloc(workspace_size);
  gsl_function F = {&drift_integrand, &cosmo};

  /* Loop over all the reasonable time-bins */
  for (int bin = 6; bin < 24; ++bin) {

    const int num_steps = (1LL << bin);
    const integertime_t time_step_size = max_nr_timesteps / num_steps;

    double min_err = 0.;
    double max_err = 0.;
    double sum_err = 0.;
    integertime_t ti_min = 0;
    integertime_t ti_max = 0;

    message("Testing %d steps", num_steps);

    /* Cycle through the time-steps */
    for (integertime_t ti = 0; ti < max_nr_timesteps; ti += time_step_size) {

      const integertime_t ti_beg = ti;
      const integertime_t ti_end =
          min(ti + time_step_size, max_nr_timesteps - 1);

      const double a_beg = cosmology_get_scale_factor(&cosmo, ti_beg);
      const double a_end = cosmology_get_scale_factor(&cosmo, ti_end);

      /* Get the drift factor from SWIFT */
      const double swift_drift_fac =
          cosmology_get_drift_factor(&cosmo, ti, ti + time_step_size);

      /* Get the exact drift factor */
      double exact_drift_fac = 0., abserr;
      gsl_integration_qag(&F, a_beg, a_end, 0, 1.0e-12, workspace_size,
                          GSL_INTEG_GAUSS61, workspace, &exact_drift_fac,
                          &abserr);

      const double rel_err = 0.5 * (swift_drift_fac - exact_drift_fac) /
                             (swift_drift_fac + exact_drift_fac);

      if (rel_err > max_err) {
        max_err = rel_err;
        ti_max = ti;
      }

      if (rel_err < min_err) {
        min_err = rel_err;
        ti_min = ti;
      }

      sum_err += fabs(rel_err);
    }

    message("Max  error: %14e at a=[%.9f %.9f] ", max_err,
            cosmology_get_scale_factor(&cosmo, ti_max),
            cosmology_get_scale_factor(&cosmo, ti_max + time_step_size));
    message("Min  error: %14e at a=[%.9f %.9f]", min_err,
            cosmology_get_scale_factor(&cosmo, ti_min),
            cosmology_get_scale_factor(&cosmo, ti_min + time_step_size));
    message("Sum  error: %14e", sum_err);
    message("Mean error: %14e", sum_err / num_steps);

    if (max_err > 1e-4 || min_err < -1e-4)
      error("Error too large to be acceptable");
  }

#ifdef HAVE_LIBGSL

  return 0;
#else
  return 0;
#endif
}
