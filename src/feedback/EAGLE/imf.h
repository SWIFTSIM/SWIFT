/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_STARS_IMF_H
#define SWIFT_EAGLE_STARS_IMF_H

/* Some standard headers. */
#include <string.h>

/* Local includes. */
#include "exp10.h"
#include "inline.h"
#include "interpolate.h"
#include "minmax.h"
#include "yield_tables.h"

/**
 * @brief the different weightings allowed for the IMF integration
 */
enum eagle_imf_integration_type {
  eagle_imf_integration_no_weight,   /*<! No weighting */
  eagle_imf_integration_mass_weight, /*<! Weighted by mass */
  eagle_imf_integration_yield_weight /*<! Weigthed by stellar yields */
} __attribute__((packed));

/**
 * @brief determine which IMF mass bins the upper and lower input mass bounds
 * belong to
 *
 * @param log10_min_mass Lower mass bound
 * @param log10_max_mass Upper mass bound
 * @param i_min (return) Index of IMF mass bin containing log10_min_mass
 * @param i_max (return) Index of IMF mass bin containing log10_max_mass
 * @param feedback_props the #feedback_props data struct
 */
INLINE static void determine_imf_bins(
    double log10_min_mass, double log10_max_mass, int *i_min, int *i_max,
    const struct feedback_props *feedback_props) {

#ifdef SWIFT_DEBUG_CHECKS
  if (log10_min_mass > log10_max_mass)
    error("Lower bound higher than larger bound.");
#endif

  const int N_bins = eagle_feedback_N_imf_bins;
  const double *const imf_bins_log10 = feedback_props->imf_mass_bin_log10;

  /* Check whether lower mass is within the IMF mass bin range */
  log10_min_mass = max(log10_min_mass, imf_bins_log10[0]);
  log10_min_mass = min(log10_min_mass, imf_bins_log10[N_bins - 1]);

  /* Check whether upper mass is within the IMF mass bin range */
  log10_max_mass = max(log10_max_mass, imf_bins_log10[0]);
  log10_max_mass = min(log10_max_mass, imf_bins_log10[N_bins - 1]);

  *i_min = 0;
  while ((*i_min < N_bins - 2) && imf_bins_log10[*i_min + 1] < log10_min_mass) {
    (*i_min)++;
  }

  *i_max = 1;
  while ((*i_max < N_bins - 1) && imf_bins_log10[*i_max] < log10_max_mass) {
    (*i_max)++;
  }
}

/**
 * @brief Integrate the IMF between a minimum and maximum mass using the
 * trapezoidal rule. The IMF may be weighted by various quantities, as specified
 * by the variable, mode, including an input array, stellar_yields.
 *
 * @param log10_min_mass log10 mass lower integration bound
 * @param log10_max_mass log10 mass upper integration bound
 * @param mode Type of weighting for the IMF integration.
 * @param stellar_yields Array of weights based on yields. Used only for
 * yield-weighted integration.
 * @param feedback_props the #feedback_props data structure
 */
INLINE static double integrate_imf(
    const double log10_min_mass, const double log10_max_mass,
    const enum eagle_imf_integration_type mode,
    const double stellar_yields[eagle_feedback_N_imf_bins],
    const struct feedback_props *feedback_props) {

  /* Pull out some common terms */
  const double *imf = feedback_props->imf;
  const double *imf_mass_bin = feedback_props->imf_mass_bin;
  const double *imf_mass_bin_log10 = feedback_props->imf_mass_bin_log10;

  /* IMF mass bin spacing in log10 space. Assumes uniform spacing. */
  const double imf_log10_mass_bin_size =
      imf_mass_bin_log10[1] - imf_mass_bin_log10[0];

  /* Determine bins to integrate over based on integration bounds */
  int i_min, i_max;
  determine_imf_bins(log10_min_mass, log10_max_mass, &i_min, &i_max,
                     feedback_props);

  /* Array for the integrand */
  double integrand[eagle_feedback_N_imf_bins];

  /* Add up the contribution from each of the IMF mass bins */
  switch (mode) {

    case eagle_imf_integration_no_weight:

      /* Integrate IMF on its own */
      for (int i = i_min; i < i_max + 1; i++) {
        integrand[i] = imf[i] * imf_mass_bin[i];
      }
      break;

    case eagle_imf_integration_mass_weight:

      /* Integrate IMF weighted by mass */
      for (int i = i_min; i < i_max + 1; i++) {
        integrand[i] = imf[i] * imf_mass_bin[i] * imf_mass_bin[i];
      }
      break;

    case eagle_imf_integration_yield_weight:

#ifdef SWIFT_DEBUG_CHECKS
      if (stellar_yields == NULL)
        error(
            "Yield array not passed in despite asking for yield-weighted IMf "
            "integration.");
#endif

      /* Integrate IMF weighted by yields */
      for (int i = i_min; i < i_max + 1; i++) {
        integrand[i] = stellar_yields[i] * imf[i] * imf_mass_bin[i];
      }
      break;

    default:
      error("Invalid mode for IMF integration");
  }

  /* Integrate using trapezoidal rule */
  double result = 0.;
  for (int i = i_min; i < i_max + 1; i++) {
    result += integrand[i];
  }

  /* Update end bins since contribution was overcounted when summing up all
   * entries */
  result -= 0.5 * (integrand[i_min] + integrand[i_max]);

  /* Correct first bin */
  const double first_bin_offset =
      (log10_min_mass - imf_mass_bin_log10[i_min]) / imf_log10_mass_bin_size;

  if (first_bin_offset < 0.5) {
    result -= first_bin_offset * integrand[i_min];
  } else {
    result -= 0.5 * integrand[i_min];
    result -= (first_bin_offset - 0.5) * integrand[i_min + 1];
  }

  /* Correct last bin */
  const double last_bin_offset =
      (log10_max_mass - imf_mass_bin_log10[i_max - 1]) /
      imf_log10_mass_bin_size;

  if (last_bin_offset < 0.5) {
    result -= 0.5 * integrand[i_max];
    result -= (0.5 - last_bin_offset) * integrand[i_max - 1];
  } else {
    result -= (1.0 - last_bin_offset) * integrand[i_max];
  }

  /* The IMF is tabulated in log10, multiply by log10(mass bin size) to get
   * result of integrating IMF */
  return result * imf_log10_mass_bin_size * M_LN10;
}

/* ================================================================
 *  New: IMF models and options-based initialization (YAML-friendly)
 * ================================================================ */

/**
 * @brief Supported IMF model choices when reading from YAML or options.
 */
enum eagle_imf_model {
  eagle_imf_model_chabrier = 0,
  eagle_imf_model_kroupa = 1,
  eagle_imf_model_salpeter = 2,
  eagle_imf_model_custom  = 3
} __attribute__((packed));

/**
 * @brief Options describing an IMF. Any value <= 0 is treated as "unset" and
 * falls back to the model defaults. All units are in Msun and dimensionless
 * slopes (phi ~ m^-alpha).
 */
struct eagle_imf_options {
  enum eagle_imf_model model;      /*<! Named model */

  /* Common optional parameters (interpreted per model) */
  double high_mass_slope;          /*<! alpha_high for m > pivot */
  double low_mass_slope;           /*<! alpha_low  for m <= pivot (Kroupa/custom) */
  double pivot_mass_msun;          /*<! break mass (default 1.0 for Chabrier, 0.5 for Kroupa) */

  /* Chabrier low-mass lognormal parameters (optional) */
  double chabrier_m_c_msun;        /*<! characteristic mass ~ 0.079 */
  double chabrier_sigma_log10;     /*<! dispersion in log10, ~ 0.69 */
};

/**
 * @brief Internal helper: allocate arrays and precompute binning for IMF.
 * Returns the log10 mass bin size via out pointer.
 */
INLINE static void eagle_imf_allocate_arrays(struct feedback_props *feedback_props,
                                            double *imf_log10_mass_bin_size) {
  const double dlog10 =
      (feedback_props->log10_imf_max_mass_msun -
       feedback_props->log10_imf_min_mass_msun) /
      (double)(eagle_feedback_N_imf_bins - 1);

  if (swift_memalign("imf-tables", (void **)&feedback_props->imf,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_N_imf_bins * sizeof(double)) != 0)
    error("Failed to allocate IMF bins table");

  if (swift_memalign("imf-tables", (void **)&feedback_props->imf_mass_bin,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_N_imf_bins * sizeof(double)) != 0)
    error("Failed to allocate IMF bins table");

  if (swift_memalign("imf-tables", (void **)&feedback_props->imf_mass_bin_log10,
                     SWIFT_STRUCT_ALIGNMENT,
                     eagle_feedback_N_imf_bins * sizeof(double)) != 0)
    error("Failed to allocate IMF bins table");

  *imf_log10_mass_bin_size = dlog10;
}

/**
 * @brief Internal helper: normalize IMF so that ∫ m φ(m) dm = 1 across [m_min, m_max].
 */
INLINE static void eagle_imf_normalize(struct feedback_props *feedback_props) {
  const float norm = integrate_imf(
      feedback_props->log10_imf_min_mass_msun,
      feedback_props->log10_imf_max_mass_msun,
      eagle_imf_integration_mass_weight, /* yields */ NULL, feedback_props);
  for (int i = 0; i < eagle_feedback_N_imf_bins; i++) feedback_props->imf[i] /= norm;
}

/**
 * @brief Initialize a Chabrier (2003) IMF with optional overrides.
 * - Low-mass: lognormal with m_c and sigma in log10.
 * - High-mass: power-law with slope alpha_high, matched continuously at pivot.
 */
INLINE static void init_imf_chabrier(struct feedback_props *feedback_props,
                                     double alpha_high_opt,
                                     double pivot_mass_opt,
                                     double m_c_opt,
                                     double sigma_log10_opt) {
  const double alpha_high = (alpha_high_opt > 0.0) ? alpha_high_opt : 2.3; /* default Chabrier */
  const double pivot_mass = (pivot_mass_opt > 0.0) ? pivot_mass_opt : 1.0;
  const double m_c        = (m_c_opt        > 0.0) ? m_c_opt        : 0.079;
  const double sig_log10  = (sigma_log10_opt> 0.0) ? sigma_log10_opt: 0.69;

  double dlog10;
  eagle_imf_allocate_arrays(feedback_props, &dlog10);

  /* Low-mass lognormal normalization constant used historically in SWIFT */
  const double log10_mc = log10(m_c);

  /* Value of the lognormal at the pivot to match continuity */
  const double log10_pivot = log10(pivot_mass);
  const double phi_low_at_pivot =
      0.852464 *
      exp((log10_pivot - log10_mc) * (log10_pivot - log10_mc) /
          (-2.0 * sig_log10 * sig_log10)) /
      pivot_mass;

  /* High-mass normalization to ensure continuity at pivot */
  const double A_high = phi_low_at_pivot * pow(pivot_mass, alpha_high);

  for (int i = 0; i < eagle_feedback_N_imf_bins; i++) {
    const double log10_m = feedback_props->log10_imf_min_mass_msun + i * dlog10;
    const double m = exp10(log10_m);

    feedback_props->imf_mass_bin[i] = m;
    feedback_props->imf_mass_bin_log10[i] = log10_m;

    if (m > pivot_mass) {
      feedback_props->imf[i] = A_high * pow(m, -alpha_high);
    } else {
      feedback_props->imf[i] = 0.852464 *
                               exp((log10_m - log10_mc) * (log10_m - log10_mc) /
                                   (-2.0 * sig_log10 * sig_log10)) /
                               m;
    }
  }

  eagle_imf_normalize(feedback_props);
}

/**
 * @brief Initialize a Kroupa (2001) broken power-law IMF.
 * Defaults: alpha_low=1.3 below pivot=0.5 Msun, alpha_high=2.3 above.
 */
INLINE static void init_imf_kroupa(struct feedback_props *feedback_props,
                                   double alpha_low_opt,
                                   double alpha_high_opt,
                                   double pivot_mass_opt) {
  const double alpha_low  = (alpha_low_opt  > 0.0) ? alpha_low_opt  : 1.3;
  const double alpha_high = (alpha_high_opt > 0.0) ? alpha_high_opt : 2.3;
  const double pivot_mass = (pivot_mass_opt > 0.0) ? pivot_mass_opt : 0.5;

  double dlog10;
  eagle_imf_allocate_arrays(feedback_props, &dlog10);

  /* Choose A_low arbitrarily; continuity determines A_high; mass-normalization comes later. */
  const double A_low = 1.0;
  const double A_high = A_low * pow(pivot_mass, alpha_high - alpha_low);

  for (int i = 0; i < eagle_feedback_N_imf_bins; i++) {
    const double log10_m = feedback_props->log10_imf_min_mass_msun + i * dlog10;
    const double m = exp10(log10_m);

    feedback_props->imf_mass_bin[i] = m;
    feedback_props->imf_mass_bin_log10[i] = log10_m;

    if (m > pivot_mass) {
      feedback_props->imf[i] = A_high * pow(m, -alpha_high);
    } else {
      feedback_props->imf[i] = A_low * pow(m, -alpha_low);
    }
  }

  eagle_imf_normalize(feedback_props);
}

/**
 * @brief Initialize a Salpeter (1955) single power-law IMF.
 * Default slope alpha=2.35.
 */
INLINE static void init_imf_salpeter(struct feedback_props *feedback_props,
                                     double alpha_opt) {
  const double alpha = (alpha_opt > 0.0) ? alpha_opt : 2.35;

  double dlog10;
  eagle_imf_allocate_arrays(feedback_props, &dlog10);

  const double A = 1.0; /* arbitrary; normalize by mass afterwards */

  for (int i = 0; i < eagle_feedback_N_imf_bins; i++) {
    const double log10_m = feedback_props->log10_imf_min_mass_msun + i * dlog10;
    const double m = exp10(log10_m);

    feedback_props->imf_mass_bin[i] = m;
    feedback_props->imf_mass_bin_log10[i] = log10_m;

    feedback_props->imf[i] = A * pow(m, -alpha);
  }

  eagle_imf_normalize(feedback_props);
}

/**
 * @brief Initialize a custom IMF. Behavior:
 * - If both low_mass_slope and high_mass_slope are set (>0) and pivot set (>0):
 *   broken power-law with continuity at pivot.
 * - Else if only high_mass_slope is set (>0): Chabrier low-mass lognormal + that high-mass slope at pivot=1 Msun.
 * - Else: falls back to classic Chabrier.
 */
INLINE static void init_imf_custom(struct feedback_props *feedback_props,
                                   const struct eagle_imf_options *opt) {
  const double alpha_low  = (opt && opt->low_mass_slope  > 0.0) ? opt->low_mass_slope  : -1.0;
  const double alpha_high = (opt && opt->high_mass_slope > 0.0) ? opt->high_mass_slope : -1.0;
  const double pivot_mass = (opt && opt->pivot_mass_msun > 0.0) ? opt->pivot_mass_msun : -1.0;

  if (alpha_low > 0.0 && alpha_high > 0.0 && pivot_mass > 0.0) {
    init_imf_kroupa(feedback_props, alpha_low, alpha_high, pivot_mass);
  } else if (alpha_high > 0.0) {
    /* Chabrier-like with custom high-mass slope */
    init_imf_chabrier(feedback_props, alpha_high,
                      (opt && opt->pivot_mass_msun > 0.0) ? opt->pivot_mass_msun : 1.0,
                      (opt && opt->chabrier_m_c_msun > 0.0) ? opt->chabrier_m_c_msun : 0.079,
                      (opt && opt->chabrier_sigma_log10 > 0.0) ? opt->chabrier_sigma_log10 : 0.69);
  } else {
    init_imf_chabrier(feedback_props, -1.0, -1.0, -1.0, -1.0);
  }
}

/**
 * @brief Initialize IMF from options (e.g. read from YAML). Defaults to
 * Chabrier when options are NULL or fields are unset.
 */
INLINE static void init_imf_from_options(struct feedback_props *feedback_props,
                                         const struct eagle_imf_options *opt) {
  enum eagle_imf_model model = (opt) ? opt->model : eagle_imf_model_chabrier;
  switch (model) {
    case eagle_imf_model_chabrier:
      init_imf_chabrier(feedback_props,
                        (opt ? opt->high_mass_slope     : -1.0),
                        (opt ? opt->pivot_mass_msun     : -1.0),
                        (opt ? opt->chabrier_m_c_msun   : -1.0),
                        (opt ? opt->chabrier_sigma_log10: -1.0));
      break;
    case eagle_imf_model_kroupa:
      init_imf_kroupa(feedback_props,
                      (opt ? opt->low_mass_slope  : -1.0),
                      (opt ? opt->high_mass_slope : -1.0),
                      (opt ? opt->pivot_mass_msun : -1.0));
      break;
    case eagle_imf_model_salpeter:
      init_imf_salpeter(feedback_props, (opt ? opt->high_mass_slope : -1.0));
      break;
    case eagle_imf_model_custom:
    default:
      init_imf_custom(feedback_props, opt);
      break;
  }
}

/**
 * @brief Backward-compatible wrapper: initialize IMF with default (Chabrier) high-mass slope.
 */
INLINE static void init_imf(struct feedback_props *feedback_props) {
  init_imf_chabrier(feedback_props, /*alpha_high=*/-1.0, /*pivot=*/-1.0,
                    /*m_c=*/-1.0, /*sigma_log10=*/-1.0);
}

/**
 * @brief Calculate mass (in solar masses) of stars that died from the star
 * particle's birth up to its current age (in Gyr).
 *
 * Calculation uses the tables of Portinari et al. 1998, A&A, 334, 505
 *
 * @param age_Gyr age of star in Gyr.
 * @param Z Star's metallicity (metal mass fraction).
 * @param feedback_props the #feedback_props data structure.
 * @return Mass of stars died up to that age in solar masses.
 */
INLINE static double dying_mass_msun(
    const double age_Gyr, const double Z,
    const struct feedback_props *feedback_props) {

  /* Pull out some common terms */
  const double *lifetime_Z = feedback_props->lifetimes.metallicity;
  const double *lifetime_m = feedback_props->lifetimes.mass;
  double **const dying_times = feedback_props->lifetimes.dyingtime;
  const int n_Z = eagle_feedback_lifetime_N_metals;
  const int n_m = eagle_feedback_lifetime_N_masses;

  /* Early abort? */
  if (age_Gyr <= 0.) {
    return feedback_props->imf_max_mass_msun;
  }

  const double log10_age_yr = log10(age_Gyr * 1.0e9);

  /* Calculate index along the metallicity axis */
  int Z_index;
  double Z_offset;
  if (Z <= lifetime_Z[0]) {

    /* Before start of the table */
    Z_index = 0;
    Z_offset = 0.;

  } else if (Z >= lifetime_Z[n_Z - 1]) {

    /* After end of the table */
    Z_index = n_Z - 2;
    Z_offset = 1.;

  } else {

    /* Normal case: Somewhere inside the table */
    Z_index = 0;
    while (Z_index < n_Z - 1 && lifetime_Z[Z_index + 1] <= Z) {
      Z_index++;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (Z_index >= n_Z) error("Z_index is beyond the range  of the table");
#endif

    Z_offset = (Z - lifetime_Z[Z_index]) /
               (lifetime_Z[Z_index + 1] - lifetime_Z[Z_index]);
  }

  /* Check whether we are not beyond the age table for the low metallicity end
   */
  int time_index_lowZ = -1;
  double time_offset_lowZ = 0.;
  if (log10_age_yr >= dying_times[Z_index][0]) {

    /* Before start of the table */
    time_index_lowZ = 0;
    time_offset_lowZ = 0.;

  } else if (log10_age_yr <= dying_times[Z_index][n_m - 1]) {

    /* After end of the table */
    time_index_lowZ = n_m - 2;
    time_offset_lowZ = 1.;
  }

  /* Check whether we are not beyond the age table for the high metallicity end
   */
  int time_index_highZ = -1;
  double time_offset_highZ = 0.;
  if (log10_age_yr >= dying_times[Z_index + 1][0]) {

    /* Before start of the table */
    time_index_highZ = 0;
    time_offset_highZ = 0.;

  } else if (log10_age_yr <= dying_times[Z_index + 1][n_m - 1]) {

    /* After end of the table */
    time_index_highZ = n_m - 2;
    time_offset_highZ = 1.0;
  }

  /* Search the table starting from the largest times until we reach
     a solution for the low-metallicity bound */
  int i = n_m - 1;
  while (i >= 0 && time_index_lowZ == -1) {

    if (dying_times[Z_index][i] >= log10_age_yr && time_index_lowZ == -1) {

      /* record index */
      time_index_lowZ = i;

      /* record distance from table element */
      time_offset_lowZ =
          (log10_age_yr - dying_times[Z_index][time_index_lowZ]) /
          (dying_times[Z_index][time_index_lowZ + 1] -
           dying_times[Z_index][time_index_lowZ]);

      break;
    }
    i--;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (time_index_lowZ == -1) error("Could not find low-metallicity bound!");
#endif

  /* Search the table starting from the largest times until we reach
     a solution for the high-metallicity bound */
  i = n_m - 1;
  while (i >= 0 && time_index_highZ == -1) {

    if (dying_times[Z_index + 1][i] >= log10_age_yr && time_index_highZ == -1) {

      /* record index */
      time_index_highZ = i;

      /* record distance from table element */
      time_offset_highZ =
          (log10_age_yr - dying_times[Z_index + 1][time_index_highZ]) /
          (dying_times[Z_index + 1][time_index_highZ + 1] -
           dying_times[Z_index + 1][time_index_highZ]);

      break;
    }
    i--;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (time_index_highZ == -1) error("Could not find high-metallicity bound!");
#endif

  /* And now interpolate the solution */
  const double mass_low_Z =
      interpolate_1d(lifetime_m, time_index_lowZ, time_offset_lowZ);
  const double mass_high_Z =
      interpolate_1d(lifetime_m, time_index_highZ, time_offset_highZ);

  double mass = (1. - Z_offset) * mass_low_Z + Z_offset * mass_high_Z;

  /* Check that we haven't killed too many stars */
  mass = min(mass, feedback_props->imf_max_mass_msun);

  return mass;
}

/**
 * @brief Calculate lifetime of stellar population in Gyr for a given mass.
 *
 * Calculation uses the tables of Portinari et al. 1998, A&A, 334, 505
 *
 * @param mass in solar masses.
 * @param Z Metallicity (metal mass fraction).
 * @param feedback_props the #feedback_props data structure.
 * @return The life time in Giga-years.
 */
INLINE static float lifetime_in_Gyr(
    const float mass, const float Z,
    const struct feedback_props *feedback_props) {

  /* Pull out some common terms */
  const double *lifetime_Z = feedback_props->lifetimes.metallicity;
  const double *lifetime_m = feedback_props->lifetimes.mass;
  double **const dying_times = feedback_props->lifetimes.dyingtime;
  const int n_Z = eagle_feedback_lifetime_N_metals;
  const int n_m = eagle_feedback_lifetime_N_masses;

  /* Calculate index along the mass axis */
  int m_index;
  float m_offset;
  if (mass <= lifetime_m[0]) {

    /* Before start of the table */
    m_index = 0;
    m_offset = 0.f;

  } else if (mass >= lifetime_m[n_m - 1]) {

    /* After end of the table */
    m_index = n_m - 2;
    m_offset = 1.f;

  } else {

    /* Normal case: Somewhere inside the table */
    for (m_index = 0; m_index < n_m - 1; m_index++)
      if (lifetime_m[m_index + 1] > mass) break;

    m_offset = (mass - lifetime_m[m_index]) /
               (lifetime_m[m_index + 1] - lifetime_m[m_index]);
  }

  /* Calculate index along the metallicity axis */
  int Z_index;
  float Z_offset;
  if (Z <= lifetime_Z[0]) {

    /* Before start of the table */
    Z_index = 0;
    Z_offset = 0.f;

  } else if (Z >= lifetime_Z[n_Z - 1]) {

    /* After end of the table */
    Z_index = n_Z - 2;
    Z_offset = 1.f;

  } else {

    for (Z_index = 0; Z_index < n_Z - 1; Z_index++)
      if (lifetime_Z[Z_index + 1] > Z) break;

    /* Normal case: Somewhere inside the table */
    Z_offset = (Z - lifetime_Z[Z_index]) /
               (lifetime_Z[Z_index + 1] - lifetime_Z[Z_index]);
  }

  /* Interpolation of the table to get the time */
  const float log_time_years =
      interpolate_2d(dying_times, Z_index, m_index, Z_offset, m_offset);

  /* Convert to Giga-years */
  const float time_Gyr = exp10f(log_time_years - 9.f);

  return time_Gyr;
}

#endif
