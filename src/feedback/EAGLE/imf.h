/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#include "interpolate.h"
#include "minmax.h"

#include <string.h>

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
 * @param star_properties the #stars_props data struct
 */
inline static void determine_imf_bins(
    double log10_min_mass, double log10_max_mass, int *i_min, int *i_max,
    const struct feedback_props *feedback_props) {

#ifdef SWIFT_DEBUG_CHECKS
  if (log10_min_mass > log10_max_mass)
    error("Lower bound higher than larger bound.");
#endif

  const int N_bins = feedback_props->n_imf_mass_bins;
  const float *imf_bins_log10 = feedback_props->imf_mass_bin_log10;

  /* Check whether lower mass is within the IMF mass bin range */
  log10_min_mass = max(log10_min_mass, imf_bins_log10[0]);
  log10_min_mass = min(log10_min_mass, imf_bins_log10[N_bins - 1]);

  /* Check whether upper mass is within the IMF mass bin range */
  log10_max_mass = max(log10_min_mass, imf_bins_log10[0]);
  log10_max_mass = min(log10_min_mass, imf_bins_log10[N_bins - 1]);

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
 * @param star_properties the #stars_props data structure
 */
inline static float integrate_imf(const float log10_min_mass,
                                  const float log10_max_mass,
                                  const enum eagle_imf_integration_type mode,
                                  const float *stellar_yields,
                                  const struct feedback_props *feedback_props) {

  /* Pull out some common terms */
  const int N_bins = feedback_props->n_imf_mass_bins;
  const float *imf = feedback_props->imf;
  const float *imf_mass_bin = feedback_props->imf_mass_bin;
  const float *imf_mass_bin_log10 =
      feedback_props->imf_mass_bin_log10;

  /* IMF mass bin spacing in log10 space. Assumes uniform spacing. */
  const float imf_log10_mass_bin_size =
      imf_mass_bin_log10[1] - imf_mass_bin_log10[0];

  /* Determine bins to integrate over based on integration bounds */
  int i_min, i_max;
  determine_imf_bins(log10_min_mass, log10_max_mass, &i_min, &i_max,
                     feedback_props);

  /* Array for the integrand */
  float integrand[N_bins];

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
  float result = 0.f;
  for (int i = i_min; i < i_max + 1; i++) {
    result += integrand[i];
  }

  /* Update end bins since contribution was overcounted when summing up all
   * entries */
  result -= 0.5 * (integrand[i_min] + integrand[i_max]);

  /* Correct first bin */
  const float first_bin_offset =
      (log10_min_mass - imf_mass_bin_log10[i_min]) / imf_log10_mass_bin_size;

  if (first_bin_offset < 0.5f) {
    result -= first_bin_offset * integrand[i_min];
  } else {
    result -= 0.5f * integrand[i_min];
    result -= (first_bin_offset - 0.5f) * integrand[i_min + 1];
  }

  /* Correct last bin */
  const float last_bin_offset =
      (log10_max_mass - imf_mass_bin_log10[i_max - 1]) /
      imf_log10_mass_bin_size;

  if (last_bin_offset < 0.5) {
    result -= 0.5f * integrand[i_max];
    result -= (0.5f - last_bin_offset) * integrand[i_max - 1];
  } else {
    result -= (1.f - last_bin_offset) * integrand[i_max];
  }

  /* The IMF is tabulated in log10, multiply by log10(mass bin size) to get
   * result of integrating IMF */
  return result * imf_log10_mass_bin_size * ((float)M_LN10);
}

/**
 * @brief Allocate space for IMF table and compute values to populate this
 * table.
 *
 * @param star_properties #stars_props data structure */
inline static void init_imf(struct feedback_props *restrict feedback_props) {

  /* Define number of imf mass bins */
  feedback_props->n_imf_mass_bins = 200;

  /* Define max and min imf masses */
  feedback_props->imf_max_mass_msun = 100.f;
  feedback_props->imf_min_mass_msun = 0.1;
  feedback_props->log10_imf_max_mass_msun =
      log10(feedback_props->imf_max_mass_msun);
  feedback_props->log10_imf_min_mass_msun =
      log10(feedback_props->imf_min_mass_msun);

  /* Compute size of mass bins in log10 space */
  const float imf_log10_mass_bin_size =
      (feedback_props->log10_imf_max_mass_msun -
       feedback_props->log10_imf_min_mass_msun) /
      (float)(feedback_props->n_imf_mass_bins - 1);

  /* Allocate IMF array */
  if (swift_memalign(
          "imf-tables", (void **)&feedback_props->imf,
          SWIFT_STRUCT_ALIGNMENT,
          feedback_props->n_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");

  /* Allocate array to store IMF mass bins */
  if (swift_memalign(
          "imf-tables", (void **)&feedback_props->imf_mass_bin,
          SWIFT_STRUCT_ALIGNMENT,
          feedback_props->n_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");

  /* Allocate array to store IMF mass bins in log10 space */
  if (swift_memalign(
          "imf-tables", (void **)&feedback_props->imf_mass_bin_log10,
          SWIFT_STRUCT_ALIGNMENT,
          feedback_props->n_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");

  /* Set IMF from Chabrier 2003 */
  if (strcmp(feedback_props->IMF_Model, "Chabrier") == 0) {
    for (int i = 0; i < feedback_props->n_imf_mass_bins; i++) {

      const float log10_mass_msun =
          feedback_props->log10_imf_min_mass_msun +
          i * imf_log10_mass_bin_size;

      const float mass_msun = exp10f(log10_mass_msun);

      if (mass_msun > 1.0)
        feedback_props->imf[i] = 0.237912 * pow(mass_msun, -2.3);
      else
        feedback_props->imf[i] =
            0.852464 *
            exp((log10(mass_msun) - log10(0.079)) *
                (log10(mass_msun) - log10(0.079)) / (-2.0 * pow(0.69, 2))) /
            mass_msun;

      feedback_props->imf_mass_bin[i] = mass_msun;
      feedback_props->imf_mass_bin_log10[i] = log10_mass_msun;
    }
  } else {
    error("Invalid IMF model %s. Valid models are: 'Chabrier'\n",
          feedback_props->IMF_Model);
  }

  /* Normalize the IMF */
  const float norm =
      integrate_imf(feedback_props->log10_imf_min_mass_msun,
                    feedback_props->log10_imf_max_mass_msun,
                    eagle_imf_integration_mass_weight,
                    /*(stellar_yields=)*/ NULL, feedback_props);

  for (int i = 0; i < feedback_props->n_imf_mass_bins; i++)
    feedback_props->imf[i] /= norm;
}

/**
 * @brief Calculate mass (in solar masses) of stars that died from the star
 * particle's birth up to its current age (in Gyr).
 *
 * Calculation follows Portinari et al. 1998
 *
 * @param age_Gyr age of star in Gyr
 * @param metallicity Star's metallicity
 * @param star_properties the #stars_props data structure
 */
inline static double dying_mass_msun(const float age_Gyr,
                                     const float metallicity,
                                     const struct feedback_props *star_props) {

  float metallicity_offset, time_offset_low_metallicity = 0,
                            time_offset_high_metallicity = 0,
                            mass_low_metallicity, mass_high_metallicity;

  int metallicity_index, time_index_low_metallicity = -1,
                         time_index_high_metallicity = -1;

  if (age_Gyr <= 0) {
    return star_props->imf_max_mass_msun;
  }

  const float log10_age_yr = log10f(age_Gyr * 1.e9);

  if (metallicity <= star_props->lifetimes.metallicity[0]) {
    metallicity_index = 0;
    metallicity_offset = 0.0;
  } else if (metallicity >=
             star_props->lifetimes
                 .metallicity[star_props->lifetimes.n_z - 1]) {
    metallicity_index = star_props->lifetimes.n_z - 2;
    metallicity_offset = 1.0;
  } else {
    for (metallicity_index = 0;
         metallicity_index < star_props->lifetimes.n_z - 1;
         metallicity_index++)
      if (star_props->lifetimes.metallicity[metallicity_index + 1] >
          metallicity)
        break;

    metallicity_offset =
        (metallicity -
         star_props->lifetimes.metallicity[metallicity_index]) /
        (star_props->lifetimes.metallicity[metallicity_index + 1] -
         star_props->lifetimes.metallicity[metallicity_index]);
  }

  if (log10_age_yr >=
      star_props->lifetimes.dyingtime[metallicity_index][0]) {
    time_index_low_metallicity = 0;
    time_offset_low_metallicity = 0.0;
  } else if (log10_age_yr <=
             star_props->lifetimes
                 .dyingtime[metallicity_index]
                           [star_props->lifetimes.n_mass - 1]) {
    time_index_low_metallicity = star_props->lifetimes.n_mass - 2;
    time_offset_low_metallicity = 1.0;
  }

  if (log10_age_yr >=
      star_props->lifetimes.dyingtime[metallicity_index + 1][0]) {
    time_index_high_metallicity = 0;
    time_offset_high_metallicity = 0.0;
  } else if (log10_age_yr <=
             star_props->lifetimes
                 .dyingtime[metallicity_index + 1]
                           [star_props->lifetimes.n_mass - 1]) {
    time_index_high_metallicity = star_props->lifetimes.n_mass - 2;
    time_offset_high_metallicity = 1.0;
  }

  int i = star_props->lifetimes.n_mass;
  while (i >= 0 && (time_index_low_metallicity == -1 ||
                    time_index_high_metallicity == -1)) {
    i--;
    if (star_props->lifetimes.dyingtime[metallicity_index][i] >=
            log10_age_yr &&
        time_index_low_metallicity == -1) {
      time_index_low_metallicity = i;
      time_offset_low_metallicity =
          (log10_age_yr -
           star_props->lifetimes
               .dyingtime[metallicity_index][time_index_low_metallicity]) /
          (star_props->lifetimes
               .dyingtime[metallicity_index][time_index_low_metallicity + 1] -
           star_props->lifetimes
               .dyingtime[metallicity_index][time_index_low_metallicity]);
    }
    if (star_props->lifetimes.dyingtime[metallicity_index + 1][i] >=
            log10_age_yr &&
        time_index_high_metallicity == -1) {
      time_index_high_metallicity = i;
      time_offset_high_metallicity =
          (log10_age_yr -
           star_props->lifetimes
               .dyingtime[metallicity_index + 1][time_index_high_metallicity]) /
          (star_props->lifetimes
               .dyingtime[metallicity_index + 1]
                         [time_index_high_metallicity + 1] -
           star_props->lifetimes
               .dyingtime[metallicity_index + 1][time_index_high_metallicity]);
    }
  }

  mass_low_metallicity =
      interpolate_1d(star_props->lifetimes.mass,
                     time_index_low_metallicity, time_offset_low_metallicity);
  mass_high_metallicity =
      interpolate_1d(star_props->lifetimes.mass,
                     time_index_high_metallicity, time_offset_high_metallicity);

  float mass = (1.0 - metallicity_offset) * mass_low_metallicity +
               metallicity_offset * mass_high_metallicity;

  /* Check that we haven't killed too many stars */
  if (mass > star_props->imf_max_mass_msun)
    mass = star_props->imf_max_mass_msun;

  return mass;
}

/**
 * @brief Calculate lifetime of star poputlation in Gyr. Approach based on
 * Portinari et al. 1998
 *
 * @param mass
 * @param metallicity
 * @param star_properties the #stars_props data struct
 */
inline static float lifetime_in_Gyr(float mass, float metallicity,
                                    const struct feedback_props *restrict
                                        feedback_props) {

  double time_Gyr = 0, mass_offset, metallicity_offset;

  int mass_index, metallicity_index;

  if (mass <= feedback_props->lifetimes.mass[0]) {
    mass_index = 0;
    mass_offset = 0.0;
  } else if (mass >=
             feedback_props->lifetimes
                 .mass[feedback_props->lifetimes.n_mass - 1]) {
    mass_index = feedback_props->lifetimes.n_mass - 2;
    mass_offset = 1.0;
  } else {
    for (mass_index = 0;
         mass_index < feedback_props->lifetimes.n_mass - 1;
         mass_index++)
      if (feedback_props->lifetimes.mass[mass_index + 1] > mass)
        break;

    mass_offset =
        (mass - feedback_props->lifetimes.mass[mass_index]) /
        (feedback_props->lifetimes.mass[mass_index + 1] -
         feedback_props->lifetimes.mass[mass_index]);
  }

  if (metallicity <= feedback_props->lifetimes.metallicity[0]) {
    metallicity_index = 0;
    metallicity_offset = 0.0;
  } else if (metallicity >=
             feedback_props->lifetimes
                 .metallicity[feedback_props->lifetimes.n_z - 1]) {
    metallicity_index = feedback_props->lifetimes.n_z - 2;
    metallicity_offset = 1.0;
  } else {
    for (metallicity_index = 0;
         metallicity_index < feedback_props->lifetimes.n_z - 1;
         metallicity_index++)
      if (feedback_props->lifetimes
              .metallicity[metallicity_index + 1] > metallicity)
        break;

    metallicity_offset =
        (metallicity -
         feedback_props->lifetimes.metallicity[metallicity_index]) /
        (feedback_props->lifetimes
             .metallicity[metallicity_index + 1] -
         feedback_props->lifetimes.metallicity[metallicity_index]);
  }

  /* time in Gyr */
  time_Gyr =
      exp(M_LN10 * interpolate_2d(feedback_props->lifetimes.dyingtime,
                                  metallicity_index, mass_index,
                                  metallicity_offset, mass_offset)) /
      1.0e9;

  return time_Gyr;
}

#endif
