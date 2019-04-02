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

/**
 * @brief determine which IMF mass bins the upper and lower input mass bounds belong to
 *
 * @param log10_min_mass Lower mass bound
 * @param log10_max_mass Upper mass bound
 * @param low_imf_mass_bin_index pointer to index of IMF mass bin containing log10_min_mass
 * @param high_imf_mass_bin_index pointer to index of IMF mass bin containing log10_max_mass
 * @param star_properties the #stars_props data struct */
inline static void determine_imf_bins(
    double log10_min_mass, double log10_max_mass, int *low_imf_mass_bin_index, int *high_imf_mass_bin_index,
    const struct stars_props *restrict star_properties) {
  int i1, i2;

  /* Check whether lower mass is within the IMF mass bin range */
  if (log10_min_mass < star_properties->feedback.imf_mass_bin_log10[0])
    log10_min_mass = star_properties->feedback.imf_mass_bin_log10[0];
  if (log10_min_mass >
      star_properties->feedback.imf_mass_bin_log10[star_properties->feedback.n_imf_mass_bins - 1])
    log10_min_mass =
        star_properties->feedback.imf_mass_bin_log10[star_properties->feedback.n_imf_mass_bins - 1];

  /* Check whether upper mass is within the IMF mass bin range */
  if (log10_max_mass < star_properties->feedback.imf_mass_bin_log10[0])
    log10_max_mass = star_properties->feedback.imf_mass_bin_log10[0];
  if (log10_max_mass >
      star_properties->feedback.imf_mass_bin_log10[star_properties->feedback.n_imf_mass_bins - 1])
    log10_max_mass =
        star_properties->feedback.imf_mass_bin_log10[star_properties->feedback.n_imf_mass_bins - 1];

  /* Find mass bin to which lower mass belongs */
  for (i1 = 0; i1 < star_properties->feedback.n_imf_mass_bins - 2 &&
               star_properties->feedback.imf_mass_bin_log10[i1 + 1] < log10_min_mass;
       i1++);
    
  /* Find mass bin to which upper mass belongs */
  for (i2 = 1; i2 < star_properties->feedback.n_imf_mass_bins - 1 &&
               star_properties->feedback.imf_mass_bin_log10[i2] < log10_max_mass;
       i2++);
    
  *low_imf_mass_bin_index = i1;
  *high_imf_mass_bin_index = i2;
}

/**
 * @brief Integrate the IMF between a minimum and maximum mass using the trapezoidal rule. The IMF may be weighted by various quantities, as specified by the variable, mode, including an input array, stellar_yields.
 *
 * @param log10_min_mass log10 mass lower integration bound
 * @param log10_max_mass log10 mass upper integration bound
 * @param mode Flag to specify weighting used for integrating IMF
 * @param stellar_yields Array of weights based on yields. Used with mode=2
 * @param star_properties the #stars_props data structure
 */
inline static float integrate_imf(
    float log10_min_mass, float log10_max_mass, int mode,
    float *stellar_yields, const struct stars_props *restrict star_properties) {

  double result;

  int low_imf_mass_bin_index, high_imf_mass_bin_index;

  float imf_log10_mass_bin_size, bin_offset;

  imf_log10_mass_bin_size = star_properties->feedback.imf_mass_bin_log10[1] -
        star_properties->feedback.imf_mass_bin_log10[0]; /* dlog(m) */

  /* Determine bins to integrate over based on integration bounds */
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index, &high_imf_mass_bin_index,
                     star_properties);

  // ALEXEI: rewrite this function so that adding on the fly and don't need to use the array.
  float integrand[star_properties->feedback.n_imf_mass_bins];

  /* Add up the contribution from each of the IMF mass bins */
  if (mode == 0)
    /* Integrate IMF on its own */
    for (int i = low_imf_mass_bin_index; i < high_imf_mass_bin_index + 1; i++)
      integrand[i] =
          star_properties->feedback.imf[i] *
          star_properties->feedback.imf_mass_bin[i]; 
  else if (mode == 1)
    /* Integrate IMF weighted by mass */
    for (int i = low_imf_mass_bin_index; i < high_imf_mass_bin_index + 1; i++)
      integrand[i] =
          star_properties->feedback.imf[i] *
          star_properties->feedback.imf_mass_bin[i] *
          star_properties->feedback.imf_mass_bin[i]; 
  else if (mode == 2)
    /* Integrate IMF weighted by yields */
    for (int i = low_imf_mass_bin_index; i < high_imf_mass_bin_index + 1; i++)
      integrand[i] =
          stellar_yields[i] * star_properties->feedback.imf[i] *
          star_properties->feedback.imf_mass_bin[i]; 
  else {
    error("invalid mode in integrate_imf = %d\n", mode);
  }

  /* integrate using trapezoidal rule (ALEXEI: Remove when adding on the fly)*/
  result = 0;
  for (int i = low_imf_mass_bin_index; i < high_imf_mass_bin_index + 1; i++) result += integrand[i];

  /* Update end bins since contribution was overcounted when summing up all entries */
  result = result - 0.5 * integrand[low_imf_mass_bin_index] - 0.5 * integrand[high_imf_mass_bin_index];

  /* correct first bin (ALEXEI: why are we doing this? (just to clarify comment))*/
  bin_offset = (log10_min_mass - star_properties->feedback.imf_mass_bin_log10[low_imf_mass_bin_index]) / imf_log10_mass_bin_size;

  if (bin_offset < 0.5)
    result -= bin_offset * integrand[low_imf_mass_bin_index];
  else {
    result -= 0.5 * integrand[low_imf_mass_bin_index];
    result -= (bin_offset - 0.5) * integrand[low_imf_mass_bin_index + 1];
  }

  /* correct last bin */
  bin_offset = (log10_max_mass - star_properties->feedback.imf_mass_bin_log10[high_imf_mass_bin_index - 1]) / imf_log10_mass_bin_size;

  if (bin_offset < 0.5) {
    result -= 0.5 * integrand[high_imf_mass_bin_index];
    result -= (0.5 - bin_offset) * integrand[high_imf_mass_bin_index - 1];
  } else {
    result -= (1 - bin_offset) * integrand[high_imf_mass_bin_index];
  }

  /* The IMF is tabulated in log10 so undo to get regular units */
  result *= imf_log10_mass_bin_size * log(10.0); 

  return result;
}

/**
 * @brief Allocate space for IMF table and compute values to populate this table.
 *
 * @param star_properties #stars_props data structure */
inline static void init_imf(struct stars_props *restrict star_properties) {

  float norm = 0, mass_msun, log10_mass_msun;

  /* Define number of imf mass bins */
  star_properties->feedback.n_imf_mass_bins = 200;

  /* Define max and min imf masses */
  star_properties->feedback.imf_max_mass_msun = 100.f;
  star_properties->feedback.imf_min_mass_msun = 0.1;
  star_properties->feedback.log10_imf_min_mass_msun = -1;
  star_properties->feedback.log10_imf_max_mass_msun = 2;

  /* Compute size of mass bins in log10 space */
  const float imf_log10_mass_bin_size = (star_properties->feedback.log10_imf_max_mass_msun - star_properties->feedback.log10_imf_min_mass_msun) /
                    (double)(star_properties->feedback.n_imf_mass_bins - 1);

  /* Allocate IMF array */
  if (posix_memalign((void **)&star_properties->feedback.imf,
                     SWIFT_STRUCT_ALIGNMENT,
                     star_properties->feedback.n_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");

  /* Allocate array to store IMF mass bins */
  if (posix_memalign((void **)&star_properties->feedback.imf_mass_bin,
                     SWIFT_STRUCT_ALIGNMENT,
                     star_properties->feedback.n_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");

  /* Allocate array to store IMF mass bins in log10 space */
  if (posix_memalign((void **)&star_properties->feedback.imf_mass_bin_log10,
                     SWIFT_STRUCT_ALIGNMENT,
                     star_properties->feedback.n_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");

  /* Set Power-law IMF (Salpeter for IMF_EXPONENT = 2.35) */
  if (strcmp(star_properties->feedback.IMF_Model, "PowerLaw") == 0) {
    if (star_properties->feedback.IMF_Exponent < 0) {
      error("imf_exponent is supposed to be > 0\n");
    }

    for (int i = 0; i < star_properties->feedback.n_imf_mass_bins; i++) {
      log10_mass_msun = star_properties->feedback.log10_imf_min_mass_msun + i * imf_log10_mass_bin_size;
      mass_msun = exp(M_LN10 * log10_mass_msun);

      star_properties->feedback.IMF_Exponent = 2.35;
      star_properties->feedback.imf[i] =
          pow(mass_msun, -star_properties->feedback.IMF_Exponent);

      star_properties->feedback.imf_mass_bin[i] = mass_msun;
      star_properties->feedback.imf_mass_bin_log10[i] = log10_mass_msun;
    }
  }
  /* Set IMF from Chabrier 2003 */
  else if (strcmp(star_properties->feedback.IMF_Model, "Chabrier") == 0) {
    for (int i = 0; i < star_properties->feedback.n_imf_mass_bins; i++) {
      log10_mass_msun = star_properties->feedback.log10_imf_min_mass_msun + i * imf_log10_mass_bin_size;
      mass_msun = exp(M_LN10 * log10_mass_msun);

      if (mass_msun > 1.0)
        star_properties->feedback.imf[i] = 0.237912 * pow(mass_msun, -2.3);
      else
        star_properties->feedback.imf[i] =
            0.852464 *
            exp((log10(mass_msun) - log10(0.079)) *
                (log10(mass_msun) - log10(0.079)) / (-2.0 * pow(0.69, 2))) /
            mass_msun;

      star_properties->feedback.imf_mass_bin[i] = mass_msun;
      star_properties->feedback.imf_mass_bin_log10[i] = log10_mass_msun;
    }
  } else {
    error("Invalid IMF model %s. Valid models are: PowerLaw and Chabrier\n",
          star_properties->feedback.IMF_Model);
  }

  /* Normalize the IMF */
  norm = integrate_imf(star_properties->feedback.log10_imf_min_mass_msun, star_properties->feedback.log10_imf_max_mass_msun, 1,
                       /*(stellar_yields=)*/ NULL, star_properties);

  for (int i = 0; i < star_properties->feedback.n_imf_mass_bins; i++)
    star_properties->feedback.imf[i] /= norm;
}

/**
 * @brief Calculate mass (in solar masses) of stars that died from the star particle's birth up to its current age (in Gyr)
 *
 * @param age_Gyr age of star in Gyr
 * @param metallicity Star's metallicity
 * @param star_properties the #stars_props data structure
 */
inline static double dying_mass_msun(double age_Gyr, float metallicity,
                                    const struct stars_props* restrict
                                        star_properties) {

  float mass = 0, metallicity_offset, time_offset_low_metallicity = 0, time_offset_high_metallicity = 0, log10_age_yr, mass_low_metallicity, mass_high_metallicity;

  int metallicity_index, i, time_index_low_metallicity = -1, time_index_high_metallicity = -1;

  /* Calculate star mass that dies based on one of several analytical models */
  switch (star_properties->feedback.stellar_lifetime_flag) {
    case 0:
      /* padovani & matteucci 1993 */

      if (age_Gyr > 0.039765318659064693)
        mass =
            exp(M_LN10 *
                (7.764 - (1.79 - (1.338 - 0.1116 * (9 + log10(age_Gyr))) *
                                     (1.338 - 0.1116 * (9 + log10(age_Gyr)))) /
                             0.2232));
      else if (age_Gyr > 0.003)
        mass = pow((age_Gyr - 0.003) / 1.2, -1 / 1.85);
      else
        mass = star_properties->feedback.imf_max_mass_msun;
      break;

    case 1:
      /* maeder & meynet 1989 */

      if (age_Gyr >= 8.4097378)
        mass = exp(M_LN10 * (1 - log10(age_Gyr)) / 0.6545);
      else if (age_Gyr >= 0.35207776)
        mass = exp(M_LN10 * (1.35 - log10(age_Gyr)) / 3.7);
      else if (age_Gyr >= 0.050931493)
        mass = exp(M_LN10 * (0.77 - log10(age_Gyr)) / 2.51);
      else if (age_Gyr >= 0.010529099)
        mass = exp(M_LN10 * (0.17 - log10(age_Gyr)) / 1.78);
      else if (age_Gyr >= 0.0037734787)
        mass = exp(M_LN10 * (-0.94 - log10(age_Gyr)) / 0.86);
      else if (age_Gyr > 0.003)
        mass = pow((age_Gyr - 0.003) / 1.2, -0.54054053);
      else
        mass = star_properties->feedback.imf_max_mass_msun;
      break;

    case 2:
      /* portinari et al. 1998 */

      if (age_Gyr <= 0) {
        mass = star_properties->feedback.imf_max_mass_msun;
        break;
      }

      log10_age_yr = log10(age_Gyr * 1.e9);

      if (metallicity <= star_properties->feedback.lifetimes.metallicity[0]) {
        metallicity_index = 0;
        metallicity_offset = 0.0;
      } else if (metallicity >=
                 star_properties->feedback.lifetimes
                     .metallicity[star_properties->feedback.lifetimes.n_z - 1]) {
        metallicity_index = star_properties->feedback.lifetimes.n_z - 2;
        metallicity_offset = 1.0;
      } else {
        for (metallicity_index = 0; metallicity_index < star_properties->feedback.lifetimes.n_z - 1;
             metallicity_index++)
          if (star_properties->feedback.lifetimes.metallicity[metallicity_index + 1] >
              metallicity)
            break;

        metallicity_offset = (metallicity -
                   star_properties->feedback.lifetimes.metallicity[metallicity_index]) /
                  (star_properties->feedback.lifetimes.metallicity[metallicity_index + 1] -
                   star_properties->feedback.lifetimes.metallicity[metallicity_index]);
      }

      if (log10_age_yr >= star_properties->feedback.lifetimes.dyingtime[metallicity_index][0]) {
        time_index_low_metallicity = 0;
        time_offset_low_metallicity = 0.0;
      } else if (log10_age_yr <=
                 star_properties->feedback.lifetimes
                     .dyingtime[metallicity_index]
                               [star_properties->feedback.lifetimes.n_mass - 1]) {
        time_index_low_metallicity = star_properties->feedback.lifetimes.n_mass - 2;
        time_offset_low_metallicity = 1.0;
      }

      if (log10_age_yr >=
          star_properties->feedback.lifetimes.dyingtime[metallicity_index + 1][0]) {
        time_index_high_metallicity = 0;
        time_offset_high_metallicity = 0.0;
      } else if (log10_age_yr <=
                 star_properties->feedback.lifetimes
                     .dyingtime[metallicity_index + 1]
                               [star_properties->feedback.lifetimes.n_mass - 1]) {
        time_index_high_metallicity = star_properties->feedback.lifetimes.n_mass - 2;
        time_offset_high_metallicity = 1.0;
      }

      i = star_properties->feedback.lifetimes.n_mass;
      while (i >= 0 && (time_index_low_metallicity == -1 || time_index_high_metallicity == -1)) {
        i--;
        if (star_properties->feedback.lifetimes.dyingtime[metallicity_index][i] >=
                log10_age_yr &&
            time_index_low_metallicity == -1) {
          time_index_low_metallicity = i;
          time_offset_low_metallicity =
              (log10_age_yr -
               star_properties->feedback.lifetimes.dyingtime[metallicity_index][time_index_low_metallicity]) /
              (star_properties->feedback.lifetimes
                   .dyingtime[metallicity_index][time_index_low_metallicity + 1] -
               star_properties->feedback.lifetimes.dyingtime[metallicity_index][time_index_low_metallicity]);
        }
        if (star_properties->feedback.lifetimes.dyingtime[metallicity_index + 1][i] >=
                log10_age_yr &&
            time_index_high_metallicity == -1) {
          time_index_high_metallicity = i;
          time_offset_high_metallicity =
              (log10_age_yr - star_properties->feedback.lifetimes
                                .dyingtime[metallicity_index + 1][time_index_high_metallicity]) /
              (star_properties->feedback.lifetimes
                   .dyingtime[metallicity_index + 1][time_index_high_metallicity + 1] -
               star_properties->feedback.lifetimes
                   .dyingtime[metallicity_index + 1][time_index_high_metallicity]);
        }
      }

      mass_low_metallicity =
          interpolate_1d(star_properties->feedback.lifetimes.mass, time_index_low_metallicity, time_offset_low_metallicity);
      mass_high_metallicity =
          interpolate_1d(star_properties->feedback.lifetimes.mass, time_index_high_metallicity, time_offset_high_metallicity);

      mass = (1.0 - metallicity_offset) * mass_low_metallicity + metallicity_offset * mass_high_metallicity;
      break;

    default:
      error("stellar lifetimes not defined\n");
  }

  /* Check that we haven't killed too many stars */
  if (mass > star_properties->feedback.imf_max_mass_msun) mass = star_properties->feedback.imf_max_mass_msun;

  return mass;
}

/**
 * @brief Calculate lifetime of star poputlation in Gyr. Lifetime model is specified by stellar_lifetime_flag read in from yml file. 
 *
 * @param mass
 * @param metallicity
 * @param star_properties the #stars_props data struct
 */
inline static float lifetime_in_Gyr(float mass, float metallicity,
                                    const struct stars_props* restrict
                                        star_properties) {

  double time_Gyr = 0, mass_offset, metallicity_offset;

  int mass_index, metallicity_index;

  switch (star_properties->feedback.stellar_lifetime_flag) {
    case 0:
      /* PM93 (Padovani & Matteucci 1993) */

      if (mass <= 0.6)
        time_Gyr = 160.0;
      else if (mass <= 6.6)
        time_Gyr =
            pow(10.0, (0.334 - sqrt(1.790 - 0.2232 * (7.764 - log10(mass)))) /
                          0.1116);
      else
        time_Gyr = 1.2 * pow(mass, -1.85) + 0.003;
      break;

    case 1:
      /* MM89 (Maeder & Meynet 1989) */

      if (mass <= 1.3)
        time_Gyr = pow(10.0, -0.6545 * log10(mass) + 1.0);
      else if (mass <= 3.0)
        time_Gyr = pow(10.0, -3.7 * log10(mass) + 1.35);
      else if (mass <= 7.0)
        time_Gyr = pow(10.0, -2.51 * log10(mass) + 0.77);
      else if (mass <= 15.0)
        time_Gyr = pow(10.0, -1.78 * log10(mass) + 0.17);
      else if (mass <= 60.0)
        time_Gyr = pow(10.0, -0.86 * log10(mass) - 0.94);
      else
        time_Gyr = 1.2 * pow(mass, -1.85) + 0.003;
      break;

    case 2:
      /* P98 (Portinari et al. 1998) */

      if (mass <= star_properties->feedback.lifetimes.mass[0]) {
        mass_index = 0;
        mass_offset = 0.0;
      } else if (mass >= star_properties->feedback.lifetimes
                             .mass[star_properties->feedback.lifetimes.n_mass - 1]) {
        mass_index = star_properties->feedback.lifetimes.n_mass - 2;
        mass_offset = 1.0;
      } else {
        for (mass_index = 0; mass_index < star_properties->feedback.lifetimes.n_mass - 1;
             mass_index++)
          if (star_properties->feedback.lifetimes.mass[mass_index + 1] > mass) break;

        mass_offset = (mass - star_properties->feedback.lifetimes.mass[mass_index]) /
                 (star_properties->feedback.lifetimes.mass[mass_index + 1] -
                  star_properties->feedback.lifetimes.mass[mass_index]);
      }

      if (metallicity <= star_properties->feedback.lifetimes.metallicity[0]) {
        metallicity_index = 0;
        metallicity_offset = 0.0;
      } else if (metallicity >=
                 star_properties->feedback.lifetimes
                     .metallicity[star_properties->feedback.lifetimes.n_z - 1]) {
        metallicity_index = star_properties->feedback.lifetimes.n_z - 2;
        metallicity_offset = 1.0;
      } else {
        for (metallicity_index = 0; metallicity_index < star_properties->feedback.lifetimes.n_z - 1;
             metallicity_index++)
          if (star_properties->feedback.lifetimes.metallicity[metallicity_index + 1] >
              metallicity)
            break;

        metallicity_offset = (metallicity -
                   star_properties->feedback.lifetimes.metallicity[metallicity_index]) /
                  (star_properties->feedback.lifetimes.metallicity[metallicity_index + 1] -
                   star_properties->feedback.lifetimes.metallicity[metallicity_index]);
      }

      /* time in Gyr */
      time_Gyr =
          exp(M_LN10 * interpolate_2d(star_properties->feedback.lifetimes.dyingtime,
                                   metallicity_index, mass_index, metallicity_offset, mass_offset)) /
          1.0e9;

      break;

    default:
      error("stellar lifetimes not defined");
  }

  return time_Gyr;
}

#endif
