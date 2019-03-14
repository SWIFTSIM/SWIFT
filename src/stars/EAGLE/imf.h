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

static const float imf_max_mass_msun = 100.f;
static const float imf_min_mass_msun = 0.1;
static const float log_imf_min_solar_mass = -1;
static const float log_imf_max_solar_mass = 2;

// should be somewhere else maybe?
static const int N_imf_mass_bins = 200;
static const int gamma_SNIa = 2;

/**
 * @brief determine which IMF mass bins the upper and lower input mass bounds belong to
 *
 * @param log_min_mass Lower mass bound
 * @param log_max_mass Upper mass bound
 * @param ilow pointer to index of IMF mass bin containing log_min_mass
 * @param ihigh pointer to index of IMF mass bin containing log_max_mass
 * @param star_properties the #stars_props data struct */
// do we need doubles in signature?
inline static void determine_imf_bins(
    double log_min_mass, double log_max_mass, int *ilow, int *ihigh,
    const struct stars_props *restrict star_properties) {
  int i1, i2;

  /* Check whether lower mass is within the IMF mass bin range */
  if (log_min_mass < star_properties->imf_mass_bin_log10[0])
    log_min_mass = star_properties->imf_mass_bin_log10[0];
  if (log_min_mass >
      star_properties->imf_mass_bin_log10[N_imf_mass_bins - 1])
    log_min_mass =
        star_properties->imf_mass_bin_log10[N_imf_mass_bins - 1];

  /* Check whether upper mass is within the IMF mass bin range */
  if (log_max_mass < star_properties->imf_mass_bin_log10[0])
    log_max_mass = star_properties->imf_mass_bin_log10[0];
  if (log_max_mass >
      star_properties->imf_mass_bin_log10[N_imf_mass_bins - 1])
    log_max_mass =
        star_properties->imf_mass_bin_log10[N_imf_mass_bins - 1];

  /* Find mass bin to which lower mass belongs (ALEXEI: Can we do this in a better way?) */
  for (i1 = 0; i1 < N_imf_mass_bins - 2 &&
               star_properties->imf_mass_bin_log10[i1 + 1] < log_min_mass;
       i1++);
    
  /* Find mass bin to which upper mass belongs (ALEXEI: Can we do this in a better way?) */
  for (i2 = 1; i2 < N_imf_mass_bins - 1 &&
               star_properties->imf_mass_bin_log10[i2] < log_max_mass;
       i2++);
    
  *ilow = i1;
  *ihigh = i2;
}

/**
 * @brief Integrate the IMF between a minimum and maximum mass using the trapezoidal rule. The IMF may be weighted by various quantities, as specified by the variable, mode, including an input array, stellar_yields.
 *
 * @param log_min_mass log10 mass lower integration bound
 * @param log_max_mass log10 mass upper integration bound
 * @param m2 ALEXEI: I have no idea. Check if it is ever used
 * @param mode Flag to specify weighting used for integrating IMF
 * @param stellar_yields Array of weights based on yields. Used with mode=2
 * @param star_properties the #stars_props data structure
 */
inline static float integrate_imf(
    float log_min_mass, float log_max_mass, float m2, int mode,
    float *stellar_yields, const struct stars_props *restrict star_properties) {

  // ALEXEI: come up with better names for indices
  double result, u, f;

  int ilow, ihigh, index;

  float dlm, dm;

  dlm = star_properties->imf_mass_bin_log10[1] -
        star_properties->imf_mass_bin_log10[0]; /* dlog(m) */

  /* Determine bins to integrate over based on integration bounds */
  determine_imf_bins(log_min_mass, log_max_mass, &ilow, &ihigh,
                     star_properties);

  // ALEXEI: rewrite this function so that adding on the fly and don't need to use the array.
  float integrand[N_imf_mass_bins];

  /* Add up the contribution from each of the IMF mass bins */
  if (mode == 0)
    /* Integrate IMF on its own */
    for (index = ilow; index < ihigh + 1; index++)
      integrand[index] =
          star_properties->imf_by_number[index] *
          star_properties->imf_mass_bin[index]; 
  else if (mode == 1)
    /* Integrate IMF weighted by mass */
    for (index = ilow; index < ihigh + 1; index++)
      integrand[index] =
          star_properties->imf_by_number[index] *
          star_properties->imf_mass_bin[index] *
          star_properties->imf_mass_bin[index]; 
  else if (mode == 2)
    /* Integrate IMF weighted by yields */
    for (index = ilow; index < ihigh + 1; index++)
      integrand[index] =
          stellar_yields[index] * star_properties->imf_by_number[index] *
          star_properties
              ->imf_mass_bin[index]; 
  else if (mode == 3)
    /* ALEXEI: should we keep this? */
    for (index = ilow; index < ihigh + 1; index++) {
      u = m2 / star_properties->imf_mass_bin[index];
      f = pow(2.0, gamma_SNIa + 1) * (gamma_SNIa + 1) * pow(u, gamma_SNIa);
      integrand[index] =
          f * star_properties->imf_by_number[index] /
          star_properties->imf_mass_bin[index]; /* integrate number * f(u) / M
                                                   ... type Ia SN */
    }
  else {
    error("invalid mode in integrate_imf = %d\n", mode);
  }

  /* integrate using trapezoidal rule (ALEXEI: Remove when adding on the fly)*/
  result = 0;
  for (index = ilow; index < ihigh + 1; index++) result += integrand[index];

  /* Update end bins since contribution was overcounted when summing up all entries */
  result = result - 0.5 * integrand[ilow] - 0.5 * integrand[ihigh];

  /* correct first bin (ALEXEI: why are we doing this? (just to clarify comment))*/
  dm = (log_min_mass - star_properties->imf_mass_bin_log10[ilow]) / dlm;

  if (dm < 0.5)
    result -= dm * integrand[ilow];
  else {
    result -= 0.5 * integrand[ilow];
    result -= (dm - 0.5) * integrand[ilow + 1];
  }

  /* correct last bin */
  dm = (log_max_mass - star_properties->imf_mass_bin_log10[ihigh - 1]) / dlm;

  if (dm < 0.5) {
    result -= 0.5 * integrand[ihigh];
    result -= (0.5 - dm) * integrand[ihigh - 1];
  } else {
    result -= (1 - dm) * integrand[ihigh];
  }

  /* The IMF is tabulated in log10 so undo to get regular units */
  result *= dlm * log(10.0); 

  return result;
}

/**
 * @brief Allocate and compute IMF
 *
 * @param star_properties #stars_props data structure */
inline static void init_imf(struct stars_props *restrict star_properties) {

  // ALEXEI: use better names for solar_mass, log_solar_mass
  float norm = 0, solar_mass, log_solar_mass;

  /* Compute size of mass bins in log10 space */
  const float dlm = (log_imf_max_solar_mass - log_imf_min_solar_mass) /
                    (double)(N_imf_mass_bins - 1);

  /* Allocate IMF array */
  if (posix_memalign((void **)&star_properties->imf_by_number,
                     SWIFT_STRUCT_ALIGNMENT,
                     N_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");

  /* Allocate array to store IMF mass bins */
  if (posix_memalign((void **)&star_properties->imf_mass_bin,
                     SWIFT_STRUCT_ALIGNMENT,
                     N_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");

  /* Allocate array to store IMF mass bins in log10 space */
  if (posix_memalign((void **)&star_properties->imf_mass_bin_log10,
                     SWIFT_STRUCT_ALIGNMENT,
                     N_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");

  /* Set Power-law IMF (Salpeter for IMF_EXPONENT = 2.35) */
  if (strcmp(star_properties->IMF_Model, "PowerLaw") == 0) {
    if (star_properties->IMF_Exponent < 0) {
      error("imf_exponent is supposed to be > 0\n");
    }

    for (int i = 0; i < N_imf_mass_bins; i++) {
      log_solar_mass = log_imf_min_solar_mass + i * dlm;
      solar_mass = exp(M_LN10 * log_solar_mass);

      star_properties->IMF_Exponent = 2.35;
      star_properties->imf_by_number[i] =
          pow(solar_mass, -star_properties->IMF_Exponent);

      star_properties->imf_mass_bin[i] = solar_mass;
      star_properties->imf_mass_bin_log10[i] = log_solar_mass;
    }
  }
  /* Set IMF from Chabrier 2003 */
  else if (strcmp(star_properties->IMF_Model, "Chabrier") == 0) {
    for (int i = 0; i < N_imf_mass_bins; i++) {
      log_solar_mass = log_imf_min_solar_mass + i * dlm;
      solar_mass = exp(M_LN10 * log_solar_mass);

      if (solar_mass > 1.0)
        star_properties->imf_by_number[i] = 0.237912 * pow(solar_mass, -2.3);
      else
        star_properties->imf_by_number[i] =
            0.852464 *
            exp((log10(solar_mass) - log10(0.079)) *
                (log10(solar_mass) - log10(0.079)) / (-2.0 * pow(0.69, 2))) /
            solar_mass;

      star_properties->imf_mass_bin[i] = solar_mass;
      star_properties->imf_mass_bin_log10[i] = log_solar_mass;
    }
  } else {
    error("Invalid IMF model %s. Valid models are: PowerLaw and Chabrier\n",
          star_properties->IMF_Model);
  }

  /* Normalize the IMF */
  norm = integrate_imf(log_imf_min_solar_mass, log_imf_max_solar_mass, 0.0, 1,
                       /*(stellar_yields=)*/ NULL, star_properties);

  for (int i = 0; i < N_imf_mass_bins; i++)
    star_properties->imf_by_number[i] /= norm;
}

/**
 * @brief Calculate mass (in solar masses) of stars that died from the star particle's birth up to its current age (in Gyr)
 *
 * @param age_Gyr age of star in Gyr
 * @param metallicity Star's metallicity
 * @param star_properties the #stars_props data structure
 */
inline static float dying_mass_msun(float age_Gyr, float metallicity,
                                    const struct stars_props* restrict
                                        star_properties) {

  // ALEXEI: check out units for all these quantities and name them accordingly
  float mass = 0, d_metal, d_time1 = 0, d_time2 = 0, log_age_yr, mass1, mass2;

  int metal_index, i, index_time1 = -1, index_time2 = -1;

  /* Calculate star mass that dies based on one of several analytical models */
  // ALEXEI: What do we do with the constants here?
  switch (star_properties->stellar_lifetime_flag) {
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
        mass = imf_max_mass_msun;
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
        mass = imf_max_mass_msun;
      break;

    case 2:
      /* portinari et al. 1998 */

      if (age_Gyr <= 0) {
        mass = imf_max_mass_msun;
        break;
      }

      log_age_yr = log10(age_Gyr * 1.e9);

      // ALEXEI: Simply copied from EAGLE, can we simplify this whole thing somehow to be more readable?
      if (metallicity <= star_properties->lifetimes.metallicity[0]) {
        metal_index = 0;
        d_metal = 0.0;
      } else if (metallicity >=
                 star_properties->lifetimes
                     .metallicity[star_properties->lifetimes.n_z - 1]) {
        metal_index = star_properties->lifetimes.n_z - 2;
        d_metal = 1.0;
      } else {
        for (metal_index = 0; metal_index < star_properties->lifetimes.n_z - 1;
             metal_index++)
          if (star_properties->lifetimes.metallicity[metal_index + 1] >
              metallicity)
            break;

        d_metal = (metallicity -
                   star_properties->lifetimes.metallicity[metal_index]) /
                  (star_properties->lifetimes.metallicity[metal_index + 1] -
                   star_properties->lifetimes.metallicity[metal_index]);
      }

      if (log_age_yr >= star_properties->lifetimes.dyingtime[metal_index][0]) {
        index_time1 = 0;
        d_time1 = 0.0;
      } else if (log_age_yr <=
                 star_properties->lifetimes
                     .dyingtime[metal_index]
                               [star_properties->lifetimes.n_mass - 1]) {
        index_time1 = star_properties->lifetimes.n_mass - 2;
        d_time1 = 1.0;
      }

      if (log_age_yr >=
          star_properties->lifetimes.dyingtime[metal_index + 1][0]) {
        index_time2 = 0;
        d_time2 = 0.0;
      } else if (log_age_yr <=
                 star_properties->lifetimes
                     .dyingtime[metal_index + 1]
                               [star_properties->lifetimes.n_mass - 1]) {
        index_time2 = star_properties->lifetimes.n_mass - 2;
        d_time2 = 1.0;
      }

      i = star_properties->lifetimes.n_mass;
      while (i >= 0 && (index_time1 == -1 || index_time2 == -1)) {
        i--;
        if (star_properties->lifetimes.dyingtime[metal_index][i] >=
                log_age_yr &&
            index_time1 == -1) {
          index_time1 = i;
          d_time1 =
              (log_age_yr -
               star_properties->lifetimes.dyingtime[metal_index][index_time1]) /
              (star_properties->lifetimes
                   .dyingtime[metal_index][index_time1 + 1] -
               star_properties->lifetimes.dyingtime[metal_index][index_time1]);
        }
        if (star_properties->lifetimes.dyingtime[metal_index + 1][i] >=
                log_age_yr &&
            index_time2 == -1) {
          index_time2 = i;
          d_time2 =
              (log_age_yr - star_properties->lifetimes
                                .dyingtime[metal_index + 1][index_time2]) /
              (star_properties->lifetimes
                   .dyingtime[metal_index + 1][index_time2 + 1] -
               star_properties->lifetimes
                   .dyingtime[metal_index + 1][index_time2]);
        }
      }

      mass1 =
          interpolate_1d(star_properties->lifetimes.mass, index_time1, d_time1);
      mass2 =
          interpolate_1d(star_properties->lifetimes.mass, index_time2, d_time2);

      mass = (1.0 - d_metal) * mass1 + d_metal * mass2;
      break;

    default:
      error("stellar lifetimes not defined\n");
  }

  /* Check that we haven't killed too many stars */
  if (mass > imf_max_mass_msun) mass = imf_max_mass_msun;

  return mass;
}

/**
 * @brief
 *
 * @param mass
 * @param metallicity
 * @param star_properties the #stars_props data struct
 */
inline static float lifetime_in_Gyr(float mass, float metallicity,
                                    const struct stars_props* restrict
                                        star_properties) {

  double time = 0, d_mass, d_metal;

  int mass_index, metal_index;

  switch (star_properties->stellar_lifetime_flag) {
    case 0:
      /* PM93 (Padovani & Matteucci 1993) */

      if (mass <= 0.6)
        time = 160.0;
      else if (mass <= 6.6)
        time =
            pow(10.0, (0.334 - sqrt(1.790 - 0.2232 * (7.764 - log10(mass)))) /
                          0.1116);
      else
        time = 1.2 * pow(mass, -1.85) + 0.003;
      break;

    case 1:
      /* MM89 (Maeder & Meynet 1989) */

      if (mass <= 1.3)
        time = pow(10.0, -0.6545 * log10(mass) + 1.0);
      else if (mass <= 3.0)
        time = pow(10.0, -3.7 * log10(mass) + 1.35);
      else if (mass <= 7.0)
        time = pow(10.0, -2.51 * log10(mass) + 0.77);
      else if (mass <= 15.0)
        time = pow(10.0, -1.78 * log10(mass) + 0.17);
      else if (mass <= 60.0)
        time = pow(10.0, -0.86 * log10(mass) - 0.94);
      else
        time = 1.2 * pow(mass, -1.85) + 0.003;
      break;

    case 2:
      /* P98 (Portinari et al. 1998) */

      if (mass <= star_properties->lifetimes.mass[0]) {
        mass_index = 0;
        d_mass = 0.0;
      } else if (mass >= star_properties->lifetimes
                             .mass[star_properties->lifetimes.n_mass - 1]) {
        mass_index = star_properties->lifetimes.n_mass - 2;
        d_mass = 1.0;
      } else {
        for (mass_index = 0; mass_index < star_properties->lifetimes.n_mass - 1;
             mass_index++)
          if (star_properties->lifetimes.mass[mass_index + 1] > mass) break;

        d_mass = (mass - star_properties->lifetimes.mass[mass_index]) /
                 (star_properties->lifetimes.mass[mass_index + 1] -
                  star_properties->lifetimes.mass[mass_index]);
      }

      if (metallicity <= star_properties->lifetimes.metallicity[0]) {
        metal_index = 0;
        d_metal = 0.0;
      } else if (metallicity >=
                 star_properties->lifetimes
                     .metallicity[star_properties->lifetimes.n_z - 1]) {
        metal_index = star_properties->lifetimes.n_z - 2;
        d_metal = 1.0;
      } else {
        for (metal_index = 0; metal_index < star_properties->lifetimes.n_z - 1;
             metal_index++)
          if (star_properties->lifetimes.metallicity[metal_index + 1] >
              metallicity)
            break;

        d_metal = (metallicity -
                   star_properties->lifetimes.metallicity[metal_index]) /
                  (star_properties->lifetimes.metallicity[metal_index + 1] -
                   star_properties->lifetimes.metallicity[metal_index]);
      }

      /* time in Gyr */
      time =
          exp(M_LN10 * interpolate_2d(star_properties->lifetimes.dyingtime,
                                   metal_index, mass_index, d_metal, d_mass)) /
          1.0e9;

      break;

    default:
      error("stellar lifetimes not defined");
  }

  return time;
}

#endif
