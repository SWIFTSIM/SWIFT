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

static const float log_imf_min_solar_mass = -1;
static const float log_imf_max_solar_mass = 2;

// should be somewhere else maybe?
static const int N_imf_mass_bins = 200;
static const int gamma_SNIa = 2;

// do we need doubles in signature?
inline static void determine_imf_bins(
    double log_min_dying_mass, double log_max_dying_mass, int *ilow, int *ihigh,
    const struct stars_props *restrict star_properties) {
  int i1, i2;

  if (log_min_dying_mass < star_properties->imf_mass_bin_log10[0])
    log_min_dying_mass = star_properties->imf_mass_bin_log10[0];

  if (log_min_dying_mass >
      star_properties->imf_mass_bin_log10[N_imf_mass_bins - 1])
    log_min_dying_mass =
        star_properties->imf_mass_bin_log10[N_imf_mass_bins - 1];

  if (log_max_dying_mass < star_properties->imf_mass_bin_log10[0])
    log_max_dying_mass = star_properties->imf_mass_bin_log10[0];

  if (log_max_dying_mass >
      star_properties->imf_mass_bin_log10[N_imf_mass_bins - 1])
    log_max_dying_mass =
        star_properties->imf_mass_bin_log10[N_imf_mass_bins - 1];

  for (i1 = 0; i1 < N_imf_mass_bins - 2 &&
               star_properties->imf_mass_bin_log10[i1 + 1] < log_min_dying_mass;
       i1++);
    
  for (i2 = 1; i2 < N_imf_mass_bins - 1 &&
               star_properties->imf_mass_bin_log10[i2] < log_max_dying_mass;
       i2++);
    
  *ilow = i1;
  *ihigh = i2;
}

inline static float integrate_imf(
    float log_min_mass, float log_max_mass, float m2, int mode,
    float *stellar_yields, const struct stars_props *restrict star_properties) {

  double result, u, f;

  int ilow, ihigh, index;

  float dlm, dm;

  dlm = star_properties->imf_mass_bin_log10[1] -
        star_properties->imf_mass_bin_log10[0]; /* dlog(m) */

  determine_imf_bins(log_min_mass, log_max_mass, &ilow, &ihigh,
                     star_properties);

  float integrand[N_imf_mass_bins];

  if (mode == 0)
    for (index = ilow; index < ihigh + 1; index++)
      integrand[index] =
          star_properties->imf_by_number[index] *
          star_properties->imf_mass_bin[index]; /* integrate by number */
  else if (mode == 1)
    for (index = ilow; index < ihigh + 1; index++)
      integrand[index] =
          star_properties->imf_by_number[index] *
          star_properties->imf_mass_bin[index] *
          star_properties->imf_mass_bin[index]; /* integrate by mass */
  else if (mode == 2)
    for (index = ilow; index < ihigh + 1; index++){
      integrand[index] =
          stellar_yields[index] * star_properties->imf_by_number[index] *
          star_properties
              ->imf_mass_bin[index]; /* integrate number * yield weighted */
    }
  else if (mode == 3)
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

  /* integrate using trapezoidal rule */
  result = 0;
  for (index = ilow; index < ihigh + 1; index++) result += integrand[index];

  result = result - 0.5 * integrand[ilow] - 0.5 * integrand[ihigh];

  /* correct first bin */
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

  result *= dlm * log(10.0); /* log(10) since mass function tabulated as
                                function of log_10(mass) */

  return result;
}

inline static void init_imf(struct stars_props *restrict star_properties) {

  // ALEXEI: use better names for solar_mass, log_solar_mass
  float norm = 0, solar_mass, log_solar_mass;
  const float dlm = (log_imf_max_solar_mass - log_imf_min_solar_mass) /
                    (double)(N_imf_mass_bins - 1);

  if (posix_memalign((void **)&star_properties->imf_by_number,
                     SWIFT_STRUCT_ALIGNMENT,
                     N_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");
  if (posix_memalign((void **)&star_properties->imf_mass_bin,
                     SWIFT_STRUCT_ALIGNMENT,
                     N_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");
  if (posix_memalign((void **)&star_properties->imf_mass_bin_log10,
                     SWIFT_STRUCT_ALIGNMENT,
                     N_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");
  if (posix_memalign((void **)&star_properties->imf_by_number1,
                     SWIFT_STRUCT_ALIGNMENT,
                     N_imf_mass_bins * sizeof(float)) != 0)
    error("Failed to allocate IMF bins table");

  float dummy_stellar_fields;

  /* Power-law IMF (Salpeter for IMF_EXPONENT = 2.35) */
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
  /* Chabrier 2003 */
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

  norm = integrate_imf(log_imf_min_solar_mass, log_imf_max_solar_mass, 0.0, 1,
                       &dummy_stellar_fields, star_properties);

  for (int i = 0; i < N_imf_mass_bins; i++)
    star_properties->imf_by_number[i] /= norm;
}

#endif
