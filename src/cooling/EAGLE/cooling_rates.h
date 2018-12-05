/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "interpolate.h"

/**
 * @brief calculates heating due to helium reionization
 *
 * @param z redshift
 * @param delta_z change in redshift over timestep
 * @param cooling the #cooling_function_data struct
 */
__attribute__((always_inline)) INLINE double
eagle_helium_reionization_extraheat(double z, double delta_z,
                                    const struct cooling_function_data *restrict
                                        cooling) {

#ifdef SWIFT_DEBUG_CHECKS
  if (delta_z > 0.f) error("Invalid value for delta_z. Should be negative.");
#endif

  /* Recover the values we need */
  const double z_centre = cooling->He_reion_z_centre;
  const double z_sigma = cooling->He_reion_z_sigma;
  const double heat_cgs = cooling->He_reion_ev_pH;

  // MATTHIEU: to do: Optimize this.

  double extra_heat;

  /* Integral of the Gaussian between z and z - delta_z */
  extra_heat = erf((z - delta_z - z_centre) / (M_SQRT2 * z_sigma));
  extra_heat -= erf((z - z_centre) / (M_SQRT2 * z_sigma));

  /* Multiply by the normalisation factor */
  extra_heat *= heat_cgs * 0.5;

  return extra_heat;
}

/**
 * @brief Interpolates temperature from internal energy based on table and
 * calculates the size of the internal energy cell for the specified
 * internal energy. Returns log base 10 of temperature.
 *
 * @param log_10_u Log base 10 of internal energy
 * @param dT_du Pointer to rate of change of log_10(temperature) with internal
 * energy
 * @param z_i Redshift index
 * @param n_h_i Hydrogen number density index
 * @param He_i Helium fraction index
 * @param d_z Redshift offset
 * @param d_n_h Hydrogen number density offset
 * @param d_He Helium fraction offset
 * @param cooling #cooling_function_data structure
 */
__attribute__((always_inline)) INLINE double eagle_convert_u_to_temp(
    double log_10_u, float redshift, float *dT_du, int n_h_i, int He_i,
    float d_n_h, float d_He,
    const struct cooling_function_data *restrict cooling) {

  int u_i;
  float d_u, log_10_T, log_10_T_high, log_10_T_low;

  get_index_1d(cooling->Therm, cooling->N_Temp, log_10_u, &u_i, &d_u);

  /* Interpolate temperature table to return temperature for current
   * internal energy (use 3D interpolation for high redshift table,
   * otherwise 4D) */
  if (redshift > cooling->Redshifts[cooling->N_Redshifts - 1]) {
    log_10_T = interpolate_3d(cooling->table.temperature, n_h_i, He_i, u_i,
                              d_n_h, d_He, d_u, cooling->N_nH, cooling->N_He,
                              cooling->N_Temp);
  } else {
    log_10_T = interpolate_4d(cooling->table.temperature, 0, n_h_i, He_i, u_i,
                              cooling->dz, d_n_h, d_He, d_u, 2, cooling->N_nH,
                              cooling->N_He, cooling->N_Temp);
  }

  /* Interpolate temperature table to return temperature for internal energy
   * at grid point above current internal energy for computing dT_du used for
   * calculation of dlambda_du in cooling.c (use 3D interpolation for high
   * redshift table, otherwise 4D) */
  if (redshift > cooling->Redshifts[cooling->N_Redshifts - 1]) {
    log_10_T_high = interpolate_3d(cooling->table.temperature, n_h_i, He_i, u_i,
                                   d_n_h, d_He, 1.0, cooling->N_nH,
                                   cooling->N_He, cooling->N_Temp);
  } else {
    log_10_T_high = interpolate_4d(
        cooling->table.temperature, 0, n_h_i, He_i, u_i, cooling->dz, d_n_h,
        d_He, 1.0, 2, cooling->N_nH, cooling->N_He, cooling->N_Temp);
  }

  /* Interpolate temperature table to return temperature for internal energy
   * at grid point below current internal energy for computing dT_du used for
   * calculation of dlambda_du in cooling.c (use 3D interpolation for high
   * redshift table, otherwise 4D) */
  if (redshift > cooling->Redshifts[cooling->N_Redshifts - 1]) {
    log_10_T_low = interpolate_3d(cooling->table.temperature, n_h_i, He_i, u_i,
                                  d_n_h, d_He, 0.0, cooling->N_nH,
                                  cooling->N_He, cooling->N_Temp);
  } else {
    log_10_T_low = interpolate_4d(
        cooling->table.temperature, 0, n_h_i, He_i, u_i, cooling->dz, d_n_h,
        d_He, 0.0, 2, cooling->N_nH, cooling->N_He, cooling->N_Temp);
  }

  /* Calculate dT/du */
  float delta_u =
      exp(cooling->Therm[u_i + 1] * M_LN10) - exp(cooling->Therm[u_i] * M_LN10);
  *dT_du = (exp(M_LN10 * log_10_T_high) - exp(M_LN10 * log_10_T_low)) / delta_u;

  return log_10_T;
}

/**
 * @brief Calculates cooling rate for given internal energy by interpolating
 * EAGLE cooling tables.
 *
 * The tables depend on redshift, temperature, hydrogen number
 * density, helium fraction and metal abundance. Since only the
 * temperature changes when cooling a given particle, the redshift,
 * hydrogen number density and helium fraction indices and offsets
 * passed in.
 *
 * If the arguement dlambda_du is non-NULL, the routine also
 * calculates derivative of cooling rate with respect to internal
 * energy.
 *
 * If the argument element_lambda is non-NULL, the routine also
 * returns the cooling rate per element in the array.
 *
 * @param log10_u_cgs Log base 10 of internal energy per unit mass in CGS units.
 * @param redshift The current redshift
 * @param n_H_cgs The Hydrogen number density in CGS units.
 * @param solar_ratio Array of ratios of particle metal abundances
 * to solar metal abundances
 *
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param cooling Cooling data structure
 * @param phys_const Physical constants structure
 *
 * @param dlambda_du (return) Derivative of the cooling rate with respect to u.
 * @param element_lambda (return) Cooling rate from each element
 *
 * @return The cooling rate
 */
INLINE static double eagle_metal_cooling_rate(
    double log10_u_cgs, double redshift, double n_H_cgs,
    const float solar_ratio[chemistry_element_count + 2], int n_h_i,
    float d_n_h, int He_i, float d_He,
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *phys_const, double *dlambda_du,
    double *element_lambda) {

  double cooling_rate = 0.0, temp_lambda, h_plus_he_electron_abundance;
  double solar_electron_abundance;

  /* used for calculating dlambda_du */
  double temp_lambda_high = 0, temp_lambda_low = 0;
  double h_plus_he_electron_abundance_high = 0;
  double h_plus_he_electron_abundance_low = 0;
  double solar_electron_abundance_high = 0;
  double solar_electron_abundance_low = 0;
  double elem_cool_low = 0, elem_cool_high = 0;
  float dT_du;

  /* counter, temperature index, value, and offset */
  int i, temp_i;
  double temp;
  float d_temp;

  /* interpolate to get temperature of particles, find where we are in
   * the temperature table. */

  temp = eagle_convert_u_to_temp(log10_u_cgs, redshift, &dT_du, n_h_i, He_i,
                                 d_n_h, d_He, cooling);
  get_index_1d(cooling->Temp, cooling->N_Temp, temp, &temp_i, &d_temp);

  const float delta_T = exp(M_LN10 * cooling->Temp[temp_i + 1]) -
                        exp(M_LN10 * cooling->Temp[temp_i]);

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

  /* contribution to cooling and electron abundance from H, He. */
  if (redshift > cooling->Redshifts[cooling->N_Redshifts - 1]) {

    /* If we're using the high redshift tables then we don't interpolate
     * in redshift */
    temp_lambda = interpolate_3d(cooling->table.H_plus_He_heating, n_h_i, He_i,
                                 temp_i, d_n_h, d_He, d_temp, cooling->N_nH,
                                 cooling->N_He, cooling->N_Temp);

    h_plus_he_electron_abundance = interpolate_3d(
        cooling->table.H_plus_He_electron_abundance, n_h_i, He_i, temp_i, d_n_h,
        d_He, d_temp, cooling->N_nH, cooling->N_He, cooling->N_Temp);

    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du. Pass in NULL pointer for
     * dlambda_du in order to skip */
    if (dlambda_du != NULL) {
      temp_lambda_high = interpolate_3d(
          cooling->table.H_plus_He_heating, n_h_i, He_i, temp_i, d_n_h, d_He,
          1.0, cooling->N_nH, cooling->N_He, cooling->N_Temp);
      temp_lambda_low = interpolate_3d(
          cooling->table.H_plus_He_heating, n_h_i, He_i, temp_i, d_n_h, d_He,
          0.0, cooling->N_nH, cooling->N_He, cooling->N_Temp);
      h_plus_he_electron_abundance_high = interpolate_3d(
          cooling->table.H_plus_He_electron_abundance, n_h_i, He_i, temp_i,
          d_n_h, d_He, 1.0, cooling->N_nH, cooling->N_He, cooling->N_Temp);
      h_plus_he_electron_abundance_low = interpolate_3d(
          cooling->table.H_plus_He_electron_abundance, n_h_i, He_i, temp_i,
          d_n_h, d_He, 0.0, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    }
  } else {
    /* Using normal tables, have to interpolate in redshift */
    temp_lambda = interpolate_4d(
        cooling->table.H_plus_He_heating, 0, n_h_i, He_i, temp_i, cooling->dz,
        d_n_h, d_He, d_temp, 2, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    h_plus_he_electron_abundance =
        interpolate_4d(cooling->table.H_plus_He_electron_abundance, 0, n_h_i,
                       He_i, temp_i, cooling->dz, d_n_h, d_He, d_temp, 2,
                       cooling->N_nH, cooling->N_He, cooling->N_Temp);

    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      temp_lambda_high = interpolate_4d(
          cooling->table.H_plus_He_heating, 0, n_h_i, He_i, temp_i, cooling->dz,
          d_n_h, d_He, 1.0, 2, cooling->N_nH, cooling->N_He, cooling->N_Temp);
      temp_lambda_low = interpolate_4d(
          cooling->table.H_plus_He_heating, 0, n_h_i, He_i, temp_i, cooling->dz,
          d_n_h, d_He, 0.0, 2, cooling->N_nH, cooling->N_He, cooling->N_Temp);
      h_plus_he_electron_abundance_high =
          interpolate_4d(cooling->table.H_plus_He_electron_abundance, 0, n_h_i,
                         He_i, temp_i, cooling->dz, d_n_h, d_He, 1.0, 2,
                         cooling->N_nH, cooling->N_He, cooling->N_Temp);
      h_plus_he_electron_abundance_low =
          interpolate_4d(cooling->table.H_plus_He_electron_abundance, 0, n_h_i,
                         He_i, temp_i, cooling->dz, d_n_h, d_He, 0.0, 2,
                         cooling->N_nH, cooling->N_He, cooling->N_Temp);
    }
  }
  cooling_rate += temp_lambda;
  if (dlambda_du != NULL)
    *dlambda_du += (temp_lambda_high - temp_lambda_low) / delta_T * dT_du;

  /* If we're testing cooling rate contributions write to array */
  if (element_lambda != NULL) element_lambda[0] = temp_lambda;

  /* ------------------ */
  /* Compton cooling    */
  /* ------------------ */

  /* inverse Compton cooling is not in collisional table
   * before reionisation so add now */

  if (redshift > cooling->Redshifts[cooling->N_Redshifts - 1] ||
      redshift > cooling->reionisation_redshift) {

    const double zp1 = 1. + redshift;
    const double zp1p4 = zp1 * zp1 * zp1 * zp1;
    const double T_CMB = cooling->T_CMB_0 * zp1;

    temp_lambda = -cooling->compton_rate_cgs * (temp - T_CMB) * zp1p4 *
                  h_plus_he_electron_abundance / n_H_cgs;

    cooling_rate += temp_lambda;
    if (element_lambda != NULL) element_lambda[1] = temp_lambda;
  }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

  /* for each element the cooling rate is multiplied by the ratio of H, He
   * electron abundance to solar electron abundance then by the ratio of the
   * particle metal abundance to solar metal abundance. */

  if (redshift > cooling->Redshifts[cooling->N_Redshifts - 1]) {
    /* If we're using the high redshift tables then we don't interpolate
     * in redshift */
    solar_electron_abundance =
        interpolate_2d(cooling->table.electron_abundance, n_h_i, temp_i, d_n_h,
                       d_temp, cooling->N_nH, cooling->N_Temp);
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      solar_electron_abundance_high =
          interpolate_2d(cooling->table.electron_abundance, n_h_i, temp_i,
                         d_n_h, 1.0, cooling->N_nH, cooling->N_Temp);
      solar_electron_abundance_low =
          interpolate_2d(cooling->table.electron_abundance, n_h_i, temp_i,
                         d_n_h, 0.0, cooling->N_nH, cooling->N_Temp);
    }

    for (i = 0; i < cooling->N_Elements; i++) {
      temp_lambda = interpolate_3d(cooling->table.metal_heating, n_h_i, temp_i,
                                   i, d_n_h, d_temp, 0.0, cooling->N_nH,
                                   cooling->N_Temp, cooling->N_Elements) *
                    (h_plus_he_electron_abundance / solar_electron_abundance) *
                    solar_ratio[i + 2];
      cooling_rate += temp_lambda;

      /* compute values at temperature gridpoints above and below input
       * temperature for calculation of dlambda_du */
      if (dlambda_du != NULL) {
        elem_cool_high = interpolate_3d(
            cooling->table.metal_heating, n_h_i, temp_i, i, d_n_h, 1.0, 0.0,
            cooling->N_nH, cooling->N_Temp, cooling->N_Elements);
        elem_cool_low = interpolate_3d(
            cooling->table.metal_heating, n_h_i, temp_i, i, d_n_h, 0.0, 0.0,
            cooling->N_nH, cooling->N_Temp, cooling->N_Elements);
        *dlambda_du += (elem_cool_high * h_plus_he_electron_abundance_high /
                            solar_electron_abundance_high -
                        elem_cool_low * h_plus_he_electron_abundance_low /
                            solar_electron_abundance_low) /
                       delta_T * dT_du * solar_ratio[i + 2];
      }
      if (element_lambda != NULL) element_lambda[i + 2] = temp_lambda;
    }
  } else {
    /* Using normal tables, have to interpolate in redshift */
    solar_electron_abundance = interpolate_3d(
        cooling->table.electron_abundance, 0, n_h_i, temp_i, cooling->dz, d_n_h,
        d_temp, 2, cooling->N_nH, cooling->N_Temp);
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      solar_electron_abundance_high = interpolate_3d(
          cooling->table.electron_abundance, 0, n_h_i, temp_i, cooling->dz,
          d_n_h, 1.0, 2, cooling->N_nH, cooling->N_Temp);
      solar_electron_abundance_low = interpolate_3d(
          cooling->table.electron_abundance, 0, n_h_i, temp_i, cooling->dz,
          d_n_h, 0.0, 2, cooling->N_nH, cooling->N_Temp);
    }

    for (i = 0; i < cooling->N_Elements; i++) {
      temp_lambda =
          interpolate_4d(cooling->table.metal_heating, 0, n_h_i, temp_i, i,
                         cooling->dz, d_n_h, d_temp, 0.0, 2, cooling->N_nH,
                         cooling->N_Temp, cooling->N_Elements) *
          (h_plus_he_electron_abundance / solar_electron_abundance) *
          solar_ratio[i + 2];
      cooling_rate += temp_lambda;

      /* compute values at temperature gridpoints above and below input
       * temperature for calculation of dlambda_du */
      if (dlambda_du != NULL) {
        elem_cool_high =
            interpolate_4d(cooling->table.metal_heating, 0, n_h_i, temp_i, i,
                           cooling->dz, d_n_h, 1.0, 0.0, 2, cooling->N_nH,
                           cooling->N_Temp, cooling->N_Elements);
        elem_cool_low =
            interpolate_4d(cooling->table.metal_heating, 0, n_h_i, temp_i, i,
                           cooling->dz, d_n_h, 0.0, 0.0, 2, cooling->N_nH,
                           cooling->N_Temp, cooling->N_Elements);
        *dlambda_du += (elem_cool_high * h_plus_he_electron_abundance_high /
                            solar_electron_abundance_high -
                        elem_cool_low * h_plus_he_electron_abundance_low /
                            solar_electron_abundance_low) /
                       delta_T * dT_du * solar_ratio[i + 2];
      }
      if (element_lambda != NULL) element_lambda[i + 2] = temp_lambda;
    }
  }

  return cooling_rate;
}

/**
 * @brief Wrapper function used to calculate cooling rate and dLambda_du.
 * Table indices and offsets for redshift, hydrogen number density and
 * helium fraction are passed it so as to compute them only once per particle.
 *
 * @param log_u_cgs Natural log of internal energy per unit mass in CGS units.
 * @param redshift The current redshift.
 * @param n_H_cgs Hydrogen number density in CGS units.
 * @param abundance_ratio Ratio of element abundance to solar.
 *
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param cooling #cooling_function_data structure
 * @param phys_const #phys_const structure
 *
 * @param dLambdaNet_du (return) Derivative of the cooling rate with respect to
 * u.
 *
 * @return The cooling rate
 */
INLINE static double eagle_cooling_rate(
    double log_u_cgs, double redshift, double n_H_cgs,
    const float abundance_ratio[chemistry_element_count + 2], int n_h_i,
    float d_n_h, int He_i, float d_He,
    const struct cooling_function_data *restrict cooling,
    const struct phys_const *phys_const, double *dLambdaNet_du) {

  return eagle_metal_cooling_rate(
      log_u_cgs / M_LN10, redshift, n_H_cgs, abundance_ratio, n_h_i, d_n_h,
      He_i, d_He, cooling, phys_const, dLambdaNet_du, /*element_lambda=*/NULL);
}
