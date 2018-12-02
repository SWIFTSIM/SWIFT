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

/*
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
 * @param cosmo #cosmology structure
 */
__attribute__((always_inline)) INLINE double eagle_convert_u_to_temp(
    double log_10_u, float *dT_du, int n_h_i, int He_i, float d_n_h, float d_He,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo) {

  int u_i;
  float d_u, log_10_T, log_10_T_high, log_10_T_low;

  get_index_1d(cooling->Therm, cooling->N_Temp, log_10_u, &u_i, &d_u);

  /* Interpolate temperature table to return temperature for current
   * internal energy (use 3D interpolation for high redshift table,
   * otherwise 4D) */
  if (cosmo->z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
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
  if (cosmo->z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
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
  if (cosmo->z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
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
 * EAGLE cooling tables which depend on redshift, temperature,
 * hydrogen number density, helium fraction and metal abundance. Since
 * only the temperature changes when cooling a given particle, the
 * redshift, hydrogen number density and helium fraction indices and
 * offsets passed in. Also calculates derivative of cooling rate with
 * respect to internal energy, which is used in Newton's method for
 * integrating the cooling equation.
 *
 *
 * @param log_10_u Log base 10 of internal energy
 * @param dlambda_du Pointer to value to be set to derivative
 * of cooling rate with respect to internal energy. If set to
 * NULL do not compute derivative.
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param phys_const Physical constants structure
 * @param element_lambda Pointer to array for printing contribution
 * to cooling rate from each of the metals. This is used only for
 * testing and is set to non-NULL when this function is called
 * in eagle_print_metal_cooling_rate. Setting element_lambda to NULL
 * will skip writing to this array (as is done in eagle_cooling_rate,
 * when running SWIFT).
 * @param solar_ratio Array of ratios of particle metal abundances
 * to solar metal abundances
 */
INLINE static double eagle_metal_cooling_rate(
    double log_10_u, double *dlambda_du, int n_h_i, float d_n_h, int He_i,
    float d_He, const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo, const struct phys_const *phys_const,
    double *element_lambda, float *solar_ratio) {
  double z = cosmo->z;
  double cooling_rate = 0.0, temp_lambda, h_plus_he_electron_abundance;
  double solar_electron_abundance;

  /* convert Hydrogen mass fraction in Hydrogen number density */
  const float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];
  const double n_h = hydro_get_physical_density(p, cosmo) * XH /
                     phys_const->const_proton_mass *
                     cooling->number_density_scale;

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

  temp = eagle_convert_u_to_temp(log_10_u, &dT_du, n_h_i, He_i, d_n_h, d_He,
                                 cooling, cosmo);
  get_index_1d(cooling->Temp, cooling->N_Temp, temp, &temp_i, &d_temp);
  float delta_T = exp(M_LN10 * cooling->Temp[temp_i + 1]) -
                  exp(M_LN10 * cooling->Temp[temp_i]);

  /* ------------------ */
  /* Metal-free cooling */
  /* ------------------ */

  /* contribution to cooling and electron abundance from H, He. */
  if (z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
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

  if (z > cooling->Redshifts[cooling->N_Redshifts - 1] ||
      z > cooling->reionisation_redshift) {

    temp_lambda = -cooling->compton_rate_cgs *
                  (temp - cooling->T_CMB_0 * (1 + z)) * (1 + z) * (1 + z) *
                  (1 + z) * (1 + z) * h_plus_he_electron_abundance / n_h;
    cooling_rate += temp_lambda;
    if (element_lambda != NULL) element_lambda[1] = temp_lambda;
  }

  /* ------------- */
  /* Metal cooling */
  /* ------------- */

  /* for each element the cooling rate is multiplied by the ratio of H, He
   * electron abundance to solar electron abundance then by the ratio of the
   * particle metal abundance to solar metal abundance. */

  if (z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
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
 * @param logu Natural log of internal energy
 * @param dLambdaNet_du Pointer to gradient of cooling rate (set in
 * eagle_metal_cooling_rate)
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param p Particle structure
 * @param cooling #cooling_function_data structure
 * @param cosmo #cosmology structure
 * @param phys_const #phys_const structure
 * @param abundance_ratio Ratio of element abundance to solar
 */
INLINE static double eagle_cooling_rate(
    double logu, double *dLambdaNet_du, int n_h_i, float d_n_h, int He_i,
    float d_He, const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo, const struct phys_const *phys_const,
    float *abundance_ratio) {

  /* set element_lambda to NULL so will not print file of
   * contributions to cooling from each element */
  double *element_lambda = NULL;
  double lambda_net = 0.0;

  /* calculate cooling rate and set dLambdaNet_du */
  lambda_net = eagle_metal_cooling_rate(
      logu / M_LN10, dLambdaNet_du, n_h_i, d_n_h, He_i, d_He, p, cooling, cosmo,
      phys_const, element_lambda, abundance_ratio);

  return lambda_net;
}
