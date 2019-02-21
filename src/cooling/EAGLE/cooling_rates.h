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
#ifndef SWIFT_EAGLE_COOLING_RATES_H
#define SWIFT_EAGLE_COOLING_RATES_H

#include "../config.h"

/* Local includes. */
#include "cooling_tables.h"
#include "exp10.h"
#include "interpolate.h"

/**
 * @brief Compute ratio of mass fraction to solar mass fraction
 * for each element carried by a given particle.
 *
 * The solar abundances are taken from the tables themselves.
 *
 * The EAGLE chemistry model does not track S and Ca. We assume
 * that their abundance with respect to solar is the same as
 * the ratio for Si.
 * We optionally apply a correction if the user asked for a different
 * ratio.
 *
 * We also re-order the elements such that they match the order of the
 * tables. This is [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe].
 *
 * The solar abundances table (from the cooling struct) is arranged as
 * [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe].
 *
 * @param p Pointer to #part struct.
 * @param cooling #cooling_function_data struct.
 * @param ratio_solar (return) Array of ratios to solar abundances.
 */
__attribute__((always_inline)) INLINE void abundance_ratio_to_solar(
    const struct part *p, const struct cooling_function_data *cooling,
    float ratio_solar[chemistry_element_count + 2]) {

  ratio_solar[0] = p->chemistry_data.metal_mass_fraction[chemistry_element_H] *
                   cooling->SolarAbundances_inv[0 /* H */];

  ratio_solar[1] = p->chemistry_data.metal_mass_fraction[chemistry_element_He] *
                   cooling->SolarAbundances_inv[1 /* He */];

  ratio_solar[2] = p->chemistry_data.metal_mass_fraction[chemistry_element_C] *
                   cooling->SolarAbundances_inv[2 /* C */];

  ratio_solar[3] = p->chemistry_data.metal_mass_fraction[chemistry_element_N] *
                   cooling->SolarAbundances_inv[3 /* N */];

  ratio_solar[4] = p->chemistry_data.metal_mass_fraction[chemistry_element_O] *
                   cooling->SolarAbundances_inv[4 /* O */];

  ratio_solar[5] = p->chemistry_data.metal_mass_fraction[chemistry_element_Ne] *
                   cooling->SolarAbundances_inv[5 /* Ne */];

  ratio_solar[6] = p->chemistry_data.metal_mass_fraction[chemistry_element_Mg] *
                   cooling->SolarAbundances_inv[6 /* Mg */];

  ratio_solar[7] = p->chemistry_data.metal_mass_fraction[chemistry_element_Si] *
                   cooling->SolarAbundances_inv[7 /* Si */];

  /* For S, we use the same ratio as Si */
  ratio_solar[8] = p->chemistry_data.metal_mass_fraction[chemistry_element_Si] *
                   cooling->SolarAbundances_inv[7 /* Si */] *
                   cooling->S_over_Si_ratio_in_solar;

  /* For Ca, we use the same ratio as Si */
  ratio_solar[9] = p->chemistry_data.metal_mass_fraction[chemistry_element_Si] *
                   cooling->SolarAbundances_inv[7 /* Si */] *
                   cooling->Ca_over_Si_ratio_in_solar;

  ratio_solar[10] =
      p->chemistry_data.metal_mass_fraction[chemistry_element_Fe] *
      cooling->SolarAbundances_inv[10 /* Fe */];
}

/**
 * @brief Computes the extra heat from Helium reionisation at a given redshift.
 *
 * We follow the implementation of Wiersma et al. 2009, MNRAS, 399, 574-600,
 * section. 2. The calculation returns energy in CGS.
 *
 * Note that delta_z is negative.
 *
 * @param z The current redshift.
 * @param delta_z The change in redhsift over the course of this time-step.
 * @param cooling The #cooling_function_data used in the run.
 * @return Helium reionization energy in CGS units.
 */
__attribute__((always_inline)) INLINE double
eagle_helium_reionization_extraheat(
    double z, double delta_z, const struct cooling_function_data *cooling) {

#ifdef SWIFT_DEBUG_CHECKS
  if (delta_z > 0.f) error("Invalid value for delta_z. Should be negative.");
#endif

  /* Recover the values we need */
  const double z_centre = cooling->He_reion_z_centre;
  const double z_sigma = cooling->He_reion_z_sigma;
  const double heat_cgs = cooling->He_reion_heat_cgs;

  double extra_heat = 0.;

  /* Integral of the Gaussian between z and z - delta_z */
  extra_heat += erf((z - delta_z - z_centre) / (M_SQRT2 * z_sigma));
  extra_heat -= erf((z - z_centre) / (M_SQRT2 * z_sigma));

  /* Multiply by the normalisation factor */
  extra_heat *= heat_cgs * 0.5;

  return extra_heat;
}

/**
 * @brief Computes the log_10 of the temperature corresponding to a given
 * internal energy, hydrogen number density, Helium fraction and redshift.
 *
 * Note that the redshift is implicitly passed in via the currently loaded
 * tables in the #cooling_function_data.
 *
 * For the low-z case, we interpolate the flattened 4D table 'u_to_temp' that
 * is arranged in the following way:
 * - 1st dim: redshift, length = eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 3rd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * For the high-z case, we interpolate the flattened 3D table 'u_to_temp' that
 * is arranged in the following way:
 * - 1st dim: Hydrogen density, length = eagle_cooling_N_density
 * - 2nd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 3rd dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * @param log_10_u_cgs Log base 10 of internal energy in cgs.
 * @param redshift Current redshift.
 * @param n_H_index Index along the Hydrogen density dimension.
 * @param He_index Index along the Helium fraction dimension.
 * @param d_n_H Offset between Hydrogen density and table[n_H_index].
 * @param d_He Offset between helium fraction and table[He_index].
 * @param cooling #cooling_function_data structure.
 *
 * @param compute_dT_du Do we want to compute dT/du ?
 * @param dT_du (return) The value of dT/du
 *
 * @return log_10 of the temperature.
 */
__attribute__((always_inline)) INLINE double eagle_convert_u_to_temp(
    const double log_10_u_cgs, const float redshift, const int compute_dT_du,
    float *dT_du, int n_H_index, int He_index, float d_n_H, float d_He,
    const struct cooling_function_data *restrict cooling) {

  /* Get index of u along the internal energy axis */
  int u_index;
  float d_u;
  get_index_1d(cooling->Therm, eagle_cooling_N_temperature, log_10_u_cgs,
               &u_index, &d_u);

  /* Interpolate temperature table to return temperature for current
   * internal energy (use 3D interpolation for high redshift table,
   * otherwise 4D) */
  float log_10_T;
  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    log_10_T = interpolation_3d(cooling->table.temperature,   /* */
                                n_H_index, He_index, u_index, /* */
                                d_n_H, d_He, d_u,             /* */
                                eagle_cooling_N_density,      /* */
                                eagle_cooling_N_He_frac,      /* */
                                eagle_cooling_N_temperature); /* */
  } else {

    log_10_T =
        interpolation_4d(cooling->table.temperature,                  /* */
                         /*z_index=*/0, n_H_index, He_index, u_index, /* */
                         cooling->dz, d_n_H, d_He, d_u,               /* */
                         eagle_cooling_N_loaded_redshifts,            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */
  }

  if (compute_dT_du) {

    float log_10_T_high, log_10_T_low;

    /* Interpolate temperature table to return temperature for internal energy
     * at grid point above current internal energy for computing dT_du used for
     * calculation of dlambda_du in cooling.c (use 3D interpolation for high
     * redshift table, otherwise 4D) */
    if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

      log_10_T_high = interpolation_3d(cooling->table.temperature,   /* */
                                       n_H_index, He_index, u_index, /* */
                                       d_n_H, d_He, /*delta_u=*/1.f, /* */
                                       eagle_cooling_N_density,      /* */
                                       eagle_cooling_N_He_frac,      /* */
                                       eagle_cooling_N_temperature); /* */

    } else {

      log_10_T_high =
          interpolation_4d(cooling->table.temperature,                  /* */
                           /*z_index=*/0, n_H_index, He_index, u_index, /* */
                           cooling->dz, d_n_H, d_He, /*delta_u=*/1.f,   /* */
                           eagle_cooling_N_loaded_redshifts,            /* */
                           eagle_cooling_N_density,                     /* */
                           eagle_cooling_N_He_frac,                     /* */
                           eagle_cooling_N_temperature);                /* */
    }

    /* Interpolate temperature table to return temperature for internal energy
     * at grid point below current internal energy for computing dT_du used for
     * calculation of dlambda_du in cooling.c (use 3D interpolation for high
     * redshift table, otherwise 4D) */
    if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

      log_10_T_low = interpolation_3d(cooling->table.temperature,   /* */
                                      n_H_index, He_index, u_index, /* */
                                      d_n_H, d_He, /*delta_u=*/0.f, /* */
                                      eagle_cooling_N_density,      /* */
                                      eagle_cooling_N_He_frac,      /* */
                                      eagle_cooling_N_temperature); /* */

    } else {

      log_10_T_low =
          interpolation_4d(cooling->table.temperature,                  /* */
                           /*z_index=*/0, n_H_index, He_index, u_index, /* */
                           cooling->dz, d_n_H, d_He, /*delta_u=*/0.f,   /* */
                           eagle_cooling_N_loaded_redshifts,            /* */
                           eagle_cooling_N_density,                     /* */
                           eagle_cooling_N_He_frac,                     /* */
                           eagle_cooling_N_temperature);                /* */
    }

    /* Calculate dT/du */
    const float delta_u = exp(cooling->Therm[u_index + 1] * M_LN10) -
                          exp(cooling->Therm[u_index] * M_LN10);
    *dT_du =
        (exp(M_LN10 * log_10_T_high) - exp(M_LN10 * log_10_T_low)) / delta_u;
  }

  /* Special case for temperatures below the start of the table */
  if (u_index == 0 && d_u == 0.f) {

    /* The temperature is multiplied by u / 10^T[0]
     * where T[0] is the first entry in the table */
    log_10_T += log_10_u_cgs - cooling->Temp[0];
  }

  return log_10_T;
}

/**
 * @brief Compute the Compton cooling rate from the CMB at a given
 * redshift, electron abundance, temperature and Hydrogen density.
 *
 * Uses an analytic formula.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param redshift The current redshift.
 * @param n_H_cgs The Hydrogen number density in CGS units.
 * @param temperature The temperature.
 * @param electron_abundance The electron abundance.
 */
__attribute__((always_inline)) INLINE double eagle_Compton_cooling_rate(
    const struct cooling_function_data *cooling, const double redshift,
    const double n_H_cgs, const double temperature,
    const double electron_abundance) {

  const double zp1 = 1. + redshift;
  const double zp1p2 = zp1 * zp1;
  const double zp1p4 = zp1p2 * zp1p2;

  /* CMB temperature at this redshift */
  const double T_CMB = cooling->T_CMB_0 * zp1;

  /* Compton cooling rate */
  return cooling->compton_rate_cgs * (temperature - T_CMB) * zp1p4 *
         electron_abundance / n_H_cgs;
}

/**
 * @brief Computes the cooling rate corresponding to a given internal energy,
 * hydrogen number density, Helium fraction, redshift and metallicity from
 * all the possible channels.
 *
 * 1) Metal-free cooling:
 * We interpolate the flattened 4D table 'H_and_He_net_heating' that is
 * arranged in the following way:
 * - 1st dim: redshift, length = eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 3rd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * 2) Electron abundance
 * We compute the electron abundance by interpolating the flattened 4d table
 * 'H_and_He_electron_abundance' that is arranged in the following way:
 * - 1st dim: redshift, length = eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 3rd dim: Helium fraction, length = eagle_cooling_N_He_frac
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * 3) Compton cooling is applied via the analytic formula.
 *
 * 4) Solar electron abudance
 * We compute the solar electron abundance by interpolating the flattened 3d
 * table 'solar_electron_abundance' that is arranged in the following way:
 * - 1st dim: redshift, length = eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 3rd dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * 5) Metal-line cooling
 * For each tracked element we interpolate the flattened 4D table
 * 'table_metals_net_heating' that is arrange in the following way:
 * - 1st dim: element, length = eagle_cooling_N_metal
 * - 2nd dim: redshift, length = eagle_cooling_N_loaded_redshifts
 * - 3rd dim: Hydrogen density, length = eagle_cooling_N_density
 * - 4th dim: Internal energy, length = eagle_cooling_N_temperature
 *
 * Note that this is a fake 4D interpolation as we do not interpolate
 * along the 1st dimension. We just do this once per element.
 *
 * Since only the temperature changes when cooling a given particle,
 * the redshift, hydrogen number density and helium fraction indices
 * and offsets passed in.
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
 * @param n_H_index Particle hydrogen number density index
 * @param d_n_H Particle hydrogen number density offset
 * @param He_index Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param cooling Cooling data structure
 *
 * @param dlambda_du (return) Derivative of the cooling rate with respect to u.
 * @param element_lambda (return) Cooling rate from each element
 *
 * @return The cooling rate
 */
INLINE static double eagle_metal_cooling_rate(
    double log10_u_cgs, double redshift, double n_H_cgs,
    const float solar_ratio[chemistry_element_count + 2], int n_H_index,
    float d_n_H, int He_index, float d_He,
    const struct cooling_function_data *restrict cooling, double *dlambda_du,
    double *element_lambda) {

#ifdef TO_BE_DONE
  /* used for calculating dlambda_du */
  double temp_lambda_high = 0, temp_lambda_low = 0;
  double h_plus_he_electron_abundance_high = 0;
  double h_plus_he_electron_abundance_low = 0;
  double solar_electron_abundance_high = 0;
  double solar_electron_abundance_low = 0;
  double elem_cool_low = 0, elem_cool_high = 0;
#endif

  /* We only need dT_du if dLambda_du is non-NULL */
  const int compute_dT_du = (dlambda_du != NULL) ? 1 : 0;

  /* Temperature */
  float dT_du = -1.f;
  const double log_10_T =
      eagle_convert_u_to_temp(log10_u_cgs, redshift, compute_dT_du, &dT_du,
                              n_H_index, He_index, d_n_H, d_He, cooling);

  /* Get index along temperature dimension of the tables */
  int T_index;
  float d_T;
  get_index_1d(cooling->Temp, eagle_cooling_N_temperature, log_10_T, &T_index,
               &d_T);

#ifdef TO_BE_DONE
  /* Difference between entries on the temperature table around u */
  const float delta_T = exp(M_LN10 * cooling->Temp[T_index + 1]) -
                        exp(M_LN10 * cooling->Temp[T_index]);
#endif

  /**********************/
  /* Metal-free cooling */
  /**********************/

  double Lambda_free;

  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    /* If we're using the high redshift tables then we don't interpolate
     * in redshift */
    Lambda_free = interpolation_3d(cooling->table.H_plus_He_heating, /* */
                                   n_H_index, He_index, T_index,     /* */
                                   d_n_H, d_He, d_T,                 /* */
                                   eagle_cooling_N_density,          /* */
                                   eagle_cooling_N_He_frac,          /* */
                                   eagle_cooling_N_temperature);     /* */

#ifdef TO_BE_DONE
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du. Pass in NULL pointer for
     * dlambda_du in order to skip */
    if (dlambda_du != NULL) {
      temp_lambda_high = interpolation_3d(
          cooling->table.H_plus_He_heating, n_H_index, He_index, T_index, d_n_h,
          d_He, 1.f, cooling->N_nH, cooling->N_He, cooling->N_Temp);
      temp_lambda_low = interpolation_3d(
          cooling->table.H_plus_He_heating, n_H_index, He_index, T_index, d_n_h,
          d_He, 0.f, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    }
#endif

  } else {

    /* Using normal tables, have to interpolate in redshift */
    Lambda_free =
        interpolation_4d(cooling->table.H_plus_He_heating,            /* */
                         /*z_index=*/0, n_H_index, He_index, T_index, /* */
                         cooling->dz, d_n_H, d_He, d_T,               /* */
                         eagle_cooling_N_loaded_redshifts,            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */

#ifdef TO_BE_DONE
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      temp_lambda_high =
          interpolation_4d(cooling->table.H_plus_He_heating, 0, n_H_index,
                           He_index, T_index, cooling->dz, d_n_h, d_He, 1.f, 2,
                           cooling->N_nH, cooling->N_He, cooling->N_Temp);
      temp_lambda_low =
          interpolation_4d(cooling->table.H_plus_He_heating, 0, n_H_index,
                           He_index, T_index, cooling->dz, d_n_h, d_He, 0.f, 2,
                           cooling->N_nH, cooling->N_He, cooling->N_Temp);
    }
#endif
  }

#ifdef TO_BE_DONE
  if (dlambda_du != NULL) {
    *dlambda_du += (temp_lambda_high - temp_lambda_low) / delta_T * dT_du;
  }
#endif

  /* If we're testing cooling rate contributions write to array */
  if (element_lambda != NULL) {
    element_lambda[0] = Lambda_free;
  }

  /**********************/
  /* Electron abundance */
  /**********************/

  double H_plus_He_electron_abundance;

  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    H_plus_He_electron_abundance =
        interpolation_3d(cooling->table.H_plus_He_electron_abundance, /* */
                         n_H_index, He_index, T_index,                /* */
                         d_n_H, d_He, d_T,                            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */
#ifdef TO_BE_DONE
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du. Pass in NULL pointer for
     * dlambda_du in order to skip */

    h_plus_he_electron_abundance_high =
        interpolation_3d(cooling->table.H_plus_He_electron_abundance, n_H_index,
                         He_index, T_index, d_n_h, d_He, 1.f, cooling->N_nH,
                         cooling->N_He, cooling->N_Temp);
    h_plus_he_electron_abundance_low =
        interpolation_3d(cooling->table.H_plus_He_electron_abundance, n_H_index,
                         He_index, T_index, d_n_h, d_He, 0.f, cooling->N_nH,
                         cooling->N_He, cooling->N_Temp);

#endif

  } else {

    H_plus_He_electron_abundance =
        interpolation_4d(cooling->table.H_plus_He_electron_abundance, /* */
                         /*z_index=*/0, n_H_index, He_index, T_index, /* */
                         cooling->dz, d_n_H, d_He, d_T,               /* */
                         eagle_cooling_N_loaded_redshifts,            /* */
                         eagle_cooling_N_density,                     /* */
                         eagle_cooling_N_He_frac,                     /* */
                         eagle_cooling_N_temperature);                /* */

#ifdef TO_BE_DONE
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    h_plus_he_electron_abundance_high =
        interpolation_4d(cooling->table.H_plus_He_electron_abundance, 0,
                         n_H_index, He_index, T_index, cooling->dz, d_n_h, d_He,
                         1.f, 2, cooling->N_nH, cooling->N_He, cooling->N_Temp);
    h_plus_he_electron_abundance_low =
        interpolation_4d(cooling->table.H_plus_He_electron_abundance, 0,
                         n_H_index, He_index, T_index, cooling->dz, d_n_h, d_He,
                         0.f, 2, cooling->N_nH, cooling->N_He, cooling->N_Temp);
#endif
  }

  /**********************/
  /* Compton cooling    */
  /**********************/

  double Lambda_Compton = 0.;

  /* Do we need to add the inverse Compton cooling? */
  /* It is *not* stored in the tables before re-ionisation */
  if ((redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) ||
      (redshift > cooling->H_reion_z)) {

    const double T = exp10(log_10_T);

    /* Note the minus sign */
    Lambda_Compton -= eagle_Compton_cooling_rate(cooling, redshift, n_H_cgs, T,
                                                 H_plus_He_electron_abundance);
  }

  /* If we're testing cooling rate contributions write to array */
  if (element_lambda != NULL) {
    element_lambda[1] = Lambda_Compton;
  }

  /*******************************/
  /* Solar electron abundance    */
  /*******************************/

  double solar_electron_abundance;

  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    /* If we're using the high redshift tables then we don't interpolate
     * in redshift */
    solar_electron_abundance =
        interpolation_2d(cooling->table.electron_abundance, /* */
                         n_H_index, T_index,                /* */
                         d_n_H, d_T,                        /* */
                         eagle_cooling_N_density,           /* */
                         eagle_cooling_N_temperature);      /* */

#ifdef TO_BE_DONE
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      solar_electron_abundance_high =
          interpolation_2d(cooling->table.electron_abundance, n_H_index,
                           T_index, d_n_h, 1.f, cooling->N_nH, cooling->N_Temp);
      solar_electron_abundance_low =
          interpolation_2d(cooling->table.electron_abundance, n_H_index,
                           T_index, d_n_h, 0.f, cooling->N_nH, cooling->N_Temp);
    }
#endif

  } else {

    /* Using normal tables, have to interpolate in redshift */
    solar_electron_abundance =
        interpolation_3d(cooling->table.electron_abundance, /* */
                         /*z_index=*/0, n_H_index, T_index, /* */
                         cooling->dz, d_n_H, d_T,           /* */
                         eagle_cooling_N_loaded_redshifts,  /* */
                         eagle_cooling_N_density,           /* */
                         eagle_cooling_N_temperature);      /* */

#ifdef TO_BE_DONE
    /* compute values at temperature gridpoints above and below input
     * temperature for calculation of dlambda_du */
    if (dlambda_du != NULL) {
      solar_electron_abundance_high = interpolation_3d(
          cooling->table.electron_abundance, 0, n_H_index, T_index, cooling->dz,
          d_n_h, 1.f, 2, cooling->N_nH, cooling->N_Temp);
      solar_electron_abundance_low = interpolation_3d(
          cooling->table.electron_abundance, 0, n_H_index, T_index, cooling->dz,
          d_n_h, 0.f, 2, cooling->N_nH, cooling->N_Temp);
    }
#endif
  }

  const double electron_abundance_ratio =
      H_plus_He_electron_abundance / solar_electron_abundance;

  /**********************/
  /* Metal-line cooling */
  /**********************/

  /* for each element the cooling rate is multiplied by the ratio of H, He
   * electron abundance to solar electron abundance then by the ratio of the
   * particle metal abundance to solar metal abundance. */

  double lambda_metal[eagle_cooling_N_metal + 2] = {0.};

  if (redshift > cooling->Redshifts[eagle_cooling_N_redshifts - 1]) {

    /* Loop over the metals (ignore H and He) */
    for (int elem = 2; elem < eagle_cooling_N_metal + 2; elem++) {

      if (solar_ratio[elem] > 0.) {

        /* Note that we do not interpolate along the x-axis
         * (element dimension) */
        lambda_metal[elem] =
            interpolation_3d_no_x(cooling->table.metal_heating,   /* */
                                  elem - 2, n_H_index, T_index,   /* */
                                  /*delta_elem=*/0.f, d_n_H, d_T, /* */
                                  eagle_cooling_N_metal,          /* */
                                  eagle_cooling_N_density,        /* */
                                  eagle_cooling_N_temperature);   /* */

        lambda_metal[elem] *= electron_abundance_ratio;
        lambda_metal[elem] *= solar_ratio[elem];
      }

#ifdef TO_BE_DONE
      /* compute values at temperature gridpoints above and below input
       * temperature for calculation of dlambda_du */
      if (dlambda_du != NULL) {
        elem_cool_high = interpolation_3d_no_x(
            cooling->table.metal_heating, elem, n_H_index, T_index, 0.f, d_n_h,
            1.f, cooling->N_Elements, cooling->N_nH, cooling->N_Temp);

        elem_cool_low = interpolation_3d_no_x(
            cooling->table.metal_heating, elem, n_H_index, T_index, 0.f, d_n_h,
            0.f, cooling->N_nH, cooling->N_Temp, cooling->N_Elements);

        *dlambda_du += (elem_cool_high * h_plus_he_electron_abundance_high /
                            solar_electron_abundance_high -
                        elem_cool_low * h_plus_he_electron_abundance_low /
                            solar_electron_abundance_low) /
                       delta_T * dT_du * solar_ratio[elem + 2];
      }
#endif
    }

  } else {

    /* Loop over the metals (ignore H and He) */
    for (int elem = 2; elem < eagle_cooling_N_metal + 2; elem++) {

      if (solar_ratio[elem] > 0.) {

        /* Note that we do not interpolate along the x-axis
         * (element dimension) */
        lambda_metal[elem] = interpolation_4d_no_x(
            cooling->table.metal_heating,                /* */
            elem - 2, /*z_index=*/0, n_H_index, T_index, /* */
            /*delta_elem=*/0.f, cooling->dz, d_n_H, d_T, /* */
            eagle_cooling_N_metal,                       /* */
            eagle_cooling_N_loaded_redshifts,            /* */
            eagle_cooling_N_density,                     /* */
            eagle_cooling_N_temperature);                /* */

        lambda_metal[elem] *= electron_abundance_ratio;
        lambda_metal[elem] *= solar_ratio[elem];
      }

#ifdef TO_BE_DONE
      /* compute values at temperature gridpoints above and below input
       * temperature for calculation of dlambda_du */
      if (dlambda_du != NULL) {
        elem_cool_high = interpolation_4d_no_x(
            cooling->table.metal_heating, elem, 0, n_H_index, T_index, 0.,
            cooling->dz, d_n_h, 1.f, cooling->N_Elements, 2, cooling->N_nH,
            cooling->N_Temp);

        elem_cool_low = interpolation_4d_no_x(
            cooling->table.metal_heating, elem, 0, n_H_index, T_index, 0.,
            cooling->dz, d_n_h, 0.f, cooling->N_Elements, 2, cooling->N_nH,
            cooling->N_Temp);

        *dlambda_du += (elem_cool_high * h_plus_he_electron_abundance_high /
                            solar_electron_abundance_high -
                        elem_cool_low * h_plus_he_electron_abundance_low /
                            solar_electron_abundance_low) /
                       delta_T * dT_du * solar_ratio[elem + 2];
      }
#endif
    }
  }

  if (element_lambda != NULL) {
    for (int elem = 2; elem < eagle_cooling_N_metal + 2; ++elem) {
      element_lambda[elem] = lambda_metal[elem];
    }
  }

  /* Sum up all the contributions */
  double Lambda_net = Lambda_free + Lambda_Compton;
  for (int elem = 2; elem < eagle_cooling_N_metal + 2; ++elem) {
    Lambda_net += lambda_metal[elem];
  }

  return Lambda_net;
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
 * @param n_H_index Particle hydrogen number density index
 * @param d_n_H Particle hydrogen number density offset
 * @param He_index Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param cooling #cooling_function_data structure
 *
 * @param dLambdaNet_du (return) Derivative of the cooling rate with respect to
 * u.
 *
 * @return The cooling rate
 */
INLINE static double eagle_cooling_rate(
    double log_u_cgs, double redshift, double n_H_cgs,
    const float abundance_ratio[chemistry_element_count + 2], int n_H_index,
    float d_n_H, int He_index, float d_He,
    const struct cooling_function_data *restrict cooling,
    double *dLambdaNet_du) {

  return eagle_metal_cooling_rate(log_u_cgs / M_LN10, redshift, n_H_cgs,
                                  abundance_ratio, n_H_index, d_n_H, He_index,
                                  d_He, cooling, dLambdaNet_du,
                                  /*element_lambda=*/NULL);
}

#endif /* SWIFT_EAGLE_COOLING_RATES_H */
