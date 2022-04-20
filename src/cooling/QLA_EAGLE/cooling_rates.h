/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_QLA_EAGLE_COOLING_RATES_H
#define SWIFT_QLA_EAGLE_COOLING_RATES_H

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
 * In this Quick Lyman-alpha model, we only care about H and He,
 * which are taken from the code's constants.
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
__attribute__((always_inline)) INLINE static void abundance_ratio_to_solar(
    const struct part *p, const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    float ratio_solar[qla_eagle_cooling_N_abundances]) {

  const float XH = 1. - phys_const->const_primordial_He_fraction;
  const float XHe = phys_const->const_primordial_He_fraction;

  ratio_solar[0] = XH * cooling->SolarAbundances_inv[0 /* H */];
  ratio_solar[1] = XHe * cooling->SolarAbundances_inv[1 /* He */];

  /* No metals... */
  ratio_solar[2] = 0.f;
  ratio_solar[3] = 0.f;
  ratio_solar[4] = 0.f;
  ratio_solar[5] = 0.f;
  ratio_solar[6] = 0.f;
  ratio_solar[7] = 0.f;
  ratio_solar[8] = 0.f;
  ratio_solar[9] = 0.f;
  ratio_solar[10] = 0.f;
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
qla_eagle_helium_reionization_extraheat(
    const double z, const double delta_z,
    const struct cooling_function_data *cooling) {

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
 * - 1st dim: redshift, length = qla_eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = qla_eagle_cooling_N_density
 * - 3rd dim: Helium fraction, length = qla_eagle_cooling_N_He_frac
 * - 4th dim: Internal energy, length = qla_eagle_cooling_N_temperature
 *
 * For the high-z case, we interpolate the flattened 3D table 'u_to_temp' that
 * is arranged in the following way:
 * - 1st dim: Hydrogen density, length = qla_eagle_cooling_N_density
 * - 2nd dim: Helium fraction, length = qla_eagle_cooling_N_He_frac
 * - 3rd dim: Internal energy, length = qla_eagle_cooling_N_temperature
 *
 * @param log_10_u_cgs Log base 10 of internal energy in cgs.
 * @param redshift Current redshift.
 * @param n_H_index Index along the Hydrogen density dimension.
 * @param He_index Index along the Helium fraction dimension.
 * @param d_n_H Offset between Hydrogen density and table[n_H_index].
 * @param d_He Offset between helium fraction and table[He_index].
 * @param cooling #cooling_function_data structure.
 *
 * @return log_10 of the temperature.
 */
__attribute__((always_inline)) INLINE double qla_eagle_convert_u_to_temp(
    const double log_10_u_cgs, const float redshift, const int n_H_index,
    const int He_index, const float d_n_H, const float d_He,
    const struct cooling_function_data *cooling) {

  /* Get index of u along the internal energy axis */
  int u_index;
  float d_u;
  get_index_1d(cooling->Therm, qla_eagle_cooling_N_temperature, log_10_u_cgs,
               &u_index, &d_u);

  /* Interpolate temperature table to return temperature for current
   * internal energy (use 3D interpolation for high redshift table,
   * otherwise 4D) */
  float log_10_T;
  if (redshift > cooling->Redshifts[qla_eagle_cooling_N_redshifts - 1]) {

    log_10_T = interpolation_3d(cooling->table.temperature,       /* */
                                n_H_index, He_index, u_index,     /* */
                                d_n_H, d_He, d_u,                 /* */
                                qla_eagle_cooling_N_density,      /* */
                                qla_eagle_cooling_N_He_frac,      /* */
                                qla_eagle_cooling_N_temperature); /* */
  } else {

    log_10_T =
        interpolation_4d(cooling->table.temperature,                  /* */
                         /*z_index=*/0, n_H_index, He_index, u_index, /* */
                         cooling->dz, d_n_H, d_He, d_u,               /* */
                         qla_eagle_cooling_N_loaded_redshifts,        /* */
                         qla_eagle_cooling_N_density,                 /* */
                         qla_eagle_cooling_N_He_frac,                 /* */
                         qla_eagle_cooling_N_temperature);            /* */
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
__attribute__((always_inline)) INLINE double qla_eagle_Compton_cooling_rate(
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
 * - 1st dim: redshift, length = qla_eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = qla_eagle_cooling_N_density
 * - 3rd dim: Helium fraction, length = qla_eagle_cooling_N_He_frac
 * - 4th dim: Internal energy, length = qla_eagle_cooling_N_temperature
 *
 * 2) Electron abundance
 * We compute the electron abundance by interpolating the flattened 4d table
 * 'H_and_He_electron_abundance' that is arranged in the following way:
 * - 1st dim: redshift, length = qla_eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = qla_eagle_cooling_N_density
 * - 3rd dim: Helium fraction, length = qla_eagle_cooling_N_He_frac
 * - 4th dim: Internal energy, length = qla_eagle_cooling_N_temperature
 *
 * 3) Compton cooling is applied via the analytic formula.
 *
 * 4) Solar electron abudance
 * We compute the solar electron abundance by interpolating the flattened 3d
 * table 'solar_electron_abundance' that is arranged in the following way:
 * - 1st dim: redshift, length = qla_eagle_cooling_N_loaded_redshifts
 * - 2nd dim: Hydrogen density, length = qla_eagle_cooling_N_density
 * - 3rd dim: Internal energy, length = qla_eagle_cooling_N_temperature
 *
 * 5) Metal-line cooling
 * For each tracked element we interpolate the flattened 4D table
 * 'table_metals_net_heating' that is arrange in the following way:
 * - 1st dim: element, length = qla_eagle_cooling_N_metal
 * - 2nd dim: redshift, length = qla_eagle_cooling_N_loaded_redshifts
 * - 3rd dim: Hydrogen density, length = qla_eagle_cooling_N_density
 * - 4th dim: Internal energy, length = qla_eagle_cooling_N_temperature
 *
 * Note that this is a fake 4D interpolation as we do not interpolate
 * along the 1st dimension. We just do this once per element.
 *
 * Since only the temperature changes when cooling a given particle,
 * the redshift, hydrogen number density and helium fraction indices
 * and offsets passed in.
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
 * @param element_lambda (return) Cooling rate from each element
 *
 * @return The cooling rate
 */
INLINE static double qla_eagle_metal_cooling_rate(
    const double log10_u_cgs, const double redshift, const double n_H_cgs,
    const float solar_ratio[qla_eagle_cooling_N_abundances],
    const int n_H_index, const float d_n_H, const int He_index,
    const float d_He, const struct cooling_function_data *cooling,
    double *element_lambda) {

  /* Temperature */
  const double log_10_T = qla_eagle_convert_u_to_temp(
      log10_u_cgs, redshift, n_H_index, He_index, d_n_H, d_He, cooling);

  /* Get index along temperature dimension of the tables */
  int T_index;
  float d_T;
  get_index_1d(cooling->Temp, qla_eagle_cooling_N_temperature, log_10_T,
               &T_index, &d_T);

  /**********************/
  /* Metal-free cooling */
  /**********************/

  double Lambda_free;

  if (redshift > cooling->Redshifts[qla_eagle_cooling_N_redshifts - 1]) {

    /* If we're using the high redshift tables then we don't interpolate
     * in redshift */
    Lambda_free = interpolation_3d(cooling->table.H_plus_He_heating, /* */
                                   n_H_index, He_index, T_index,     /* */
                                   d_n_H, d_He, d_T,                 /* */
                                   qla_eagle_cooling_N_density,      /* */
                                   qla_eagle_cooling_N_He_frac,      /* */
                                   qla_eagle_cooling_N_temperature); /* */
  } else {

    /* Using normal tables, have to interpolate in redshift */
    Lambda_free =
        interpolation_4d(cooling->table.H_plus_He_heating,            /* */
                         /*z_index=*/0, n_H_index, He_index, T_index, /* */
                         cooling->dz, d_n_H, d_He, d_T,               /* */
                         qla_eagle_cooling_N_loaded_redshifts,        /* */
                         qla_eagle_cooling_N_density,                 /* */
                         qla_eagle_cooling_N_He_frac,                 /* */
                         qla_eagle_cooling_N_temperature);            /* */
  }

  /* If we're testing cooling rate contributions write to array */
  if (element_lambda != NULL) {
    element_lambda[0] = Lambda_free;
  }

  /**********************/
  /* Electron abundance */
  /**********************/

  double H_plus_He_electron_abundance;

  if (redshift > cooling->Redshifts[qla_eagle_cooling_N_redshifts - 1]) {

    H_plus_He_electron_abundance =
        interpolation_3d(cooling->table.H_plus_He_electron_abundance, /* */
                         n_H_index, He_index, T_index,                /* */
                         d_n_H, d_He, d_T,                            /* */
                         qla_eagle_cooling_N_density,                 /* */
                         qla_eagle_cooling_N_He_frac,                 /* */
                         qla_eagle_cooling_N_temperature);            /* */

  } else {

    H_plus_He_electron_abundance =
        interpolation_4d(cooling->table.H_plus_He_electron_abundance, /* */
                         /*z_index=*/0, n_H_index, He_index, T_index, /* */
                         cooling->dz, d_n_H, d_He, d_T,               /* */
                         qla_eagle_cooling_N_loaded_redshifts,        /* */
                         qla_eagle_cooling_N_density,                 /* */
                         qla_eagle_cooling_N_He_frac,                 /* */
                         qla_eagle_cooling_N_temperature);            /* */
  }

  /**********************/
  /* Compton cooling    */
  /**********************/

  double Lambda_Compton = 0.;

  /* Do we need to add the inverse Compton cooling? */
  /* It is *not* stored in the tables before re-ionisation */
  if ((redshift > cooling->Redshifts[qla_eagle_cooling_N_redshifts - 1]) ||
      (redshift > cooling->H_reion_z)) {

    const double T = exp10(log_10_T);

    /* Note the minus sign */
    Lambda_Compton -= qla_eagle_Compton_cooling_rate(
        cooling, redshift, n_H_cgs, T, H_plus_He_electron_abundance);
  }

  /* If we're testing cooling rate contributions write to array */
  if (element_lambda != NULL) {
    element_lambda[1] = Lambda_Compton;
  }

  /*******************************/
  /* Solar electron abundance    */
  /*******************************/

  double solar_electron_abundance;

  if (redshift > cooling->Redshifts[qla_eagle_cooling_N_redshifts - 1]) {

    /* If we're using the high redshift tables then we don't interpolate
     * in redshift */
    solar_electron_abundance =
        interpolation_2d(cooling->table.electron_abundance, /* */
                         n_H_index, T_index,                /* */
                         d_n_H, d_T,                        /* */
                         qla_eagle_cooling_N_density,       /* */
                         qla_eagle_cooling_N_temperature);  /* */

  } else {

    /* Using normal tables, have to interpolate in redshift */
    solar_electron_abundance =
        interpolation_3d(cooling->table.electron_abundance,    /* */
                         /*z_index=*/0, n_H_index, T_index,    /* */
                         cooling->dz, d_n_H, d_T,              /* */
                         qla_eagle_cooling_N_loaded_redshifts, /* */
                         qla_eagle_cooling_N_density,          /* */
                         qla_eagle_cooling_N_temperature);     /* */
  }

  const double electron_abundance_ratio =
      H_plus_He_electron_abundance / solar_electron_abundance;

  /**********************/
  /* Metal-line cooling */
  /**********************/

  /* for each element the cooling rate is multiplied by the ratio of H, He
   * electron abundance to solar electron abundance then by the ratio of the
   * particle metal abundance to solar metal abundance. */

  double lambda_metal[qla_eagle_cooling_N_metal + 2] = {0.};

  if (redshift > cooling->Redshifts[qla_eagle_cooling_N_redshifts - 1]) {

    /* Loop over the metals (ignore H and He) */
    for (int elem = 2; elem < qla_eagle_cooling_N_metal + 2; elem++) {

      if (solar_ratio[elem] > 0.) {

        /* Note that we do not interpolate along the x-axis
         * (element dimension) */
        lambda_metal[elem] =
            interpolation_3d_no_x(cooling->table.metal_heating,     /* */
                                  elem - 2, n_H_index, T_index,     /* */
                                  /*delta_elem=*/0.f, d_n_H, d_T,   /* */
                                  qla_eagle_cooling_N_metal,        /* */
                                  qla_eagle_cooling_N_density,      /* */
                                  qla_eagle_cooling_N_temperature); /* */

        lambda_metal[elem] *= electron_abundance_ratio;
        lambda_metal[elem] *= solar_ratio[elem];
      }
    }

  } else {

    /* Loop over the metals (ignore H and He) */
    for (int elem = 2; elem < qla_eagle_cooling_N_metal + 2; elem++) {

      if (solar_ratio[elem] > 0.) {

        /* Note that we do not interpolate along the x-axis
         * (element dimension) */
        lambda_metal[elem] = interpolation_4d_no_x(
            cooling->table.metal_heating,                /* */
            elem - 2, /*z_index=*/0, n_H_index, T_index, /* */
            /*delta_elem=*/0.f, cooling->dz, d_n_H, d_T, /* */
            qla_eagle_cooling_N_metal,                   /* */
            qla_eagle_cooling_N_loaded_redshifts,        /* */
            qla_eagle_cooling_N_density,                 /* */
            qla_eagle_cooling_N_temperature);            /* */

        lambda_metal[elem] *= electron_abundance_ratio;
        lambda_metal[elem] *= solar_ratio[elem];

        // message("lambda[%d]=%e", elem, lambda_metal[elem]);
      }
    }
  }

  if (element_lambda != NULL) {
    for (int elem = 2; elem < qla_eagle_cooling_N_metal + 2; ++elem) {
      element_lambda[elem] = lambda_metal[elem];
    }
  }

  /* Sum up all the contributions */
  double Lambda_net = Lambda_free + Lambda_Compton;
  for (int elem = 2; elem < qla_eagle_cooling_N_metal + 2; ++elem) {
    Lambda_net += lambda_metal[elem];
  }

  return Lambda_net;
}

/**
 * @brief Wrapper function used to calculate cooling rate.
 * Table indices and offsets for redshift, hydrogen number density and
 * helium fraction are passed it so as to compute them only once per particle.
 *
 * @param log10_u_cgs Log base 10 of internal energy per unit mass in CGS units.
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
 * @return The cooling rate
 */
INLINE static double qla_eagle_cooling_rate(
    const double log10_u_cgs, const double redshift, const double n_H_cgs,
    const float abundance_ratio[qla_eagle_cooling_N_abundances],
    const int n_H_index, const float d_n_H, const int He_index,
    const float d_He, const struct cooling_function_data *cooling) {

  return qla_eagle_metal_cooling_rate(
      log10_u_cgs, redshift, n_H_cgs, abundance_ratio, n_H_index, d_n_H,
      He_index, d_He, cooling, /* element_lambda=*/NULL);
}

#endif /* SWIFT_QLA_EAGLE_COOLING_RATES_H */
