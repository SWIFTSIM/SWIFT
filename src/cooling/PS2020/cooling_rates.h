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
#ifndef SWIFT_PS2020_COOLING_RATES_H
#define SWIFT_PS2020_COOLING_RATES_H

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "cooling_tables.h"
#include "error.h"
#include "exp10.h"
#include "interpolate.h"

/**
 * @brief Returns the element index of the particle-carried chemistry field
 * corresponding to a given cooling element.
 */
__attribute__((always_inline)) INLINE int element_from_table_to_code(
    const enum colibre_cooling_element i) {

#ifdef SWIFT_DEBUG_CHECKS
  if ((i >= colibre_cooling_N_elementtypes) || (i < 0))
    error("Outside range of elements in cooling tables");
#endif

  switch (i) {
    case element_H:
      return chemistry_element_H;
    case element_He:
      return chemistry_element_He;
    case element_C:
      return chemistry_element_C;
    case element_N:
      return chemistry_element_N;
    case element_O:
      return chemistry_element_O;
    case element_Ne:
      return chemistry_element_Ne;
    case element_Mg:
      return chemistry_element_Mg;
    case element_Si:
      return chemistry_element_Si;
      /* S and Ca are not tracked individually; their abundance is
       * assumed to be the same as Si (with respect to solar) */
    case element_S:
      return chemistry_element_Si;
    case element_Ca:
      return chemistry_element_Si;
    case element_Fe:
      return chemistry_element_Fe;
      /* other elements, if used, scale with metallicity */
    case element_OA:
      return -1;
  }

  return -1;
}

/**
 * @brief Compute ratio of mass fraction to solar mass fraction
 * for each element carried by a given particle.
 *
 * The solar abundances are taken from the tables themselves.
 *
 * The PS2020 chemistry model does not track S and Ca. We assume
 * that their abundance with respect to solar is the same as
 * the ratio for Si.
 *
 * The other un-tracked elements are scaled with the total metallicity.
 *
 * We optionally apply a correction if the user asked for a different
 * ratio.
 *
 * We also re-order the elements such that they match the order of the
 * tables. This is [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe, OA].
 *
 * The solar abundances table (from the cooling struct) is arranged as
 * [H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe].
 *
 * @param p Pointer to #part struct.
 * @param cooling #cooling_function_data struct.
 * @param ratio_solar (return) Array of ratios to solar abundances.
 *
 * @return The log10 of the total metallicity with respect to solar, i.e.
 * log10(Z / Z_sun).
 */
__attribute__((always_inline)) INLINE static float abundance_ratio_to_solar(
    const struct part *p, const struct cooling_function_data *cooling,
    float ratio_solar[colibre_cooling_N_elementtypes]) {

  /* Get the particle's metal mass fractions (M_x / M) */
  const float *Z_mass_frac = chemistry_get_metal_mass_fraction_for_cooling(p);

  /* Convert mass fractions to abundances (nx/nH) and compute metal mass */
  for (int elem_nr = element_H; elem_nr < element_OA; elem_nr++) {

    enum colibre_cooling_element elem = (enum colibre_cooling_element)elem_nr;

    /* Normal elements: Get the abundance from the particle carried arrays */
    if ((elem != element_S) && (elem != element_Ca)) {

      const int indx1d = row_major_index_2d(cooling->indxZsol, elem,
                                            colibre_cooling_N_metallicity,
                                            colibre_cooling_N_elementtypes);

      const float Mfrac = Z_mass_frac[element_from_table_to_code(elem)];

      /* ratio_X = ((M_x/M) / (M_H/M)) * (m_H / m_X) * (1 / Z_sun_X) */
      ratio_solar[elem] =
          (Mfrac / Z_mass_frac[element_from_table_to_code(element_H)]) *
          cooling->atomicmass[element_H] * cooling->atomicmass_inv[elem] *
          cooling->Abundances_inv[indx1d];

      /* Special case: S scales with Si */
    } else if (elem == element_S) {

      ratio_solar[element_S] =
          ratio_solar[element_Si] * cooling->S_over_Si_ratio_in_solar;

      /* Special case: Ca scales with Si */
    } else if (elem == element_Ca) {

      ratio_solar[element_Ca] =
          ratio_solar[element_Si] * cooling->Ca_over_Si_ratio_in_solar;

    } else {
#ifdef SWIFT_DEBUG_CHECKS
      error("Invalid element!");
#endif
    }
  }

  /* Get total metallicity [Z/Z_sun] from the particle data */
  const float Z_total =
      (float)chemistry_get_total_metal_mass_fraction_for_star_formation(p);
  float ZZsol = Z_total * cooling->Zsol_inv[0];
  if (ZZsol <= 0.0f) ZZsol = FLT_MIN;
  const float logZZsol = log10f(ZZsol);

  /* All other elements (element_OA): scale with metallicity */
  ratio_solar[element_OA] = ZZsol;

  /* Get index and offset from the metallicity table conresponding to this
   * metallicity */
  int met_index;
  float d_met;

  get_index_1d(cooling->Metallicity, colibre_cooling_N_metallicity, logZZsol,
               &met_index, &d_met);

  /* At this point ratio_solar is (nx/nH) / (nx/nH)_sol.
   * To multiply with the tables, we want the individual abundance ratio
   * relative to what is used in the tables for each metallicity
   */

  /* For example: for a metallicity of 1 per cent solar, the metallicity bin
   * for logZZsol = -2 has already the reduced cooling rates for each element;
   * it should therefore NOT be multiplied by 0.01 again.
   *
   * BUT: if e.g. Carbon is twice as abundant as the solar abundance ratio,
   * i.e. nC / nH = 0.02 * (nC/nH)_sol for the overall metallicity of 0.01,
   * the Carbon cooling rate is multiplied by 2
   *
   * We only do this if we are not in the primodial metallicity bin in which
   * case the ratio to solar should be 0.
   */

  for (int i = 0; i < colibre_cooling_N_elementtypes; i++) {

    /* Are we considering a metal and are *not* in the primodial metallicity
     * bin? Or are we looking at H or He? */
    if ((met_index > 0) || (i == element_H) || (i == element_He)) {

      /* Get min/max abundances */
      const float log_nx_nH_min = cooling->LogAbundances[row_major_index_2d(
          met_index, i, colibre_cooling_N_metallicity,
          colibre_cooling_N_elementtypes)];

      const float log_nx_nH_max = cooling->LogAbundances[row_major_index_2d(
          met_index + 1, i, colibre_cooling_N_metallicity,
          colibre_cooling_N_elementtypes)];

      /* Get solar abundances */
      const float log_nx_nH_sol = cooling->LogAbundances[row_major_index_2d(
          cooling->indxZsol, i, colibre_cooling_N_metallicity,
          colibre_cooling_N_elementtypes)];

      /* Interpolate ! (linearly in log-space) */
      const float log_nx_nH =
          (log_nx_nH_min * (1.f - d_met) + log_nx_nH_max * d_met) -
          log_nx_nH_sol;

      ratio_solar[i] *= exp10f(-log_nx_nH);

    } else {

      /* Primordial bin --> Z/Z_sun is 0 for that element */
      ratio_solar[i] = 0.;
    }
  }

  /* at this point ratio_solar is (nx/nH) / (nx/nH)_table [Z],
   * the metallicity dependent abundance ratio for solar abundances.
   * We also return the total metallicity */

  return logZZsol;
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
__attribute__((always_inline)) INLINE static double
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
 * internal energy, hydrogen number density, metallicity and redshift
 *
 * @param log_10_u_cgs Log base 10 of internal energy in cgs.
 * @param redshift Current redshift.
 * @param n_H_index Index along the Hydrogen density dimension.
 * @param d_n_H Offset between Hydrogen density and table[n_H_index].
 * @param met_index Index along the metallicity dimension.
 * @param d_met Offset between metallicity and table[met_index].
 * @param red_index Index along the redshift dimension.
 * @param d_red Offset between redshift and table[red_index].
 * @param cooling #cooling_function_data structure.
 *
 * @return log_10 of the temperature.
 *
 * TO DO: outside table ranges, it uses at the moment the minimum, maximu value
 */
__attribute__((always_inline)) INLINE static float colibre_convert_u_to_temp(
    const double log_10_u_cgs, const float redshift, const int n_H_index,
    const float d_n_H, const int met_index, const float d_met,
    const int red_index, const float d_red,
    const struct cooling_function_data *cooling) {

  /* Get index of u along the internal energy axis */
  int u_index;
  float d_u;

  get_index_1d(cooling->Therm, colibre_cooling_N_internalenergy, log_10_u_cgs,
               &u_index, &d_u);

  /* Interpolate temperature table to return temperature for current
   * internal energy */
  float log_10_T;

  /* Temperature from internal energy */
  log_10_T = interpolation_4d(
      cooling->table.T_from_U, red_index, u_index, met_index, n_H_index, d_red,
      d_u, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density);

  /* General interpolation returns u for the first (or last) element
   * but we want to extrapolate in this case. We assume that the u-T relation
   * does not change outside the table range */
  if (u_index == 0 && d_u == 0.f) {

    log_10_T += log_10_u_cgs - cooling->Therm[0];

  } else if (u_index >= colibre_cooling_N_internalenergy - 2 && d_u == 1.f) {

    log_10_T +=
        log_10_u_cgs - cooling->Therm[colibre_cooling_N_internalenergy - 1];
  }

  return log_10_T;
}

/**
 * @brief Computes the log_10 of the internal energy corresponding to a given
 * temperature, hydrogen number density, metallicity and redshift
 *
 * @param log_10_T Log base 10 of temperature in K
 * @param redshift Current redshift.
 * @param n_H_index Index along the Hydrogen density dimension.
 * @param d_n_H Offset between Hydrogen density and table[n_H_index].
 * @param met_index Index along the metallicity dimension.
 * @param d_met Offset between metallicity and table[met_index].
 * @param red_index Index along the redshift dimension.
 * @param d_red Offset between redshift and table[red_index].
 * @param cooling #cooling_function_data structure.
 *
 * @return log_10 of the internal energy in cgs
 *
 * TO DO: outside table ranges, it uses at the moment the minimum, maximu value
 */
__attribute__((always_inline)) INLINE static float colibre_convert_temp_to_u(
    const double log_10_T, const float redshift, const int n_H_index,
    const float d_n_H, const int met_index, const float d_met,
    const int red_index, const float d_red,
    const struct cooling_function_data *cooling) {

  /* Get index of u along the internal energy axis */
  int T_index;
  float d_T;

  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log_10_T, &T_index,
               &d_T);

  /* Interpolate internal energy table to return internal energy for current
   * temperature */
  float log_10_U;

  /* Internal energy from temperature*/
  log_10_U = interpolation_4d(
      cooling->table.U_from_T, red_index, T_index, met_index, n_H_index, d_red,
      d_T, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_temperature, colibre_cooling_N_metallicity,
      colibre_cooling_N_density);

  /* General interpolation returns u for the first (or last) element
   * but we want to extrapolate in this case. We assume that the u-T relation
   * does not change outside the table range */
  if (T_index == 0 && d_T == 0.f) {

    log_10_U += cooling->Temp[0] - log_10_T;
  } else if (T_index >= colibre_cooling_N_temperature - 2 && d_T == 1.f) {

    log_10_U += cooling->Temp[colibre_cooling_N_temperature - 1] - log_10_T;
  }

  return log_10_U;
}

/**
 * @brief Computes the mean particle mass for a given
 * metallicity, temperature, redshift, and density.
 *
 * @param log_T_cgs Log base 10 of temperature in K
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param d_n_H Offset between Hydrogen density and table[n_H_index]
 * @param met_index Index along the metallicity dimension
 * @param d_met Offset between metallicity and table[met_index]
 * @param red_index Index along redshift dimension
 * @param d_red Offset between redshift and table[red_index]
 * @param cooling #cooling_function_data structure
 *
 * @return mean particle mass
 */
INLINE static float colibre_meanparticlemass_temperature(
    const double log_T_cgs, const double redshift, const double n_H_cgs,
    const int n_H_index, const float d_n_H, const int met_index,
    const float d_met, const int red_index, const float d_red,
    const struct cooling_function_data *cooling) {

  /* Get index of T along the temperature axis */
  int T_index;
  float d_T;

  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log_T_cgs,
               &T_index, &d_T);

  const float mu = interpolation_4d(
      cooling->table.Tmu, red_index, T_index, met_index, n_H_index, d_red, d_T,
      d_met, d_n_H, colibre_cooling_N_redshifts, colibre_cooling_N_temperature,
      colibre_cooling_N_metallicity, colibre_cooling_N_density);

  return mu;
}

/**
 * @brief Computes the electron density for a given element
 * abundance ratio, internal energy, redshift, and density.
 *
 * @param log_u_cgs Log base 10 of internal energy in cgs [erg g-1]
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param abundance_ratio Abundance ratio for each element x relative to solar
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param d_n_H Offset between Hydrogen density and table[n_H_index]
 * @param met_index Index along the metallicity dimension
 * @param d_met Offset between metallicity and table[met_index]
 * @param red_index Index along redshift dimension
 * @param d_red Offset between redshift and table[red_index]
 * @param cooling #cooling_function_data structure
 *
 * @return linear electron density in cm-3 (NOT the electron fraction)
 */
INLINE static float colibre_electron_density(
    const double log_u_cgs, const double redshift, const double n_H_cgs,
    const float abundance_ratio[colibre_cooling_N_elementtypes],
    const int n_H_index, const float d_n_H, const int met_index,
    const float d_met, const int red_index, const float d_red,
    const struct cooling_function_data *cooling) {

  /* Get index of u along the internal energy axis */
  int U_index;
  float d_U;

  get_index_1d(cooling->Therm, colibre_cooling_N_internalenergy, log_u_cgs,
               &U_index, &d_U);

  /* n_e / n_H */
  const float electron_fraction = interpolation4d_plus_summation(
      cooling->table.Uelectron_fraction, abundance_ratio, element_H,
      colibre_cooling_N_electrontypes - 4, red_index, U_index, met_index,
      n_H_index, d_red, d_U, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density, colibre_cooling_N_electrontypes);

  return electron_fraction * n_H_cgs;
}

/**
 * @brief Computes the electron density for a given element
 * abundance ratio, temperature, redshift, and density.
 *
 * @param log_T_cgs Log base 10 of temperature
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param ZZsol Metallicity relative to the solar value from the tables
 * @param abundance_ratio Abundance ratio for each element x relative to solar
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param d_n_H Offset between Hydrogen density and table[n_H_index]
 * @param met_index Index along the metallicity dimension
 * @param d_met Offset between metallicity and table[met_index]
 * @param red_index Index along redshift dimension
 * @param d_red Offset between redshift and table[red_index]
 * @param cooling #cooling_function_data structure
 *
 * @return linear electron density in cm-3 (NOT the electron fraction)
 */
INLINE static float colibre_electron_density_temperature(
    const double log_T_cgs, const double redshift, const double n_H_cgs,
    const float ZZsol,
    const float abundance_ratio[colibre_cooling_N_elementtypes],
    const int n_H_index, const float d_n_H, const int met_index,
    const float d_met, const int red_index, const float d_red,
    const struct cooling_function_data *cooling) {

  /* Get index of u along the internal energy axis */
  int T_index;
  float d_T;

  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log_T_cgs,
               &T_index, &d_T);

  /* n_e / n_H */
  const float electron_fraction = interpolation4d_plus_summation(
      cooling->table.Telectron_fraction, abundance_ratio, element_H,
      colibre_cooling_N_electrontypes - 4, red_index, T_index, met_index,
      n_H_index, d_red, d_T, d_met, d_n_H, colibre_cooling_N_redshifts,
      colibre_cooling_N_internalenergy, colibre_cooling_N_metallicity,
      colibre_cooling_N_density, colibre_cooling_N_electrontypes);

  return electron_fraction * n_H_cgs;
}

/**
 * @brief Computes the net cooling rate (heating - cooling) for a given element
 * abundance ratio, internal energy, redshift, and density. The unit of the net
 * cooling rate is Lambda / nH**2 [erg cm^3 s-1] and all input values are in
 * cgs. The Compton cooling is not taken from the tables but calculated
 * analytically and added separately
 *
 * @param log_u_cgs Log base 10 of internal energy in cgs [erg g-1]
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param abundance_ratio Abundance ratio for each element x relative to solar
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param d_n_H Offset between Hydrogen density and table[n_H_index]
 * @param met_index Index along the metallicity dimension
 * @param d_met Offset between metallicity and table[met_index]
 * @param red_index Index along redshift dimension
 * @param d_red Offset between redshift and table[red_index]
 * @param cooling #cooling_function_data structure
 *
 * @param onlyicool if true / 1 only plot cooling channel icool
 * @param onlyiheat if true / 1 only plot cooling channel iheat
 * @param icool cooling channel to be used
 * @param iheat heating channel to be used
 *
 * Throughout the code: onlyicool = onlyiheat = icool = iheat = 0
 * These are only used for testing: examples/CoolingRates/CoolingRatesPS2020
 */
INLINE static double colibre_cooling_rate(
    const double log_u_cgs, const double redshift, const double n_H_cgs,
    const float abundance_ratio[colibre_cooling_N_elementtypes],
    const int n_H_index, const float d_n_H, const int met_index,
    const float d_met, const int red_index, const float d_red,
    const struct cooling_function_data *cooling, const int onlyicool,
    const int onlyiheat, const int icool, const int iheat) {

  /* Set weights for cooling rates */
  float weights_cooling[colibre_cooling_N_cooltypes - 2];
  for (int i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {

    if (i < colibre_cooling_N_elementtypes) {
      /* Use abundance ratios */
      weights_cooling[i] = abundance_ratio[i];
    } else if (i == cooltype_Compton) {
      /* added analytically later, do not use value from table*/
      weights_cooling[i] = 0.f;
    } else {
      /* use same abundances as in the tables */
      weights_cooling[i] = 1.f;
    }
  }

  /* If we care only about one channel, cancel all the other ones */
  if (onlyicool != 0) {
    for (int i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {
      if (i != icool) weights_cooling[i] = 0.f;
    }
  }

  /* Set weights for heating rates */
  float weights_heating[colibre_cooling_N_heattypes - 2];
  for (int i = 0; i < colibre_cooling_N_heattypes - 2; i++) {
    if (i < colibre_cooling_N_elementtypes) {
      weights_heating[i] = abundance_ratio[i];
    } else {
      weights_heating[i] = 1.f; /* use same abundances as in the tables */
    }
  }

  /* If we care only about one channel, cancel all the other ones */
  if (onlyiheat != 0) {
    for (int i = 0; i < colibre_cooling_N_heattypes - 2; i++) {
      if (i != iheat) weights_heating[i] = 0.f;
    }
  }

  /* Get index of u along the internal energy axis */
  int U_index;
  float d_U;
  get_index_1d(cooling->Therm, colibre_cooling_N_internalenergy, log_u_cgs,
               &U_index, &d_U);

  /* n_e / n_H */
  const double electron_fraction = interpolation4d_plus_summation(
      cooling->table.Uelectron_fraction, abundance_ratio, /* */
      element_H, colibre_cooling_N_electrontypes - 4,     /* */
      red_index, U_index, met_index, n_H_index,           /* */
      d_red, d_U, d_met, d_n_H,                           /* */
      colibre_cooling_N_redshifts,                        /* */
      colibre_cooling_N_internalenergy,                   /* */
      colibre_cooling_N_metallicity,                      /* */
      colibre_cooling_N_density,                          /* */
      colibre_cooling_N_electrontypes);                   /* */

  /* Lambda / n_H**2 */
  const double cooling_rate = interpolation4d_plus_summation(
      cooling->table.Ucooling, weights_cooling,   /* */
      element_H, colibre_cooling_N_cooltypes - 3, /* */
      red_index, U_index, met_index, n_H_index,   /* */
      d_red, d_U, d_met, d_n_H,                   /* */
      colibre_cooling_N_redshifts,                /* */
      colibre_cooling_N_internalenergy,           /* */
      colibre_cooling_N_metallicity,              /* */
      colibre_cooling_N_density,                  /* */
      colibre_cooling_N_cooltypes);               /* */

  /* Gamma / n_H**2 */
  const double heating_rate = interpolation4d_plus_summation(
      cooling->table.Uheating, weights_heating,   /* */
      element_H, colibre_cooling_N_heattypes - 3, /* */
      red_index, U_index, met_index, n_H_index,   /* */
      d_red, d_U, d_met, d_n_H,                   /* */
      colibre_cooling_N_redshifts,                /* */
      colibre_cooling_N_internalenergy,           /* */
      colibre_cooling_N_metallicity,              /* */
      colibre_cooling_N_density,                  /* */
      colibre_cooling_N_heattypes);               /* */

  /* Temperature from internal energy */
  const double logtemp =
      interpolation_4d(cooling->table.T_from_U,                  /* */
                       red_index, U_index, met_index, n_H_index, /* */
                       d_red, d_U, d_met, d_n_H,                 /* */
                       colibre_cooling_N_redshifts,              /* */
                       colibre_cooling_N_internalenergy,         /* */
                       colibre_cooling_N_metallicity,            /* */
                       colibre_cooling_N_density);               /* */
                                                                 /* */
  const double temp = exp10(logtemp);

  /* Compton cooling/heating */
  double Compton_cooling_rate = 0.;
  if (onlyicool == 0 || (onlyicool == 1 && icool == cooltype_Compton)) {

    const double zp1 = 1. + redshift;
    const double zp1p2 = zp1 * zp1;
    const double zp1p4 = zp1p2 * zp1p2;

    /* CMB temperature at this redshift */
    const double T_CMB = cooling->T_CMB_0 * zp1;

    /* Analytic Compton cooling rate: Lambda_Compton / n_H**2 */
    Compton_cooling_rate = cooling->compton_rate_cgs * (temp - T_CMB) * zp1p4 *
                           electron_fraction / n_H_cgs;
  }

  /* Return the net heating rate (Lambda_heat - Lambda_cool) */
  return heating_rate - cooling_rate - Compton_cooling_rate;
}

/**
 * @brief Computes the net cooling rate (cooling - heating) for a given element
 * abundance ratio, temperature, redshift, and density. The unit of the net
 * cooling rate is Lambda / nH**2 [erg cm^3 s-1] and all input values are in
 * cgs. The Compton cooling is not taken from the tables but calculated
 * analytically and added separately
 *
 * @param log_T_cgs Log base 10 of temperature in K
 * @param redshift Current redshift
 * @param n_H_cgs Hydrogen number density in cgs
 * @param abundance_ratio Abundance ratio for each element x relative to solar
 * @param n_H_index Index along the Hydrogen number density dimension
 * @param d_n_H Offset between Hydrogen density and table[n_H_index]
 * @param met_index Index along the metallicity dimension
 * @param d_met Offset between metallicity and table[met_index]
 * @param red_index Index along redshift dimension
 * @param d_red Offset between redshift and table[red_index]
 * @param cooling #cooling_function_data structure
 *
 * @param onlyicool if true / 1 only plot cooling channel icool
 * @param onlyiheat if true / 1 only plot cooling channel iheat
 * @param icool cooling channel to be used
 * @param iheat heating channel to be used
 *
 * Throughout the code: onlyicool = onlyiheat = icool = iheat = 0
 * These are only used for testing: examples/CoolingRates/CoolingRatesPS2020
 */
INLINE static double colibre_cooling_rate_temperature(
    const double log_T_cgs, const double redshift, const double n_H_cgs,
    const float abundance_ratio[colibre_cooling_N_elementtypes],
    const int n_H_index, const float d_n_H, const int met_index,
    const float d_met, const int red_index, const float d_red,
    const struct cooling_function_data *cooling, const int onlyicool,
    const int onlyiheat, const int icool, const int iheat) {

  /* Set weights for cooling rates */
  float weights_cooling[colibre_cooling_N_cooltypes - 2];
  for (int i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {

    if (i < colibre_cooling_N_elementtypes) {
      /* Use abundance ratios */
      weights_cooling[i] = abundance_ratio[i];
    } else if (i == cooltype_Compton) {
      /* added analytically later, do not use value from table*/
      weights_cooling[i] = 0.f;
    } else {
      /* use same abundances as in the tables */
      weights_cooling[i] = 1.f;
    }
  }

  /* If we care only about one channel, cancel all the other ones */
  if (onlyicool != 0) {
    for (int i = 0; i < colibre_cooling_N_cooltypes - 2; i++) {
      if (i != icool) weights_cooling[i] = 0.f;
    }
  }

  /* Set weights for heating rates */
  float weights_heating[colibre_cooling_N_heattypes - 2];
  for (int i = 0; i < colibre_cooling_N_heattypes - 2; i++) {
    if (i < colibre_cooling_N_elementtypes) {
      weights_heating[i] = abundance_ratio[i];
    } else {
      weights_heating[i] = 1.f; /* use same abundances as in the tables */
    }
  }

  /* If we care only about one channel, cancel all the other ones */
  if (onlyiheat != 0) {
    for (int i = 0; i < colibre_cooling_N_heattypes - 2; i++) {
      if (i != iheat) weights_heating[i] = 0.f;
    }
  }

  /* Get index of T along the internal energy axis */
  int T_index;
  float d_T;
  get_index_1d(cooling->Temp, colibre_cooling_N_temperature, log_T_cgs,
               &T_index, &d_T);

  /* n_e / n_H */
  const double electron_fraction = interpolation4d_plus_summation(
      cooling->table.Telectron_fraction, abundance_ratio, /* */
      element_H, colibre_cooling_N_electrontypes - 4,     /* */
      red_index, T_index, met_index, n_H_index,           /* */
      d_red, d_T, d_met, d_n_H,                           /* */
      colibre_cooling_N_redshifts,                        /* */
      colibre_cooling_N_temperature,                      /* */
      colibre_cooling_N_metallicity,                      /* */
      colibre_cooling_N_density,                          /* */
      colibre_cooling_N_electrontypes);                   /* */

  /* Lambda / n_H**2 */
  const double cooling_rate = interpolation4d_plus_summation(
      cooling->table.Tcooling, weights_cooling,   /* */
      element_H, colibre_cooling_N_cooltypes - 3, /* */
      red_index, T_index, met_index, n_H_index,   /* */
      d_red, d_T, d_met, d_n_H,                   /* */
      colibre_cooling_N_redshifts,                /* */
      colibre_cooling_N_temperature,              /* */
      colibre_cooling_N_metallicity,              /* */
      colibre_cooling_N_density,                  /* */
      colibre_cooling_N_cooltypes);               /* */

  /* Gamma / n_H**2 */
  const double heating_rate = interpolation4d_plus_summation(
      cooling->table.Theating, weights_heating,   /* */
      element_H, colibre_cooling_N_heattypes - 3, /* */
      red_index, T_index, met_index, n_H_index,   /* */
      d_red, d_T, d_met, d_n_H,                   /* */
      colibre_cooling_N_redshifts,                /* */
      colibre_cooling_N_temperature,              /* */
      colibre_cooling_N_metallicity,              /* */
      colibre_cooling_N_density,                  /* */
      colibre_cooling_N_heattypes);               /* */

  const double temp = exp10(log_T_cgs);

  double Compton_cooling_rate = 0.;
  if (onlyicool == 0 || (onlyicool == 1 && icool == cooltype_Compton)) {

    const double zp1 = 1. + redshift;
    const double zp1p2 = zp1 * zp1;
    const double zp1p4 = zp1p2 * zp1p2;

    /* CMB temperature at this redshift */
    const double T_CMB = cooling->T_CMB_0 * zp1;

    /* Analytic Compton cooling rate: Lambda_Compton / n_H**2 */
    Compton_cooling_rate = cooling->compton_rate_cgs * (temp - T_CMB) * zp1p4 *
                           electron_fraction / n_H_cgs;
  }

  /* Return the net heating rate (Lambda_heat - Lambda_cool) */
  return heating_rate - cooling_rate - Compton_cooling_rate;
}

#endif /* SWIFT_PS2020_COOLING_RATES_H */
