/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (mschaller@lorentz.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_TABLES_TREECOOL_H
#define SWIFT_COOLING_TABLES_TREECOOL_H

/**
 * @file src/cooling/TREECOOL/cooling_tables.h
 * @brief Reading and interpolation of the TREECOOL table describing the
 * evolution of the UV background with redshift.
 *
 * A TREECOOL file is a plain text table with one row per redshift and seven
 * columns:
 *
 *   log10(1 + z)
 *   Gamma_HI   Gamma_HeI   Gamma_HeII      [s^-1]
 *   eps_HI     eps_HeI     eps_HeII        [erg * s^-1]
 *
 * where the Gamma are the photo-ionization rates and the eps the photo-heating
 * rates per ion of the corresponding species. Such tables are distributed with
 * most of the widely used UV background models; an example file,
 * TREECOOL_UV_background.txt, is provided in
 * examples/Cooling/ConstantCosmoTempEvolution/.
 */

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <math.h>
#include <stdio.h>

/* Local includes. */
#include "cooling_properties.h"
#include "error.h"
#include "exp10.h"
#include "inline.h"

/**
 * @brief Read the TREECOOL table from the file given in the parameter file.
 *
 * Many TREECOOL files are padded with rows of zeros beyond the redshift at
 * which the UV background model stops. We stop reading at the first such row
 * and treat everything above that redshift as having no UV background.
 *
 * @param cooling The #cooling_function_data to fill in.
 */
INLINE static void treecool_read_table(struct cooling_function_data *cooling) {

  FILE *file = fopen(cooling->TREECOOL_file, "r");
  if (file == NULL)
    error("Cannot open the TREECOOL file '%s'", cooling->TREECOOL_file);

  int count = 0;
  double log10_1_plus_z, gamma_H0, gamma_He0, gamma_Hep;
  double epsilon_H0, epsilon_He0, epsilon_Hep;

  while (fscanf(file, "%lg %lg %lg %lg %lg %lg %lg", &log10_1_plus_z, &gamma_H0,
                &gamma_He0, &gamma_Hep, &epsilon_H0, &epsilon_He0,
                &epsilon_Hep) == 7) {

    /* Ignore the padding at the end of the table */
    if (gamma_H0 <= 0.) break;

    if (count >= treecool_cooling_max_N_redshifts)
      error(
          "The TREECOOL file '%s' contains more than %d entries. Increase "
          "treecool_cooling_max_N_redshifts.",
          cooling->TREECOOL_file, treecool_cooling_max_N_redshifts);

    if (count > 0 &&
        log10_1_plus_z <= cooling->TREECOOL_log10_1_plus_z[count - 1])
      error(
          "The redshifts in the TREECOOL file '%s' are not in increasing "
          "order (entry %d)",
          cooling->TREECOOL_file, count);

    cooling->TREECOOL_log10_1_plus_z[count] = log10_1_plus_z;
    cooling->TREECOOL_gamma_H0_cgs[count] = gamma_H0;
    cooling->TREECOOL_gamma_He0_cgs[count] = gamma_He0;
    cooling->TREECOOL_gamma_Hep_cgs[count] = gamma_Hep;
    cooling->TREECOOL_epsilon_H0_cgs[count] = epsilon_H0;
    cooling->TREECOOL_epsilon_He0_cgs[count] = epsilon_He0;
    cooling->TREECOOL_epsilon_Hep_cgs[count] = epsilon_Hep;

    ++count;
  }

  fclose(file);

  if (count < 2)
    error("The TREECOOL file '%s' contains fewer than 2 usable entries",
          cooling->TREECOOL_file);

  cooling->N_redshifts = count;
}

/**
 * @brief Interpolate one column of the TREECOOL table.
 *
 * The interpolation is linear in log10(1 + z) and in the log of the rate
 * itself, as the rates vary by many orders of magnitude across the table.
 *
 * @param table The column to interpolate.
 * @param index The index of the entry just below the requested redshift.
 * @param f_lo The weight of the entry just below the requested redshift.
 * @param f_hi The weight of the entry just above the requested redshift.
 */
__attribute__((always_inline)) INLINE static double treecool_interpolate_table(
    const double *table, const int index, const double f_lo,
    const double f_hi) {

  /* A rate that vanishes at either end cannot be interpolated in log-space.
   * Fall back to a linear interpolation in that case. */
  if (table[index] <= 0. || table[index + 1] <= 0.)
    return f_lo * table[index] + f_hi * table[index + 1];

  return exp10(f_lo * log10(table[index]) + f_hi * log10(table[index + 1]));
}

/**
 * @brief Set the photo-ionization and photo-heating rates at a given redshift.
 *
 * The rates are switched off entirely above the highest redshift covered by
 * the table and above the (optional) redshift at which the user asked for the
 * UV background to be switched on.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param redshift The current redshift.
 */
INLINE static void treecool_set_UV_background(
    struct cooling_function_data *cooling, const double redshift) {

  /* Start from no UV background at all */
  cooling->UV_background_on = 0;
  cooling->gamma_H0_cgs = 0.;
  cooling->gamma_He0_cgs = 0.;
  cooling->gamma_Hep_cgs = 0.;
  cooling->epsilon_H0_cgs = 0.;
  cooling->epsilon_He0_cgs = 0.;
  cooling->epsilon_Hep_cgs = 0.;

  /* Have the sources switched on yet? */
  if (redshift > cooling->UV_background_start_redshift) return;

  const double log10_1_plus_z = log10(1. + redshift);

  /* Above the table there is no UV background */
  if (log10_1_plus_z >
      cooling->TREECOOL_log10_1_plus_z[cooling->N_redshifts - 1])
    return;

  /* Find the entry just below the current redshift. Below the first entry of
   * the table we simply use its first interval. */
  int index = 0;
  for (int i = 0; i < cooling->N_redshifts - 1; ++i) {
    if (cooling->TREECOOL_log10_1_plus_z[i] < log10_1_plus_z)
      index = i;
    else
      break;
  }

  const double d_lo = log10_1_plus_z - cooling->TREECOOL_log10_1_plus_z[index];
  const double d_hi =
      cooling->TREECOOL_log10_1_plus_z[index + 1] - log10_1_plus_z;
  const double inv_d = 1. / (d_lo + d_hi);

  /* Note the swap: being close to the lower entry means a small d_lo, hence a
   * large weight d_hi for that entry. */
  const double f_lo = d_hi * inv_d;
  const double f_hi = d_lo * inv_d;

  cooling->gamma_H0_cgs = treecool_interpolate_table(
      cooling->TREECOOL_gamma_H0_cgs, index, f_lo, f_hi);
  cooling->gamma_He0_cgs = treecool_interpolate_table(
      cooling->TREECOOL_gamma_He0_cgs, index, f_lo, f_hi);
  cooling->gamma_Hep_cgs = treecool_interpolate_table(
      cooling->TREECOOL_gamma_Hep_cgs, index, f_lo, f_hi);
  cooling->epsilon_H0_cgs = treecool_interpolate_table(
      cooling->TREECOOL_epsilon_H0_cgs, index, f_lo, f_hi);
  cooling->epsilon_He0_cgs = treecool_interpolate_table(
      cooling->TREECOOL_epsilon_He0_cgs, index, f_lo, f_hi);
  cooling->epsilon_Hep_cgs = treecool_interpolate_table(
      cooling->TREECOOL_epsilon_Hep_cgs, index, f_lo, f_hi);

  cooling->UV_background_on = 1;
}

#endif /* SWIFT_COOLING_TABLES_TREECOOL_H */
