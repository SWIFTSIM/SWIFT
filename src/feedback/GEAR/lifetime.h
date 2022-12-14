/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_LIFETIME_GEAR_H
#define SWIFT_LIFETIME_GEAR_H

#include "hdf5_functions.h"
#include "stellar_evolution_struct.h"

/**
 * @brief Print the lifetime model.
 *
 * @param lf The #lifetime.
 */
__attribute__((always_inline)) INLINE static void lifetime_print(
    const struct lifetime* lf) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  message("Quadratic terms: %.2g %.2g %.2g", lf->quadratic[0], lf->quadratic[1],
          lf->quadratic[2]);

  message("Linear terms: %.2g %.2g %.2g", lf->linear[0], lf->linear[1],
          lf->linear[2]);

  message("Constant terms: %.2g %.2g %.2g", lf->constant[0], lf->constant[1],
          lf->constant[2]);
}

/**
 * @brief Compute the lifetime of a star.
 *
 * @param life The #lifetime model.
 * @param log_mass The star's mass (in log10(solMass)).
 * @param metallicity The star's metallicity.
 *
 * @return The star's lifetime (in log10(Myr)).
 */
__attribute__((always_inline)) INLINE static float
lifetime_get_log_lifetime_from_mass(const struct lifetime* life, float log_mass,
                                    float metallicity) {

  /* Compute quadratic term */
  const float quadratic =
      (life->quadratic[0] * metallicity + life->quadratic[1]) * metallicity +
      life->quadratic[2];
  /* Compute linear term */
  const float linear =
      (life->linear[0] * metallicity + life->linear[1]) * metallicity +
      life->linear[2];
  /* Compute constant term */
  const float constant =
      (life->constant[0] * metallicity + life->constant[1]) * metallicity +
      life->constant[2];

  /* Compute lifetime */
  return (quadratic * log_mass + linear) * log_mass + constant;
}

/**
 * @brief Compute the mass of a star with a given lifetime
 *
 * @param life The #lifetime model.
 * @param log_time The star's lifetime (in log10(Myr)).
 * @param metallicity The star's metallicity.
 *
 * @return The star's mass (in log10(solMass))
 */
__attribute__((always_inline)) INLINE static float
lifetime_get_log_mass_from_lifetime(const struct lifetime* life, float log_time,
                                    float metallicity) {

  /* Compute quadratic term */
  const float quadratic =
      (life->quadratic[0] * metallicity + life->quadratic[1]) * metallicity +
      life->quadratic[2];
  /* Compute linear term */
  const float linear =
      (life->linear[0] * metallicity + life->linear[1]) * metallicity +
      life->linear[2];
  /* Compute constant term */
  const float constant =
      (life->constant[0] * metallicity + life->constant[1]) * metallicity +
      life->constant[2];

  /* Compute the "c" with the time */
  const float c_t = constant - log_time;

  /* Use the quadratic formula to find the mass */
  if (quadratic != 0) {
    const float delta = linear * linear - 4 * quadratic * c_t;

    /* Avoid complex number should not happen in real simulation */
    if (delta < 0) {
      return -linear / (2. * quadratic);
    } else {
      return (-linear - sqrt(delta)) / (2. * quadratic);
    }
  } else {
    return -c_t / linear;
  }
}

/**
 * @brief Read lifetime parameters from tables.
 *
 * @param lt The #lifetime.
 * @param params The #swift_params.
 * @param filename The filename of the chemistry table.
 */
__attribute__((always_inline)) INLINE static void lifetime_read_from_tables(
    struct lifetime* lt, struct swift_params* params, const char* filename) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(filename, "Data/LifeTimes", &file_id, &group_id);

  /* Allocate the temporary array */
  float* tmp;
  if ((tmp = (float*)malloc(sizeof(float) * 9)) == NULL)
    error("Failed to allocate the temporary array.");

  /* Read the coefficients */
  io_read_array_dataset(group_id, "coeff_z", FLOAT, tmp, 9);

  /* Copy the coefficents */
  const int dim = 3;
  for (int i = 0; i < dim; i++) {
    lt->quadratic[i] = tmp[i];
    lt->linear[i] = tmp[i + dim];
    lt->constant[i] = tmp[i + 2 * dim];
  }

  /* Cleanup everything */
  free(tmp);
  h5_close_group(file_id, group_id);
}

/**
 * @brief Inititialize the Lifetime.
 *
 * @param lt The #lifetime.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The #swift_params.
 * @param filename The filename of the chemistry table.
 */
__attribute__((always_inline)) INLINE static void lifetime_init(
    struct lifetime* lt, const struct phys_const* phys_const,
    const struct unit_system* us, struct swift_params* params,
    const char* filename) {

  /* Read params from yields table */
  lifetime_read_from_tables(lt, params, filename);

  /* Change units from yr into Myr */
  const int dim = 3;
  lt->constant[dim - 1] -= 6;
}

/**
 * @brief Write a lifetime struct to the given FILE as a stream of bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param lt the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
__attribute__((always_inline)) INLINE static void lifetime_dump(
    const struct lifetime* lt, FILE* stream, const struct stellar_model* sm) {

  /* Nothing to do here */
}

/**
 * @brief Restore a lifetime struct from the given FILE as a stream of
 * bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param lt the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
__attribute__((always_inline)) INLINE static void lifetime_restore(
    struct lifetime* lt, FILE* stream, const struct stellar_model* sm) {

  /* Nothing to do here */
}

/**
 * @brief Clean the allocated memory.
 *
 * @param lifetime the #lifetime.
 */
__attribute__((always_inline)) INLINE static void lifetime_clean(
    struct lifetime* lifetime) {}

#endif  // SWIFT_LIFETIME_GEAR_H
