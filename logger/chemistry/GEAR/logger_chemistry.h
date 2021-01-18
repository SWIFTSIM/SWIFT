/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_LOGGER_CHEMISTRY_H
#define SWIFT_GEAR_LOGGER_CHEMISTRY_H

#include "../config.h"

/* local includes */
#include "chemistry_logger.h"
#include "logger_interpolation.h"
#include "logger_loader_io.h"
#include "logger_python_tools.h"

/* Index of the mask in the header mask array */
extern int
    chemistry_logger_local_to_global_part[chemistry_logger_field_part_count];
extern int
    chemistry_logger_local_to_global_spart[chemistry_logger_field_spart_count];

/* Size for each mask */
extern const int
    chemistry_logger_field_size_part[chemistry_logger_field_part_count];
extern const int
    chemistry_logger_field_size_spart[chemistry_logger_field_spart_count];

/**
 * @brief Create the link between the fields and their derivatives for #part.
 *
 * @param head The #header.
 */
__attribute__((always_inline)) INLINE static void
chemistry_logger_reader_link_derivatives_part(struct header *head) {}

/**
 * @brief Create the link between the fields and their derivatives for #spart.
 *
 * @param head The #header.
 */
__attribute__((always_inline)) INLINE static void
chemistry_logger_reader_link_derivatives_spart(struct header *head) {}

/**
 * @brief Interpolate a field of the #part at the given time.
 * Here we use a linear interpolation for most of the fields.
 * For the position (velocity), we use a quintic (cubic) hermite interpolation
 * based on the positions, velocities and accelerations at the time of the two
 * particles.
 *
 * @param before Pointer to the #logger_field at a time < t.
 * @param after Pointer to the #logger_field at a time > t.
 * @param otuput Pointer to the output value.
 * @param t_before Time of field_before (< t).
 * @param t_after Time of field_after (> t).
 * @param t Requested time.
 * @param field The field to reconstruct (follows the order of
 * #chemistry_logger_fields).
 */
__attribute__((always_inline)) INLINE static void
chemistry_logger_interpolate_field_part(
    const double t_before, const struct logger_field *restrict before,
    const double t_after, const struct logger_field *restrict after,
    void *restrict output, const double t, const int field) {
#ifdef SWIFT_DEBUG_CHECKS
  /* Check the times */
  if (t_before > t || t_after < t) {
    error(
        "The times for the interpolation are not correct"
        " %g < %g < %g.",
        t_before, t, t_after);
  }
#endif

  switch (field) {
    case chemistry_logger_field_part_all:
      interpolate_linear_double_ND(t_before, before, t_after, after, output, t,
                                   2 * GEAR_CHEMISTRY_ELEMENT_COUNT);
      break;

    default:
      error("Not implemented");
  }
}

/**
 * @brief Interpolate a field of the #spart at the given time.
 * Here we use a linear interpolation for most of the fields.
 * For the position (velocity), we use a quintic (cubic) hermite interpolation
 * based on the positions, velocities and accelerations at the time of the two
 * particles.
 *
 * @param before Pointer to the #logger_field at a time < t.
 * @param after Pointer to the #logger_field at a time > t.
 * @param otuput Pointer to the output value.
 * @param t_before Time of field_before (< t).
 * @param t_after Time of field_after (> t).
 * @param t Requested time.
 * @param field The field to reconstruct (follows the order of
 * #chemistry_logger_fields).
 */
__attribute__((always_inline)) INLINE static void
chemistry_logger_interpolate_field_spart(
    const double t_before, const struct logger_field *restrict before,
    const double t_after, const struct logger_field *restrict after,
    void *restrict output, const double t, const int field) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Check the times */
  if (t_before > t || t_after < t) {
    error(
        "The times for the interpolation are not correct"
        " %g < %g < %g.",
        t_before, t, t_after);
  }
#endif

  switch (field) {
    case chemistry_logger_field_spart_metal_mass_fractions:
      interpolate_linear_double_ND(t_before, before, t_after, after, output, t,
                                   GEAR_CHEMISTRY_ELEMENT_COUNT);
      break;

    default:
      error("Not implemented");
  }
}

#ifdef HAVE_PYTHON
/**
 * @brief Defines the different arrays for python.
 */
__attribute__((always_inline)) INLINE static void
chemistry_logger_generate_python_part(struct logger_python_field *fields) {
  fields[chemistry_logger_field_part_all] =
      logger_loader_python_field(2 * GEAR_CHEMISTRY_ELEMENT_COUNT, NPY_DOUBLE);
}

/**
 * @brief Defines the different arrays for python.
 */
__attribute__((always_inline)) INLINE static void
chemistry_logger_generate_python_spart(struct logger_python_field *fields) {
  fields[chemistry_logger_field_spart_metal_mass_fractions] =
      logger_loader_python_field(GEAR_CHEMISTRY_ELEMENT_COUNT, NPY_DOUBLE);
}

#endif  // HAVE_PYTHON
#endif  // SWIFT_GEAR_LOGGER_CHEMISTRY_H
