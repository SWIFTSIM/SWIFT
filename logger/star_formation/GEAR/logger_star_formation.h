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
#ifndef SWIFT_GEAR_LOGGER_STAR_FORMATION_H
#define SWIFT_GEAR_LOGGER_STAR_FORMATION_H

#include "../config.h"

/* local includes */
#include "logger_interpolation.h"
#include "logger_loader_io.h"
#include "logger_parameters.h"
#include "logger_python_tools.h"
#include "star_formation_particle_logger.h"

/* Index of the mask in the header mask array */
extern int
    star_formation_logger_local_to_global[star_formation_logger_field_count];

/* Size for each mask */
extern const int
    star_formation_logger_field_size[star_formation_logger_field_count];

/**
 * @brief Create the link between the fields and their derivatives.
 *
 * @param head The #header.
 */
__attribute__((always_inline)) INLINE static void
star_formation_logger_reader_link_derivatives(struct header *head) {}

/**
 * @brief Interpolate a field of the particle at the given time.
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
 * #star_formation_logger_fields).
 * @param params The simulation's #logger_parameters.
 */
__attribute__((always_inline)) INLINE static void
star_formation_logger_interpolate_field(
    const double t_before, const struct logger_field *restrict before,
    const double t_after, const struct logger_field *restrict after,
    void *restrict output, const double t, const int field,
    const struct logger_parameters *params) {

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
    case star_formation_logger_field_all:
      interpolate_linear_float_ND(t_before, before, t_after, after, output, t,
                                  /* Dimension */ 2, /* periodic= */ 0, params);

      /* Deal with the progenitor id */
      long long *proj_id =
          (long long *)((void *)before->field + 2 * sizeof(float));
      output += 2 * sizeof(float);
      *(long long *)output = *proj_id;
      break;

    default:
      error("Not implemented");
  }
}

#ifdef HAVE_PYTHON
__attribute__((always_inline)) INLINE static void
star_formation_logger_generate_python(struct logger_python_field *fields) {
  fields[star_formation_logger_field_all] =
      logger_loader_python_field(/* Dimension */ 2, NPY_FLOAT);

  logger_loader_python_field_add_subfield(
      &fields[star_formation_logger_field_all], "BirthDensities",
      /* offset */ 0, /* Dimension */ 1, NPY_FLOAT32);

  logger_loader_python_field_add_subfield(
      &fields[star_formation_logger_field_all], "BirthMasses",
      /* offset */ sizeof(float), /* Dimension */ 1, NPY_FLOAT32);

  logger_loader_python_field_add_subfield(
      &fields[star_formation_logger_field_all], "ProgenitorIDs",
      /* offset */ 2 * sizeof(float), /* Dimension */ 1, NPY_LONGLONG);
}

#endif  // HAVE_PYTHON
#endif  // SWIFT_GEAR_LOGGER_STAR_FORMATION_H
