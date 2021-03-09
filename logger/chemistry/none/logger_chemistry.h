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
#ifndef SWIFT_NONE_LOGGER_CHEMISTRY_H
#define SWIFT_NONE_LOGGER_CHEMISTRY_H

#include "../config.h"

/* local includes */
#include "chemistry_logger.h"
#include "logger_interpolation.h"
#include "logger_loader_io.h"
#include "logger_parameters.h"
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
 * @param params The simulation's #logger_parameters.
 */
__attribute__((always_inline)) INLINE static void
chemistry_logger_interpolate_field_part(
    const double t_before, const struct logger_field *restrict before,
    const double t_after, const struct logger_field *restrict after,
    void *restrict output, const double t, const int field,
    const struct logger_parameters *params) {}

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
 * @param params The simulation's #logger_parameters.
 */
__attribute__((always_inline)) INLINE static void
chemistry_logger_interpolate_field_spart(
    const double t_before, const struct logger_field *restrict before,
    const double t_after, const struct logger_field *restrict after,
    void *restrict output, const double t, const int field,
    const struct logger_parameters *params) {}

#ifdef HAVE_PYTHON
/**
 * @brief Defines the different arrays for python.
 */
__attribute__((always_inline)) INLINE static void
chemistry_logger_generate_python_part(struct logger_python_field *fields) {}

/**
 * @brief Defines the different arrays for python.
 */
__attribute__((always_inline)) INLINE static void
chemistry_logger_generate_python_spart(struct logger_python_field *fields) {}

#endif  // HAVE_PYTHON
#endif  // SWIFT_NONE_LOGGER_CHEMISTRY_H
