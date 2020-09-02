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
#ifndef SWIFT_GADGET2_LOGGER_HYDRO_H
#define SWIFT_GADGET2_LOGGER_HYDRO_H

#include "../config.h"

/* local includes */
#include "hydro_logger.h"
#include "logger_loader_io.h"
#include "logger_python_tools.h"

/* Index of the mask in the header mask array */
extern int hydro_logger_local_to_global[hydro_logger_field_count];

/* Size for each mask */
extern const int hydro_logger_field_size[hydro_logger_field_count];

/**
 * @brief Create the link between the fields and their derivatives.
 *
 * @param head The #header.
 */
__attribute__((always_inline)) INLINE static void
hydro_logger_reader_link_derivatives(struct header *head) {

  /* Set the first and second derivatives */
  const int pos_id =
      hydro_logger_local_to_global[hydro_logger_field_coordinates];
  const int vel_id =
      hydro_logger_local_to_global[hydro_logger_field_velocities];
  const int acc_id =
      hydro_logger_local_to_global[hydro_logger_field_accelerations];

  /* Coordinates */
  header_set_first_derivative(head, pos_id, vel_id);
  header_set_second_derivative(head, pos_id, acc_id);

  /* Velocities */
  header_set_first_derivative(head, vel_id, acc_id);
}

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
 * #hydro_logger_fields).
 */
__attribute__((always_inline)) INLINE static void
hydro_logger_interpolate_field(const double t_before,
                               const struct logger_field *restrict before,
                               const double t_after,
                               const struct logger_field *restrict after,
                               void *restrict output, const double t,
                               const int field) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Check the times */
  if (t_before > t || t_after < t) {
    error(
        "The times for the interpolation are not correct"
        " %g < %g < %g.",
        t_before, t, t_after);
  }
#endif

  /* Compute the interpolation scaling. */
  const double wa = (t - t_before) / (t_after - t_before);
  const double wb = 1. - wa;

  switch (field) {
    /* Do the position */
    case hydro_logger_field_coordinates:
      for (int i = 0; i < 3; i++) {
        double *x = (double *)output;
        const double *x_bef = (double *)before->field;
        const float *v_bef = (float *)before->first_deriv;
        const float *a_bef = (float *)before->second_deriv;

        const double *x_aft = (double *)after->field;
        const float *v_aft = (float *)after->first_deriv;
        const float *a_aft = (float *)after->second_deriv;

        /* Use quintic hermite spline. */
        if (v_bef && v_aft && a_bef && a_aft) {
          x[i] = logger_tools_quintic_hermite_spline(
              t_before, x_bef[i], v_bef[i], a_bef[i], t_after, x_aft[i],
              v_aft[i], a_aft[i], t);
        }
        /* Use cubic hermite spline. */
        else if (v_bef && v_aft) {
          x[i] = logger_tools_cubic_hermite_spline(
              t_before, x_bef[i], v_bef[i], t_after, x_aft[i], v_aft[i], t);
        }
        /* Use linear interpolation. */
        else {
          x[i] = wa * x_aft[i] + wb * x_bef[i];
        }
      }
      break;
      /* Do the velocity */
    case hydro_logger_field_velocities:
      for (int i = 0; i < 3; i++) {
        float *v = (float *)output;
        const float *v_bef = (float *)before->field;
        const float *a_bef = (float *)before->first_deriv;

        const float *v_aft = (float *)after->field;
        const float *a_aft = (float *)after->first_deriv;

        /* Use a cubic hermite spline. */
        if (a_bef && a_aft) {
          v[i] = logger_tools_cubic_hermite_spline(
              t_before, v_bef[i], a_bef[i], t_after, v_aft[i], a_aft[i], t);
        }
        /* Use linear interpolation. */
        else {
          v[i] = wa * v_aft[i] + wb * v_bef[i];
        }
      }
      break;
    case hydro_logger_field_accelerations:
      /* interpolate vectors. */
      for (int i = 0; i < 3; i++) {
        float *a = (float *)output;
        const float *a_bef = (float *)before->field;
        const float *a_aft = (float *)after->field;
        a[i] = wa * a_aft[i] + wb * a_bef[i];
      }
      break;
      /* Do the linear interpolation of float. */
    case hydro_logger_field_smoothing_lengths:
    case hydro_logger_field_entropies:
    case hydro_logger_field_densities:
    case hydro_logger_field_masses:
      ((float *)output)[0] =
          wa * ((float *)after->field)[0] + wb * ((float *)before->field)[0];
      break;
      /* Check the ids */
    case hydro_logger_field_particle_ids:
      if (*(long long *)after->field != *(long long *)before->field) {
        error("Interpolating different particles");
      }
      *(long long *)output = *(long long *)after->field;
      break;
    default:
      error("Not implemented");
  }
}

#ifdef HAVE_PYTHON
__attribute__((always_inline)) INLINE static void hydro_logger_generate_python(
    struct logger_python_field *fields) {

  fields[hydro_logger_field_coordinates] =
      logger_loader_python_field(/* Dimension */ 3, NPY_DOUBLE);
  fields[hydro_logger_field_velocities] =
      logger_loader_python_field(/* Dimension */ 3, NPY_FLOAT32);
  fields[hydro_logger_field_accelerations] =
      logger_loader_python_field(/* Dimension */ 3, NPY_FLOAT32);
  fields[hydro_logger_field_masses] =
      logger_loader_python_field(/* Dimension */ 1, NPY_FLOAT32);
  fields[hydro_logger_field_smoothing_lengths] =
      logger_loader_python_field(/* Dimension */ 1, NPY_FLOAT32);
  fields[hydro_logger_field_entropies] =
      logger_loader_python_field(/* Dimension */ 1, NPY_FLOAT32);
  fields[hydro_logger_field_particle_ids] =
      logger_loader_python_field(/* Dimension */ 1, NPY_LONGLONG);
  fields[hydro_logger_field_densities] =
      logger_loader_python_field(/* Dimension */ 1, NPY_FLOAT32);
}

#endif  // HAVE_PYTHON
#endif  /* SWIFT_GADGET2_LOGGER_HYDRO_H */
