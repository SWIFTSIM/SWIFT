/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GADGET2_HYDRO_LOGGER_H
#define SWIFT_GADGET2_HYDRO_LOGGER_H

/* Include configuration */
#include "config.h"

#ifdef WITH_LOGGER

/* Include the particles */
#include "align.h"
#include "hydro_part.h"

/* Include the hydro functions */
#include "hydro.h"

/* Include the logger functions */
#include "logger_io.h"

/*
 * List of all possible mask.
 * Outside the module, only hydro_logger_field_count is used.
 */
enum hydro_logger_fields {
  hydro_logger_field_coordinates = 0,
  hydro_logger_field_velocities,
  hydro_logger_field_accelerations,
  hydro_logger_field_masses,
  hydro_logger_field_smoothing_lengths,
  hydro_logger_field_entropies,
  hydro_logger_field_particle_ids,
  hydro_logger_field_densities,
  hydro_logger_field_count,
};

/* Name of each possible mask. */
extern const char *hydro_logger_field_names[hydro_logger_field_count];

/**
 * @brief Initialize the logger.
 *
 * WARNING: The order should be the same in all the functions and
 * #hydro_logger_fields!
 *
 * @param mask_data Data for each type of mask.
 *
 * @return Number of masks used.
 */
INLINE static int hydro_logger_writer_populate_mask_data(
    struct mask_data *mask_data) {
  mask_data[hydro_logger_field_coordinates] = logger_create_mask_entry(
      hydro_logger_field_names[hydro_logger_field_coordinates],
      3 * sizeof(double));

  mask_data[hydro_logger_field_velocities] = logger_create_mask_entry(
      hydro_logger_field_names[hydro_logger_field_velocities],
      3 * sizeof(float));

  mask_data[hydro_logger_field_accelerations] = logger_create_mask_entry(
      hydro_logger_field_names[hydro_logger_field_accelerations],
      3 * sizeof(float));

  mask_data[hydro_logger_field_masses] = logger_create_mask_entry(
      hydro_logger_field_names[hydro_logger_field_masses], sizeof(float));

  mask_data[hydro_logger_field_smoothing_lengths] = logger_create_mask_entry(
      hydro_logger_field_names[hydro_logger_field_smoothing_lengths],
      sizeof(float));

  mask_data[hydro_logger_field_entropies] = logger_create_mask_entry(
      hydro_logger_field_names[hydro_logger_field_entropies], sizeof(float));

  mask_data[hydro_logger_field_particle_ids] = logger_create_mask_entry(
      hydro_logger_field_names[hydro_logger_field_particle_ids],
      sizeof(long long));

  mask_data[hydro_logger_field_densities] = logger_create_mask_entry(
      hydro_logger_field_names[hydro_logger_field_densities], sizeof(float));

  return hydro_logger_field_count;
}

/**
 * @brief Generates the mask and compute the size of the record.
 *
 * WARNING: The order should be the same in all the functions and
 * #hydro_logger_fields!
 *
 * @param masks The list of masks (same order than in #hydro_logger_init).
 * @param part The #part that will be written.
 * @param xpart The #xpart that will be written.
 * @param write_all Are we forcing to write all the fields?
 *
 * @param buffer_size (out) The requested size for the buffer.
 * @param mask (out) The mask that will be written.
 */
INLINE static void hydro_logger_compute_size_and_mask(
    const struct mask_data *masks, const struct part *part,
    const struct xpart *xpart, const int write_all, size_t *buffer_size,
    unsigned int *mask) {

  /* Here you can decide your own writing logic */

  /* Add the coordinates. */
  *mask |= logger_add_field_to_mask(masks[hydro_logger_field_coordinates],
                                    buffer_size);

  /* Add the velocities. */
  *mask |= logger_add_field_to_mask(masks[hydro_logger_field_velocities],
                                    buffer_size);

  /* Add the accelerations. */
  *mask |= logger_add_field_to_mask(masks[hydro_logger_field_accelerations],
                                    buffer_size);

  /* Add the masses. */
  *mask |=
      logger_add_field_to_mask(masks[hydro_logger_field_masses], buffer_size);

  /* Add the smoothing lengths. */
  *mask |= logger_add_field_to_mask(masks[hydro_logger_field_smoothing_lengths],
                                    buffer_size);

  /* Add the entropies. */
  *mask |= logger_add_field_to_mask(masks[hydro_logger_field_entropies],
                                    buffer_size);

  /* Add the ID. */
  *mask |= logger_add_field_to_mask(masks[hydro_logger_field_particle_ids],
                                    buffer_size);

  /* Add the density. */
  *mask |= logger_add_field_to_mask(masks[hydro_logger_field_densities],
                                    buffer_size);
}

/**
 * @brief Write a particle to the logger.
 *
 * WARNING: The order should be the same in all the functions and
 * #hydro_logger_fields!
 *
 * @param masks The list of masks (same order than in #hydro_logger_init).
 * @param p The #part to write.
 * @param xp The #xpart to write.
 * @param mask The mask to use for this record.
 * @param buff The buffer where to write the particle.
 *
 * @return The buffer after the data.
 */
INLINE static char *hydro_logger_write_particle(
    const struct mask_data *mask_data, const struct part *p,
    const struct xpart *xp, unsigned int *mask, char *buff) {

  /* Write the coordinate. */
  if (logger_should_write_field(mask_data[hydro_logger_field_coordinates],
                                mask)) {
    memcpy(buff, p->x, 3 * sizeof(double));
    buff += 3 * sizeof(double);
  }

  /* Write the velocity. */
  if (logger_should_write_field(mask_data[hydro_logger_field_velocities],
                                mask)) {
    memcpy(buff, p->v, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Write the acceleration. */
  if (logger_should_write_field(mask_data[hydro_logger_field_accelerations],
                                mask)) {

    /* Compute the acceleration due to hydro and gravity */
    float *acc = (float *)buff;
    acc[0] = p->a_hydro[0];
    acc[1] = p->a_hydro[1];
    acc[2] = p->a_hydro[2];
    if (p->gpart) {
      acc[0] += p->gpart->a_grav[0];
      acc[1] += p->gpart->a_grav[1];
      acc[2] += p->gpart->a_grav[2];
    }

    memcpy(buff, acc, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Write the mass. */
  if (logger_should_write_field(mask_data[hydro_logger_field_masses], mask)) {
    memcpy(buff, &p->mass, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the smoothing length. */
  if (logger_should_write_field(mask_data[hydro_logger_field_smoothing_lengths],
                                mask)) {
    memcpy(buff, &p->h, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the entropy. */
  if (logger_should_write_field(mask_data[hydro_logger_field_entropies],
                                mask)) {
    memcpy(buff, &p->entropy, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the Id. */
  if (logger_should_write_field(mask_data[hydro_logger_field_particle_ids],
                                mask)) {
    memcpy(buff, &p->id, sizeof(long long));
    buff += sizeof(long long);
  }

  /* Write the density. */
  if (logger_should_write_field(mask_data[hydro_logger_field_densities],
                                mask)) {
    memcpy(buff, &p->rho, sizeof(float));
    buff += sizeof(float);
  }

  return buff;
}

#endif  // WITH_LOGGER
#endif  // SWIFT_GADGET2_HYDRO_LOGGER_H
