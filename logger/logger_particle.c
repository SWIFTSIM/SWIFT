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
#include "logger_particle.h"
#include "logger_header.h"
#include "logger_loader_io.h"
#include "logger_reader.h"
#include "logger_time.h"
#include "logger_tools.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Print the properties of a logger_particle.
 *
 * @param p The #logger_particle to print
 */
void logger_particle_print(const struct logger_particle *p) {
  message("ID:            %lli.", p->id);
  message("Mass:          %g", p->mass);
  message("Time:          %g.", p->time);
  message("Cutoff Radius: %g.", p->h);
  message("Positions:     (%g, %g, %g).", p->pos[0], p->pos[1], p->pos[2]);
  message("Velocities:    (%g, %g, %g).", p->vel[0], p->vel[1], p->vel[2]);
  message("Accelerations: (%g, %g, %g).", p->acc[0], p->acc[1], p->acc[2]);
  message("Entropy:       %g.", p->entropy);
  message("Density:       %g.", p->density);
  message("Type:          %i.", p->type);
}

/**
 * @brief Initialize a logger_particle.
 *
 * @param part The #logger_particle to initialize.
 */
void logger_particle_init(struct logger_particle *part) {
  for (size_t k = 0; k < DIM; k++) {
    part->pos[k] = 0;
    part->vel[k] = 0;
    part->acc[k] = 0;
  }

  part->entropy = -1;
  part->density = -1;
  part->h = -1;
  part->mass = -1;
  part->id = SIZE_MAX;

  part->type = 0;
}

/**
 * @brief Read a single named entry for a particle.
 *
 * @param part The #logger_particle to update.
 * @param map The mapped data.
 * @param field field to read.
 * @param size number of bits to read.
 *
 * @return mapped data after the block read.
 */
void *logger_particle_read_field(struct logger_particle *part, void *map,
                                 const char *field, const size_t size) {
  void *p = NULL;

  /* Get the correct pointer. */
  if (strcmp("positions", field) == 0) {
    p = &part->pos;
  } else if (strcmp("velocities", field) == 0) {
    p = &part->vel;
  } else if (strcmp("accelerations", field) == 0) {
    p = &part->acc;
  } else if (strcmp("entropy", field) == 0) {
    p = &part->entropy;
  } else if (strcmp("smoothing length", field) == 0) {
    p = &part->h;
  } else if (strcmp("density", field) == 0) {
    p = &part->density;
  } else if (strcmp("consts", field) == 0) {
    p = malloc(size);
  } else if (strcmp("special flags", field) == 0) {
    p = &part->type;
  } else {
    error("Type %s not defined.", field);
  }

  /* read the data. */
  map = logger_loader_io_read_data(map, size, p);

  /* Split the required fields. */
  if (strcmp("consts", field) == 0) {
    part->mass = 0;
    part->id = 0;
    memcpy(&part->mass, p, sizeof(float));
    p = (char *)p + sizeof(float);
    memcpy(&part->id, p, sizeof(size_t));
    p = (char *)p - sizeof(float);
    free(p);
  }

  return map;
}

/**
 * @brief Read a particle entry in the log file.
 *
 * @param reader The #logger_reader.
 * @param part The #logger_particle to update.
 * @param offset offset of the record to read.
 * @param time time to interpolate (not used if constant interpolation).
 * @param reader_type #logger_reader_type.
 *
 * @return position after the record.
 */
size_t logger_particle_read(struct logger_particle *part,
                            const struct logger_reader *reader, size_t offset,
                            const double time,
                            const enum logger_reader_type reader_type) {

  /* Save the offset */
  part->offset = offset;

  /* Get a few pointers. */
  const struct header *h = &reader->log.header;
  void *map = reader->log.log.map;

  const struct time_array *times = &reader->log.times;

  size_t mask = 0;
  size_t h_offset = 0;

  logger_particle_init(part);

  /* Read the record's mask. */
  map = logger_loader_io_read_mask(h, (char *)map + offset, &mask, &h_offset);

  /* Check that the mask is meaningful */
  if (mask > (unsigned int)(1 << h->masks_count)) {
    error("Found an unexpected mask %zi", mask);
  }

  /* Check if it is not a time record. */
  if (mask == h->timestamp_mask) {
    error("Unexpected timestamp while reading a particle: %lu.", mask);
  }

  /* Read all the fields. */
  for (size_t i = 0; i < h->masks_count; i++) {
    if (mask & h->masks[i].mask) {
      map = logger_particle_read_field(part, map, h->masks[i].name,
                                       h->masks[i].size);
    }
  }

  /* Get the time of current record.
     This check is required for the manipulating the file before
     the initialization of the time_array. */
  if (times->size != 0) {
    part->time = time_array_get_time(times, offset);
  } else
    part->time = -1;

  /* update the offset. */
  offset = (size_t)((char *)map - (char *)reader->log.log.map);

  /* Check if an interpolation is required. */
  if (reader_type == logger_reader_const) return offset;

  /* Start reading next record. */
  struct logger_particle part_next;

  /* Check that the offset are in the correct direction. */
  if (!header_is_forward(h)) {
    error("Cannot read a particle with non forward offsets.");
  }

  /* No next particle. */
  if (h_offset == 0) return (size_t)((char *)map - (char *)reader->log.log.map);

  /* get absolute offset of next particle. */
  h_offset += offset - header_get_record_size_from_mask(h, mask) -
              LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE;

  /* Get time of next record. */
  part_next.time = time_array_get_time(times, h_offset);

  /* Read next record. */
  h_offset = logger_particle_read(&part_next, reader, h_offset, part_next.time,
                                  logger_reader_const);

  /* Interpolate the two particles. */
  logger_particle_interpolate(part, &part_next, time);

  return offset;
}

/**
 * @brief Compute the quintic hermite spline interpolation.
 *
 * @param t0 The time at the left of the interval.
 * @param x0 The function at the left of the interval.
 * @param v0 The first derivative at the left of the interval.
 * @param a0 The second derivative at the left of the interval.
 * @param t1 The time at the right of the interval.
 * @param x1 The function at the right of the interval.
 * @param v1 The first derivative at the right of the interval.
 * @param a1 The second derivative at the right of the interval.
 * @param t The time of the interpolation.
 *
 * @return The function evaluated at t.
 */
double logger_particle_quintic_hermite_spline(double t0, double x0, float v0,
                                              float a0, double t1, double x1,
                                              float v1, float a1, double t) {

  /* Generates recurring variables  */
  /* Time differences */
  const double dt = t1 - t0;
  const double dt2 = dt * dt;
  const double dt3 = dt2 * dt;
  const double dt4 = dt3 * dt;
  const double dt5 = dt4 * dt;

  const double t_t0 = t - t0;
  const double t_t0_2 = t_t0 * t_t0;
  const double t_t0_3 = t_t0_2 * t_t0;
  const double t_t1 = t - t1;
  const double t_t1_2 = t_t1 * t_t1;

  /* Derivatives */
  const double v0_dt = v0 * dt;
  const double a0_dt2 = 0.5 * a0 * dt2;
  const double v1_dt = v1 * dt;
  const double a1_dt2 = 0.5 * a1 * dt2;

  /* Do the first 3 terms of the hermite spline */
  double x = x0 + v0 * t_t0 + 0.5 * a0 * t_t0_2;

  /* Cubic term */
  x += (x1 - x0 - v0_dt - a0_dt2) * t_t0_3 / dt3;

  /* Quartic term */
  x += (3. * x0 - 3. * x1 + v1_dt + 2. * v0_dt + a0_dt2) * t_t0_3 * t_t1 / dt4;

  /* Quintic term */
  x += (6. * x1 - 6. * x0 - 3. * v0_dt - 3. * v1_dt + a1_dt2 - a0_dt2) *
       t_t0_3 * t_t1_2 / dt5;

  return x;
}

/**
 * @brief Compute the cubic hermite spline interpolation.
 *
 * @param t0 The time at the left of the interval.
 * @param v0 The first derivative at the left of the interval.
 * @param a0 The second derivative at the left of the interval.
 * @param t1 The time at the right of the interval.
 * @param v1 The first derivative at the right of the interval.
 * @param a1 The second derivative at the right of the interval.
 * @param t The time of the interpolation.
 *
 * @return The function evaluated at t.
 */
float logger_particle_cubic_hermite_spline(double t0, float v0, float a0,
                                           double t1, float v1, float a1,
                                           double t) {

  /* Generates recurring variables  */
  /* Time differences */
  const float dt = t1 - t0;
  const float dt2 = dt * dt;
  const float dt3 = dt2 * dt;

  const float t_t0 = t - t0;
  const float t_t0_2 = t_t0 * t_t0;
  const float t_t1 = t - t1;

  /* Derivatives */
  const float a0_dt = a0 * dt;
  const float a1_dt = a1 * dt;

  /* Do the first 2 terms of the hermite spline */
  float x = v0 + a0 * t_t0;

  /* Square term */
  x += (v1 - v0 - a0_dt) * t_t0_2 / dt2;

  /* Cubic term */
  x += (2. * v0 - 2. * v1 + a1_dt + a0_dt) * t_t0_2 * t_t1 / dt3;

  return x;
}

/**
 * @brief interpolate two particles at a given time
 *
 * @param part_curr #logger_particle In: current particle (before time), Out:
 * interpolated particle
 * @param part_next #logger_particle next particle (after time)
 * @param time interpolation time
 *
 */
void logger_particle_interpolate(struct logger_particle *part_curr,
                                 const struct logger_particle *part_next,
                                 const double time) {

  /* Check that a particle is provided. */
  if (!part_curr) error("part_curr is NULL.");
  if (!part_next) error("part_next is NULL.");

#ifdef SWIFT_DEBUG_CHECKS
  /* Check the particle order. */
  if (part_next->time < part_curr->time)
    error("Wrong particle order (next before current): %g, %g", part_next->time,
          part_curr->time);
  if ((time < part_curr->time) || (part_next->time < time))
    error(
        "Cannot extrapolate (particle time: %f, "
        "interpolating time: %f, next particle time: %f).",
        part_curr->time, time, part_next->time);
#endif

  /* Compute the interpolation scaling. */
  double scaling = part_next->time - part_curr->time;

  scaling = (time - part_curr->time) / scaling;

  float ftmp;

  /* interpolate vectors. */
  for (int i = 0; i < 3; i++) {
    /* Positions */
    part_curr->pos[i] = logger_particle_quintic_hermite_spline(
        part_curr->time, part_curr->pos[i], part_curr->vel[i],
        part_curr->acc[i], part_next->time, part_next->pos[i],
        part_next->vel[i], part_next->acc[i], time);

    /* Velocities */
    part_curr->vel[i] = logger_particle_cubic_hermite_spline(
        part_curr->time, part_curr->vel[i], part_curr->acc[i], part_next->time,
        part_next->vel[i], part_next->acc[i], time);

    /* Accelerations */
    ftmp = (part_next->acc[i] - part_curr->acc[i]);
    part_curr->acc[i] += ftmp * scaling;
  }

  /* interpolate scalars. */
  ftmp = (part_next->entropy - part_curr->entropy);
  part_curr->entropy += ftmp * scaling;

  ftmp = (part_next->h - part_curr->h);
  part_curr->h += ftmp * scaling;

  ftmp = (part_next->density - part_curr->density);
  part_curr->density += ftmp * scaling;

  /* set time. */
  part_curr->time = time;
}
