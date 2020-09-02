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
#include "logger_tools.h"

#include "logger_header.h"
#include "logger_loader_io.h"
#include "logger_particle.h"
#include "logger_reader.h"

#include <stdio.h>

/**
 * @brief get the offset of the next corresponding record.
 *
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset In: initial offset, Out: offset of the next record
 * @param file_size The file size.
 *
 * @return -1 if no next record, otherwise 0
 */
int tools_get_next_record(const struct header *h, void *map, size_t *offset,
                          size_t file_size) {
  if (header_is_forward(h))
    return _tools_get_next_record_forward(h, map, offset);
  if (header_is_backward(h))
    return _tools_get_next_record_backward(h, map, offset, file_size);
  else
    error("Offsets are corrupted.");
}

/**
 * @brief internal function of #tools_get_next_record. Should not be used
 * outside.
 *
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset (Out) offset of the next record
 *
 * @return error code, -1 if no next record
 */
int _tools_get_next_record_forward(const struct header *h, void *map,
                                   size_t *offset) {
  size_t diff_offset = 0;

  /* Read the offset. */
  map =
      logger_loader_io_read_mask(h, (char *)map + *offset, NULL, &diff_offset);

  if (diff_offset == 0) return -1;

  /* Set the absolute offset. */
  *offset += diff_offset;
  return 0;
}

/**
 * @brief internal function of #tools_get_next_record. Should not be used (very
 * slow)
 *
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset In: initial offset, Out: offset of the next record
 * @param file_size The file size.
 *
 * @return error code, -1 if no next record
 */
int _tools_get_next_record_backward(const struct header *h, void *map,
                                    size_t *offset, size_t file_size) {
#ifndef SWIFT_DEBUG_CHECKS
  error("Should not be used, method too slow");
#endif
  size_t current_offset = *offset;
  size_t record_header = LOGGER_MASK_SIZE + LOGGER_OFFSET_SIZE;

  while (current_offset < file_size) {
    size_t mask = 0;
    size_t prev_offset;
    logger_loader_io_read_mask(h, (char *)map + current_offset, &mask,
                               &prev_offset);

    prev_offset = current_offset - prev_offset - record_header;
    if (*offset == prev_offset) {
      *offset = current_offset - record_header;
      return 0;
    }

    current_offset += header_get_record_size_from_mask(h, mask);
  }

  return -1;
}

/**
 * @brief switch side offset.
 *
 * From current record, switch side of the offset of the previous one.
 * @param h #header structure of the file.
 * @param file_map file mapping.
 * @param offset position of the record.
 *
 * @return position after the record.
 */
size_t tools_reverse_offset(const struct header *h, void *file_map,
                            size_t offset) {
  size_t mask = 0;
  size_t prev_offset = 0;
  const size_t cur_offset = offset;
  void *map = file_map;

  /* read mask + offset. */
  map =
      logger_loader_io_read_mask(h, (char *)map + offset, &mask, &prev_offset);

  /* write offset of zero (in case it is the last record). */
  const size_t zero = 0;
  map = (char *)map - LOGGER_OFFSET_SIZE;
  map = logger_loader_io_write_data(map, LOGGER_OFFSET_SIZE, &zero);

  /* set offset after current record. */
  map = (char *)map + header_get_record_size_from_mask(h, mask);
  size_t after_current_record = (size_t)((char *)map - (char *)file_map);

  /* first records do not have a previous partner. */
  if (prev_offset == cur_offset) return after_current_record;

  if (prev_offset > cur_offset)
    error("Unexpected offset: header %lu, current %lu.", prev_offset,
          cur_offset);

  /* modify previous offset. */
  map = (char *)file_map + cur_offset - prev_offset + LOGGER_MASK_SIZE;
  map = logger_loader_io_write_data(map, LOGGER_OFFSET_SIZE, &prev_offset);

#ifdef SWIFT_DEBUG_CHECKS
  size_t prev_mask = 0;
  map = (char *)map - LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE;
  logger_loader_io_read_mask(h, map, &prev_mask, NULL);

  /* Check if we are not mixing timestamp and particles */
  if ((prev_mask != h->timestamp_mask && mask == h->timestamp_mask) ||
      (prev_mask == h->timestamp_mask && mask != h->timestamp_mask))
    error("Unexpected mask: %lu, got %lu.", mask, prev_mask);

#endif  // SWIFT_DEBUG_CHECKS

  return after_current_record;
}

/**
 * @brief debugging function checking the offset and the mask of a record.
 *
 * Compare the mask with the one pointed by the header.
 * if the record is a particle, check the id too.
 *
 * @param reader The #logger_reader.
 * @param offset position of the record.
 *
 * @return position after the record.
 */
size_t tools_check_record_consistency(const struct logger_reader *reader,
                                      size_t offset) {
#ifndef SWIFT_DEBUG_CHECKS
  error("Should not check in non debug mode.");
#endif

  const struct header *h = &reader->log.header;
  void *file_init = reader->log.log.map;
  void *map = (char *)file_init + offset;

  size_t mask;
  size_t pointed_offset;

  /* read mask + offset. */
  map = logger_loader_io_read_mask(h, map, &mask, &pointed_offset);

  /* get absolute offset. */
  if (header_is_forward(h))
    pointed_offset += offset;
  else if (header_is_backward(h)) {
    if (offset < pointed_offset)
      error("Offset too large (%lu) at %lu with mask %lu.", pointed_offset,
            offset, mask);
    pointed_offset = offset - pointed_offset;
  } else {
    error("Offset are corrupted.");
  }

  /* set offset after current record. */
  map = (char *)map + header_get_record_size_from_mask(h, mask);

  if (pointed_offset == offset || pointed_offset == 0)
    return (size_t)((char *)map - (char *)file_init);

  /* read mask of the pointed record. */
  size_t pointed_mask = 0;
  logger_loader_io_read_mask(h, (char *)file_init + pointed_offset,
                             &pointed_mask, NULL);

  /* check if not mixing timestamp and particles. */
  if ((pointed_mask != h->timestamp_mask && mask == h->timestamp_mask) ||
      (pointed_mask == h->timestamp_mask && mask != h->timestamp_mask))
    error("Error in the offset (mask %lu at %lu != %lu at %lu).", mask, offset,
          pointed_mask, pointed_offset);

  return (size_t)((char *)map - (char *)file_init);
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
double logger_tools_quintic_hermite_spline(double t0, double x0, float v0,
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
float logger_tools_cubic_hermite_spline(double t0, float v0, float a0,
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
