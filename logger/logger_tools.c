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
#include "logger_io.h"
#include "logger_reader.h"

#include "logger_particle.h"

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
    error("Offsets are corrupted");
}

/**
 * @brief internal function of #tools_get_next_record. Should not be used outside.
 *
 * @param h #header structure of the file
 * @param map file mapping
 * @param offset In: initial offset, Out: offset of the next record
 *
 * @return error code, -1 if no next record
 */
int _tools_get_next_record_forward(const struct header *h, void *map,
                                  size_t *offset) {
  size_t diff_offset = 0;

  /* read offset */
  *offset = io_read_mask(h, map, *offset, NULL, &diff_offset);

  if (diff_offset == 0) return -1;

  /* set absolute offset */
  *offset += diff_offset - LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE;
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
    current_offset = io_read_mask(h, map, current_offset, &mask, &prev_offset);

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
 * @param map file mapping.
 * @param offset position of the record.
 *
 * @return position after the record.
 */
size_t tools_reverse_offset(const struct header *h, void *map, size_t offset) {
  size_t mask = 0;
  size_t prev_offset = 0;
  const size_t cur_offset = offset;

  /* read mask + offset */
  offset = io_read_mask(h, map, offset, &mask, &prev_offset);

  /* write offset of zero (in case it is the last record) */
  const size_t zero = 0;
  offset -= LOGGER_OFFSET_SIZE;
  offset = io_write_data(map, LOGGER_OFFSET_SIZE, &zero, offset);

  /* set offset after current record */
  offset += header_get_record_size_from_mask(h, mask);

  /* first records do not have a previous partner */
  if (prev_offset == cur_offset)
    return offset;

  if (prev_offset > cur_offset)
    error("Unexpected offset, header %lu, current %lu", prev_offset,
          cur_offset);

  /* modify previous offset */
  size_t abs_prev_offset = cur_offset - prev_offset + LOGGER_MASK_SIZE;
  abs_prev_offset = io_write_data(map, LOGGER_OFFSET_SIZE, &prev_offset, abs_prev_offset);

#ifdef SWIFT_DEBUG_CHECKS
  size_t prev_mask = 0;
  abs_prev_offset -= LOGGER_MASK_SIZE + LOGGER_OFFSET_SIZE;
  abs_prev_offset = io_read_mask(h, map, abs_prev_offset, &prev_mask, NULL);

  if (prev_mask != mask)
    error("Unexpected mask: %lu (at %lu), got %lu (at %lu)", mask, offset,
          prev_mask, prev_offset);

#endif  // SWIFT_DEBUG_CHECKS

  return offset;
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
size_t tools_check_record_consistency(const struct logger_reader *reader, size_t offset) {
#ifndef SWIFT_DEBUG_CHECKS
  error("Should not check in non debug mode");
#endif

  const struct header *h = &reader->log.header;
  void *map = reader->log.log.map;

  size_t tmp = offset;

  size_t mask;
  size_t pointed_offset;

  /* read mask + offset */
  offset = io_read_mask(h, map, offset, &mask, &pointed_offset);

  /* get absolute offset */
  if (header_is_forward(h))
    pointed_offset += tmp;
  else if (header_is_backward(h)) {
    if (tmp < pointed_offset)
      error("Offset too large (%lu) at %lu with mask %lu", pointed_offset, tmp,
            mask);
    pointed_offset = tmp - pointed_offset;
  }
  else {
    error("Offset are corrupted");
  }

  /* set offset after current record */
  offset += header_get_record_size_from_mask(h, mask);

  if (pointed_offset == tmp || pointed_offset == 0)
    return offset;

  /* read mask of the pointed record */
  size_t pointed_mask = 0;
  pointed_offset = io_read_mask(h, map, pointed_offset, &pointed_mask, NULL);

  /* check masks */
  if (pointed_mask != mask)
    error("Error in the offset (mask %lu != %lu) at %lu and %lu", mask,
          pointed_mask,
          offset - header_get_record_size_from_mask(h, mask) - LOGGER_MASK_SIZE -
	  LOGGER_OFFSET_SIZE,
          pointed_offset - LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE);

  if (pointed_mask == 128)
    return offset;

  struct logger_particle part;
  tmp = logger_particle_read(&part, reader, tmp, 0, logger_reader_const);

  size_t id = part.id;
  tmp = pointed_offset - LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE;
  tmp = logger_particle_read(&part, reader, tmp, 0, logger_reader_const);

  if (id != part.id)
    error("Offset wrong, id incorrect (%lu != %lu) at %lu", id, part.id, tmp);

  return offset;
}

