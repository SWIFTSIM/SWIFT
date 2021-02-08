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
 * @brief Compute the number of fields for a given particle type.
 *
 * @param type The type of particles.
 *
 * @return The number of fields.
 */
int tools_get_number_fields(enum part_type type) {
  int number_fields = 0;
  switch (type) {
    case swift_type_gas:
      number_fields = hydro_logger_field_count;
      number_fields += chemistry_logger_field_part_count;
      return number_fields;

    case swift_type_dark_matter:
    case swift_type_dark_matter_background:
      number_fields = gravity_logger_field_count;
      return number_fields;

    case swift_type_stars:
      number_fields = stars_logger_field_count;
      number_fields += chemistry_logger_field_spart_count;
      number_fields += star_formation_logger_field_count;
      return number_fields;

    default:
      message("Particle type %i not implemented, skipping it", type);
      return 0;
  }
}

#define copy_field_to_struct_internal(MODULE, PART, TYPE)             \
  for (int j = 0; j < MODULE##_logger_field##PART##_count; j++) {     \
                                                                      \
    /* Save the main properties */                                    \
    fields[i].module = TYPE;                                          \
    fields[i].name = MODULE##_logger_field_names##PART[j];            \
                                                                      \
    /* Get the indexes */                                             \
    const int global = MODULE##_logger_local_to_global##PART[j];      \
    int first = h->masks[global].reader.first_deriv;                  \
    int second = h->masks[global].reader.second_deriv;                \
                                                                      \
    /* Save the global indexes */                                     \
    fields[i].global_index = global;                                  \
    fields[i].global_index_first = first;                             \
    fields[i].global_index_second = second;                           \
                                                                      \
    /* Convert the first derivatives into local index */              \
    if (first != -1) {                                                \
      for (int k = 0; k < MODULE##_logger_field##PART##_count; k++) { \
        if (MODULE##_logger_local_to_global##PART[k] == first) {      \
          first = k;                                                  \
          break;                                                      \
        }                                                             \
      }                                                               \
    }                                                                 \
                                                                      \
    /* Convert the second derivatives into local index */             \
    if (second != -1) {                                               \
      for (int k = 0; k < MODULE##_logger_field##PART##_count; k++) { \
        if (MODULE##_logger_local_to_global##PART[k] == second) {     \
          second = k;                                                 \
          break;                                                      \
        }                                                             \
      }                                                               \
    }                                                                 \
                                                                      \
    /* Initialize the structure */                                    \
    fields[i].local_index = j;                                        \
    fields[i].local_index_first = first;                              \
    fields[i].local_index_second = second;                            \
    i++;                                                              \
  }

/**
 * Same function as set_links_local_global_internal before but with only two
 * arguments.
 */
#define copy_field_to_struct_single_particle_type(MODULE, TYPE) \
  copy_field_to_struct_internal(MODULE, , TYPE)

/**
 * Same function as set_links_local_global_internal before but with a cleaner
 * argument.
 */
#define copy_field_to_struct(MODULE, PART, TYPE) \
  copy_field_to_struct_internal(MODULE, _##PART, TYPE)

/**
 * @brief Construct the list of fields for a given particle type.
 *
 * @param fields (output) The list of fields (need to be already allocated).
 * @param type The type of particle.
 * @param h The #header
 */
void tools_get_list_fields(struct field_information *fields,
                           enum part_type type, const struct header *h) {
  int i = 0;
  switch (type) {
    case swift_type_gas:
      copy_field_to_struct_single_particle_type(hydro, field_module_default);
      copy_field_to_struct(chemistry, part, field_module_chemistry);
      break;

    case swift_type_dark_matter:
    case swift_type_dark_matter_background:
      copy_field_to_struct_single_particle_type(gravity, field_module_default);
      break;

    case swift_type_stars:
      copy_field_to_struct_single_particle_type(stars, field_module_default);
      copy_field_to_struct(chemistry, spart, field_module_chemistry);
      copy_field_to_struct_single_particle_type(star_formation,
                                                field_module_star_formation);
      break;

    default:
      error("Particle type %i not implemented", type);
  }
}

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
int tools_get_next_record(const struct header *h, char *map, size_t *offset,
                          size_t file_size) {
  if (header_is_forward(h))
    return _tools_get_next_record_forward(h, map, offset);
  if (header_is_backward(h))
    return _tools_get_next_record_backward(h, map, offset, file_size);
  else
    error_python("Offsets are corrupted.");
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
int _tools_get_next_record_forward(const struct header *h, char *map,
                                   size_t *offset) {
  size_t diff_offset = 0;

  /* Read the offset. */
  map = logger_loader_io_read_mask(h, map + *offset, NULL, &diff_offset);

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
int _tools_get_next_record_backward(const struct header *h, char *map,
                                    size_t *offset, size_t file_size) {
#ifndef SWIFT_DEBUG_CHECKS
  error_python("Should not be used, method too slow");
#endif
  size_t current_offset = *offset;
  size_t record_header = LOGGER_MASK_SIZE + LOGGER_OFFSET_SIZE;

  while (current_offset < file_size) {
    size_t mask = 0;
    size_t prev_offset;
    logger_loader_io_read_mask(h, map + current_offset, &mask, &prev_offset);

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
size_t tools_reverse_offset(const struct header *h, char *file_map,
                            size_t offset) {
  size_t mask = 0;
  size_t prev_offset = 0;
  const size_t cur_offset = offset;
  char *map = file_map;

  /* read mask + offset. */
  map = logger_loader_io_read_mask(h, map + offset, &mask, &prev_offset);

  /* write offset of zero (in case it is the last record). */
  const size_t zero = 0;
  map = map - LOGGER_OFFSET_SIZE;
  map = logger_loader_io_write_data(map, LOGGER_OFFSET_SIZE, &zero);

  /* set offset after current record. */
  map = map + header_get_record_size_from_mask(h, mask);
  size_t after_current_record = (size_t)(map - file_map);

  /* first records do not have a previous partner. */
  if (prev_offset == cur_offset) return after_current_record;
  if (prev_offset > cur_offset)
    error_python("Unexpected offset: header %lu, current %lu.", prev_offset,
                 cur_offset);

  /* modify previous offset. */
  map = file_map + cur_offset - prev_offset + LOGGER_MASK_SIZE;
  map = logger_loader_io_write_data(map, LOGGER_OFFSET_SIZE, &prev_offset);

#ifdef SWIFT_DEBUG_CHECKS
  size_t prev_mask = 0;
  map = map - LOGGER_MASK_SIZE - LOGGER_OFFSET_SIZE;
  logger_loader_io_read_mask(h, map, &prev_mask, NULL);

  /* Check if we are not mixing timestamp and particles */
  if ((prev_mask != h->timestamp_mask && mask == h->timestamp_mask) ||
      (prev_mask == h->timestamp_mask && mask != h->timestamp_mask))
    error_python("Unexpected mask: %lu, got %lu.", mask, prev_mask);

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
  error_python("Should not check in non debug mode.");
#endif

  const struct header *h = &reader->log.header;
  char *file_init = reader->log.log.map;
  char *map = file_init + offset;

  size_t mask;
  size_t pointed_offset;

  const size_t mask_special_flag =
      h->masks[header_get_field_index(h, "SpecialFlags")].mask;

  /* read mask + offset. */
  map = logger_loader_io_read_mask(h, map, &mask, &pointed_offset);

  /* set offset after current record. */
  map = map + header_get_record_size_from_mask(h, mask);
  const size_t offset_ret = (size_t)(map - file_init);

  /* If something happened, skip the check. */
  if (mask & mask_special_flag) {
    return offset_ret;
  }

  /* get absolute offset. */
  if (header_is_forward(h))
    pointed_offset += offset;
  else if (header_is_backward(h)) {
    if (offset < pointed_offset)
      error_python("Offset too large (%lu) at %lu with mask %lu.",
                   pointed_offset, offset, mask);
    pointed_offset = offset - pointed_offset;
  } else {
    error_python("Offset are corrupted.");
  }

  if (pointed_offset == offset || pointed_offset == 0) return offset_ret;

  /* read mask of the pointed record. */
  size_t pointed_mask = 0;
  logger_loader_io_read_mask(h, file_init + pointed_offset, &pointed_mask,
                             NULL);

  /* check if not mixing timestamp and particles. */
  if ((pointed_mask != h->timestamp_mask && mask == h->timestamp_mask) ||
      (pointed_mask == h->timestamp_mask && mask != h->timestamp_mask))
    error_python("Error in the offset (mask %lu at %lu != %lu at %lu).", mask,
                 offset, pointed_mask, pointed_offset);

  return offset_ret;
}

/**
 * @brief Remove the previous text and print the progress of current task.
 *
 * @param percentage The current progress.
 * @param remaining_time The remaining time of the process (in seconds)
 * @param message The message to display before the progress message.
 */
void tools_print_progress(float percentage, int remaining_time,
                          const char *message) {

  /* Compute the time */
  int hour = remaining_time / (60 * 60);
  remaining_time -= hour * 60 * 60;
  int min = remaining_time / 60;
  int sec = remaining_time - min * 60;

  /* Allocate the output string */
  char output[300];

  /* Write the message */
  char *current =
      output + sprintf(output, "%s: %2.1f%% done, Remaining time: ", message,
                       percentage);

  /* Write the hour */
  if (hour == 0)
    current += sprintf(current, "    ");
  else
    current += sprintf(current, "%.2ih ", hour);

  /* Write the minutes */
  if (hour == 0 && min == 0)
    current += sprintf(current, "      ");
  else
    current += sprintf(current, "%.2imin ", min);

  /* Write the seconds */
  current += sprintf(current, "%.2is", sec);

  /* Print the string */
  printf("\r%s", output);
  fflush(stdout);
}
